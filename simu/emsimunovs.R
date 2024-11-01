M_step_novs = function(weights,
                       old_sigmasq,
                       inputlist,
                       K,
                       lambda = NULL,
                       init = FALSE,
                       lambdamax = 1e10,
                       plot = F,
                       maxiter = NULL,
                       alpha2 = 0) {
  aa = Sys.time()
  
  N = inputlist$N
  Ts = inputlist$Ts
  Y = inputlist$Y
  X_i = inputlist$X_i
  
  p = inputlist$p
  Q = inputlist$Q
  Linv = inputlist$Linv
  nbeta_basis = inputlist$nbeta_basis
  Ts = length(Ts)
  
  
  
  Y_vec = rep(unlist(Y), K)
  Q_mat = kronecker(diag(1, K), cbind(matrix(1, N * 10, 1), do.call(rbind, Q)))
  
  
  prop = apply(weights, 2, mean)
  
  WLS_weights = NULL
  for (k in 1:K) {
    WLS_weights = c(WLS_weights, rep(weights[, k] / old_sigmasq[k], each = 10))
  }
  WLS_weights = WLS_weights / max(WLS_weights)
  
  
  
  fit = glmnet(
    x = Q_mat,
    y = Y_vec,
    family = "gaussian",
    weights = WLS_weights,
    alpha = alpha2,
    maxit = maxiter
  )
  
  df = matrix(0, length(fit$lambda), 1)
  aa = Sys.time()
  for (k in 1:K) {
    svdx = svd(Q_mat[(1 + N * Ts * (k - 1)):(N * Ts * k),
                     (1 + (nbeta_basis * p + 1) * (k - 1)):((nbeta_basis *
                                                               p + 1) * (k))]
               , nu = 0, nv = 0)
    d = svdx$d
    for (i in 1:length(fit$lambda)) {
      df[i] = df[i] + sum(d / (d + fit$lambda[i]))
    }
  }
  print(Sys.time() - aa)
  
  fit$df = df
  updateBIC = matrix(0, length(fit$lambda), 1)
  updateAIC = matrix(0, length(fit$lambda), 1)
  
  
  Ypre = cbind(rep(1, N * Ts * K), Q_mat) %*% rbind(fit$a0, fit$beta)
  
  for (lam in 1:length(fit$lambda)) {
    alpha = list()
    for (k in 1:K) {
      alpha[[k]] = fit$beta[(1 + (nbeta_basis * p + 1) * (k - 1)):((nbeta_basis *
                                                                      p + 1) * (k)), lam]
      alpha[[k]][1] = alpha[[k]][1] + fit$a0[lam]
    }
    
    termloss = array(((unlist(Y) - Ypre[, lam]) ^ 2), c(Ts, N, K))
    termloss_indi = t(apply(termloss, c(2, 3), mean))
    
    new_sigmasq = matrix(0, K, 1)
    for (k in 1:K) {
      new_sigmasq[k] = sum(termloss_indi[k, ] * weights[, k]) / sum(weights[, k])
    }
    
    termloss_divide = apply(termloss, c(1, 2), function(x)
      log(sum(prop *
                (
                  exp(-1 * as.vector(x) / (2 * new_sigmasq)) /
                    sqrt(new_sigmasq)
                ))))
    ll = sum(termloss_divide)
    
    updateBIC[lam] = -2 * ll + ((fit$df[lam]) * log(10 * N))
    updateAIC[lam] = -2 * ll + 2 * (fit$df[lam])
    
  }
  
  lastone = which.min(updateBIC)
  if (fit$lambda[lastone] > lambdamax) {
    lastone = which(fit$lambda < lambdamax)[1]
    if (is.na(lastone))
      lastone = length(fit$lambda)
  }
  if (init == TRUE)
    lastone = length(fit$lambda)
  cat("lambda choice: ", lastone, "of", length(fit$lambda), '\n')
  cat("lambda=",
      fit$lambda[lastone],
      " maxlambda=",
      lambdamax,
      "at",
      min(which(fit$lambda < lambdamax)),
      '\n')
  
  alpha = list()
  for (k in 1:K) {
    alpha[[k]] = fit$beta[(1 + (nbeta_basis * p + 1) * (k - 1)):((nbeta_basis *
                                                                    p + 1) * (k)), lastone]
    alpha[[k]][1] = alpha[[k]][1] + fit$a0[lastone]
  }
  
  b = alpha
  for (k in 1:K) {
    for (j in 1:p) {
      var_index = (2 + (j - 1) * nbeta_basis):(1 + (j) * nbeta_basis)
      b[[k]][var_index] = Linv %*% alpha[[k]][var_index]
    }
  }
  
  termloss = array(((unlist(Y) - Ypre[, lastone]) ^ 2), c(Ts, N, K))
  termloss_indi = t(apply(termloss, c(2, 3), mean))
  
  new_sigmasq = matrix(0, K, 1)
  for (k in 1:K) {
    new_sigmasq[k] = sum(termloss_indi[k, ] * weights[, k]) / sum(weights[, k])
  }
  
  Z = apply(weights, 1, which.max)
  
  Ypre_all = Ypre[, lastone]
  
  Rsq = matrix(0, K, K)
  for (k1 in 1:K) {
    kind = Z == k1
    Y_true = matrix(Y_vec[1:(N * 10)], 10, N)[, kind]
    for (k2 in 1:K) {
      Y_pre = matrix(Ypre_all[(1 + (k2 - 1) * N * 10):(k2 * N * 10)], 10, N)[, kind]
      RSS = mean((Y_true - Y_pre) ^ 2)
      ERR = mean((Y_true - mean(Y_true)) ^ 2)
      Rsq[k1, k2] = 1 - RSS / ERR
    }
  }
  
  termloss_divide = apply(termloss, c(1, 2), function(x)
    log(sum(prop *
              (
                exp(-1 * as.vector(x) / (2 * new_sigmasq)) /
                  sqrt(new_sigmasq)
              ))))
  ll = sum(termloss_divide)
  
  res = list(
    alpha = alpha,
    b = b,
    df = fit$df[lastone],
    sigmasq = new_sigmasq,
    Rsq = Rsq,
    BIC = updateBIC[lastone],
    AIC = updateAIC[lastone],
    prop = prop,
    weights = weights,
    ll = ll,
    fit = fit,
    bestone = lastone
  )
  return(res)
}
########################################################################################
EM_novs = function(inputlist,
                   weights_init = NULL,
                   sigmasq_init = NULL,
                   Rseed,
                   K,
                   task,
                   save = FALSE,
                   inflation = 2,
                   plot_M_step = TRUE,
                   plot_converge_res = TRUE,
                   max_inner_iter = 10000,
                   max_outer_iter = 100,
                   numbers) {
  N = inputlist$N
  Ts = inputlist$Ts
  Y = inputlist$Y
  X_i = inputlist$X_i
  
  p = inputlist$p
  Q = inputlist$Q
  Linv = inputlist$Linv
  nbeta_basis = inputlist$nbeta_basis
  
  
  if (is.null(weights_init)) {
    weights_init = matrix(0, N, K)
    init_label = sample(1:N, N, replace = FALSE)
    len = ceiling(N / K)
    for (k in 1:K)
      weights_init[init_label[(1 + len * (k - 1)):min((len * k), N)] , k] = 1
  }
  if (is.null(sigmasq_init)) {
    sigmasq_init = rep(1, K)
  }
  
  M_step_res = M_step_novs(
    weights = weights_init,
    old_sigmasq = sigmasq_init,
    inputlist = inputlist,
    K = K,
    lambda = NULL,
    init = TRUE,
    maxiter = max_inner_iter
  )#0.007148256
  lambda_seq = M_step_res$fit$lambda
  
  
  #######
  coefnorm_old = 0
  
  coefnorm_new = unlist(M_step_res$alpha)
  
  weights_old = matrix(0, N, K)
  weights_new = weights_init
  
  lambdaold = 10
  lambdanew = 100
  
  Z_old = rep(1, N)
  Z_new = apply(weights_new, 1, which.max)
  loopind = 1
  max_iter = max_outer_iter
  Zpath = matrix(0, max_iter, N)
  
  while (loopind < max_iter &
         ((sum((
           coefnorm_old - coefnorm_new
         ) ^ 2) /
         (sum((
           coefnorm_new
         ) ^ 2) + sum((
           coefnorm_old
         ) ^ 2))) > 1e-5  |
         abs(lambdaold - lambdanew) / (lambdanew + lambdaold) > 1e-2)) {
    cat("iter", loopind, " ")
    for (k in 1:K)
      cat(mean(weights_old[, k]), " ")
    cat("\n")
    Z_old = apply(weights_old, 1, which.max)
    
    cat("Zdiff: ", sum(Z_old != Z_new))
    cat(" ")
    
    cat("coef diff:", (sum((
      coefnorm_old - coefnorm_new
    ) ^ 2) /
      (sum(
        (coefnorm_new) ^ 2
      ) + sum(
        (coefnorm_old) ^ 2
      ))))
    cat("lambdadiff: ",
        abs(lambdaold - lambdanew) / (lambdanew + lambdaold),
        '\n')
    
    coefnorm_old = coefnorm_new
    weights_old = weights_new
    lambdaold = lambdanew
    
    
    weights_new = E_step_simu(inputlist = inputlist,
                              M_step_res = M_step_res,
                              K = K)
    Z_new = apply(weights_new, 1, which.max)
    Zpath[loopind, ] = Z_new
    M_step_res = M_step_novs(
      weights = weights_new,
      old_sigmasq = M_step_res$sigmasq,
      inputlist = inputlist,
      K = K,
      maxiter = max_inner_iter,
      lambdamax = inflation * M_step_res$fit$lambda[M_step_res$bestone],
      lambda = NULL
    )
    
    coefnorm_new = unlist(M_step_res$alpha)
    lambdanew = M_step_res$fit$lambda[M_step_res$bestone]
    
    cat(" \n")
    loopind = loopind + 1
  }
  diverge = FALSE
  if (loopind >= max_iter)
    diverge = TRUE
  
  
  Z_new = apply(weights_new, 1, which.max)
  BIC = M_step_res$BIC
  AIC = M_step_res$AIC
  
  cat(BIC)
  cat("\n")
  
  Zind = matrix(0, 3, 1)
  Zind[1] = as.numeric(names(which.max(table(Z_new[1:numbers[1]]))))
  Zind[2] = as.numeric(names(which.max(table(Z_new[(1 + numbers[1]):(numbers[1] +
                                                                       numbers[2])]))))
  Zind[3] = as.numeric(names(which.max(table(Z_new[(1 + numbers[1] + numbers[2]):(numbers[1] +
                                                                                    numbers[2] + numbers[3])]))))
  if (K == 3 & length(table(Z_new)) == 3 & length(table(Zind)) == 3) {
    weightstemp = weights_new
    
    sigmasqtemp = M_step_res$sigmasq
    for (j in 1:3) {
      weights_new[, j] = weightstemp[, Zind[j]]
      M_step_res$sigmasq[j] = sigmasqtemp[Zind[j]]
    }
    
    
    M_step_res = M_step_novs(
      weights = weights_new,
      old_sigmasq = M_step_res$sigmasq,
      inputlist = inputlist,
      K = K,
      lambda = NULL,
      maxiter = max_inner_iter
    )
    
    Z_new = apply(weights_new, 1, which.max)
    
    
    
  }
  
  
  
  
  
  beta_est = list()
  eval.beta = eval.basis(seq(0, 1, , 100), inputlist$beta_basis)
  for (j in 1:p) {
    beta_est[[j]] = matrix(0, 100, K)
    for (k in 1:K) {
      beta_est[[j]][, k] =
        M_step_res$b[[k]][(2 + (j - 1) * nbeta_basis):(1 + j * nbeta_basis)] %*%
        t(eval.beta)
    }
  }
  
  res = list(
    beta_est = beta_est,
    updateBIC = BIC,
    updateAIC = AIC,
    Z = Z_new,
    K = K,
    intercept = sapply(M_step_res$alpha, function(x)
      x[1]),
    Rsq = M_step_res$Rsq,
    loopind = loopind,
    Zpath = Zpath,
    diverge = diverge,
    M_step_res = M_step_res,
    weights_new = weights_new
  )
  
  return(res)
}
