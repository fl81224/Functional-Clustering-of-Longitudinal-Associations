init_label = function(Rseed, N, K) {
  set.seed(seed = Rseed)
  init_label = sample(1:N, N, replace = FALSE)
  len = ceiling(N / K)
  temp_vec = rep(1, N)
  for (k in 1:K)
    temp_vec[init_label[(1 + len * (k - 1)):min((len * k), N)]] = k
  return(temp_vec)
}
############################################################################################
E_step_simu = function(inputlist, M_step_res, K) {
  N = inputlist$N
  Ts = inputlist$Ts
  Y = inputlist$Y
  X_i = inputlist$X_i
  
  p = inputlist$p
  Q = inputlist$Q
  Linv = inputlist$Linv
  nbeta_basis = inputlist$nbeta_basis
  
  alpha = M_step_res$alpha
  sigmasq = M_step_res$sigmasq
  prop = M_step_res$prop
  
  weights = matrix(0, N, K)
  for (i in 1:N) {
    termloss = matrix(0, K, length(Ts))
    for (k in 1:K) {
      termloss[k, ] = ((Y[[i]] - Q[[i]] %*%
                          alpha[[k]][-1] - alpha[[k]][1]) ^ 2)
      weights[i, k] = prop[k] * prod(exp(-1 * as.vector(termloss[k, ]) /
                                           (2 * sigmasq[k])) /
                                       sqrt(sigmasq[k]))
    }
    weights[i, ] = weights[i, ] / sum(weights[i, ])
    if (is.nan(weights[i, 1])) {
      log_weights_temp = matrix(0, 1, K)
      for (k in 1:K) {
        log_weights_temp[k] = log(prop[k]) +
          sum((-1 * as.vector(termloss[k, ]) / (2 * sigmasq[k])) -
                log(sqrt(sigmasq[k])))
      }
      log_weights_temp = log_weights_temp - max(log_weights_temp)
      weights[i, ] = exp(log_weights_temp)
      weights[i, ] = weights[i, ] / sum(weights[i, ])
    }
  }
  
  return(weights)
}
########################################################################################
M_step_simu = function(weights,
                       old_sigmasq,
                       inputlist,
                       K,
                       lambda = NULL,
                       init = FALSE,
                       maxiter = NULL,
                       lambdamax = 1e10,
                       samediff,
                       alpha2 = 0.5,
                       plot = TRUE,
                       criteria = "BIC",
                       penalty = "grSCAD",
                       lambdamin = NULL) {
  aa = Sys.time()
  
  N = inputlist$N
  Ts = inputlist$Ts
  Y = inputlist$Y
  X_i = inputlist$X_i
  
  p = inputlist$p
  Q = inputlist$Q
  Linv = inputlist$Linv
  nbeta_basis = inputlist$nbeta_basis
  
  Ts = length(Ts)##CARE in simu ts=1:10 in rd ts=10
  
  Y_mean = mean(unlist(Y))
  Y_vec = rep(unlist(Y) - Y_mean, K)
  Q_mat = kronecker((diag(1, K)), (cbind(matrix(1, N * (
    Ts
  ), 1), do.call(rbind, Q))))
  
  prop = apply(weights, 2, mean)
  
  WLS_weights = NULL
  for (k in 1:K) {
    WLS_weights = c(WLS_weights, rep(weights[, k] / old_sigmasq[k], each = 10))
  }
  WLS_weights = WLS_weights / max(WLS_weights)
  
  Y_vec_reweight = Y_vec * sqrt(WLS_weights)
  Q_mat_reweight = Q_mat
  for (i in 1:length(WLS_weights)) {
    Q_mat_reweight[i, ] = Q_mat_reweight[i, ] * sqrt(WLS_weights[i])
  }
  
  if (samediff == "diff") {
    lambdaweight = NULL
    group = NULL
    for (k in 1:K) {
      group = c(group, ((p + 1) * (k - 1) + 1))
      lambdaweight = c(lambdaweight, sqrt(prop[k]))
      
      for (i in 1:p) {
        group = c(group, rep((p + 1) * (k - 1) + i + 1 , nbeta_basis))
        lambdaweight = c(lambdaweight, sqrt(prop[k]) * sqrt(nbeta_basis))
      }
    }
    lambdaweight = NULL
  }
  if (samediff == "same") {
    lambdaweight = NULL
    group = NULL
    for (k in 1:K) {
      group = c(group, ((p + 1) * (k - 1) + 1))
      
      for (i in 1:p) {
        group = c(group, rep((p + 1) * (0) + i + 1 , nbeta_basis))
      }
    }
  }
  
  if (is.null(lambda)) {
    fit = tryCatch({
      grpreg2(
        X = as.matrix(Q_mat_reweight),
        y = Y_vec_reweight,
        group = group,
        penalty = penalty,
        alpha = alpha2,
        max.iter = maxiter,
        family = "gaussian",
        group.multiplier = lambdaweight
      )
    },
    error = function(e)
      1)
  }
  
  
  if (!is.null(lambda)) {
    fit = tryCatch({
      grpreg2(
        X = as.matrix(Q_mat_reweight),
        y = Y_vec_reweight,
        group = group,
        penalty = penalty,
        alpha = alpha2,
        max.iter = maxiter,
        family = "gaussian",
        group.multiplier = lambdaweight,
        lambda = lambda
      )
    },
    error = function(e)
      1)
  }
  
  if (!is.list(fit)) {
    res = list(BIC = 1e10, AIC = 1e10)
    return(res)
  }
  
  updateBIC = matrix(0, length(fit$lambda), 1)
  updateAIC = matrix(0, length(fit$lambda), 1)
  llres = matrix(0, length(fit$lambda), 1)
  
  
  Ypre = cbind(rep(1, N * Ts * K), Q_mat) %*% fit$beta + Y_mean
  for (lam in 1:length(fit$lambda)) {
    # transform alpha
    alpha = list()
    lambda2 = fit$lambda[lam] * (1 - alpha2)
    for (k in 1:K) {
      alpha[[k]] = fit$beta[(2 + (nbeta_basis * p + 1) * (k - 1)):(1 + (nbeta_basis *
                                                                          p + 1) * (k)), lam] * (1 + 2 * lambda2 * fit$group.multiplier[2])
      alpha[[k]][1] = alpha[[k]][1] * (1 + 2 * lambda2) + Y_mean
    }
    
    termloss = array(((unlist(Y) - Ypre[, lam]) ^ 2), c(Ts, N, K))
    termloss_indi = t(apply(termloss, c(2, 3), mean))
    
    new_sigmasq = matrix(0, K, 1)
    for (k in 1:K) {
      new_sigmasq[k] = sum(termloss_indi[k, ] * weights[, k]) / sum(weights[, k])
      if (is.nan(new_sigmasq[k]))
        new_sigmasq[k] = 1e-8
    }
    
    
    termloss_divide = apply(termloss, c(1, 2), function(x)
      log(sum(prop *
                (
                  exp(-1 * as.vector(x) / (2 * new_sigmasq)) /
                    sqrt(new_sigmasq)
                ))))
    ll = sum(termloss_divide)
    
    llres[lam] = ll
    updateBIC[lam] = -2 * ll + ((fit$df[lam]) * log(Ts * N))
    updateAIC[lam] = -2 * ll + 2 * (fit$df[lam])
    
  }
  
  if (criteria == "BIC") {
    lastone = which.min(updateBIC)
  }
  if (criteria == "AIC") {
    lastone = which.min(updateAIC)
  }
  if (fit$lambda[lastone] > lambdamax) {
    lastone = which(fit$lambda < lambdamax)[1]
    if (is.na(lastone))
      lastone = length(fit$lambda)
  }
  if (init == TRUE)
    lastone = length(fit$lambda)
  
  alpha = list()
  for (k in 1:K) {
    alpha[[k]] = fit$beta[(2 + (nbeta_basis * p + 1) * (k - 1)):(1 + (nbeta_basis *
                                                                        p + 1) * (k)), lastone] * (1 + 2 * lambda2 * fit$group.multiplier[2])
    alpha[[k]][1] = alpha[[k]][1] * (1 + 2 * lambda2) + Y_mean
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
    if (is.nan(new_sigmasq[k]))
      new_sigmasq[k] = 1e-8
  }
  Z = apply(weights, 1, which.max)
  
  Ypre_all = Ypre[, lastone] - Y_mean
  
  Rsq = matrix(0, K, K)
  for (k1 in 1:K) {
    kind = Z == k1
    Y_true = matrix(Y_vec[1:(N * (Ts))], (Ts), N)[, kind]
    for (k2 in 1:K) {
      Y_pre = matrix(Ypre_all[(1 + (k2 - 1) * N * (Ts)):(k2 * N * (Ts))], (Ts), N)[, kind]
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
    BIClist = updateBIC,
    AIClist = updateAIC,
    bestone = lastone,
    alpha2 = alpha2
  )
  return(res)
}
###########################################################################
EM_simu = function(inputlist,
                   weights_init = NULL,
                   sigmasq_init = NULL,
                   Rseed,
                   K,
                   task,
                   save = FALSE,
                   inflation = 2,
                   samediff,
                   numbers,
                   plot_M_step = TRUE,
                   plot_converge_res = TRUE,
                   max_inner_iter = 1000,
                   max_outer_iter = 100,
                   alpha2 = alpha2) {
  N = inputlist$N
  Ts = inputlist$Ts
  Y = inputlist$Y
  X_i = inputlist$X_i
  
  p = inputlist$p
  Q = inputlist$Q
  Linv = inputlist$Linv
  nbeta_basis = inputlist$nbeta_basis
  lRP = inputlist$lRP
  
  random_init = F
  set.seed(seed = Rseed)
  if (is.null(weights_init)) {
    random_init = T
    weights_init = matrix(0, N, K)
    init_label = sample(1:N, N, replace = FALSE)
    len = ceiling(N / K)
    for (k in 1:K)
      weights_init[init_label[(1 + len * (k - 1)):min((len * k), N)] , k] = 1
  }
  if (is.null(sigmasq_init)) {
    sigmasq_init = rep(1, K)
  }
  if (alpha2 < 10) {
    M_step_res0 = M_step_simu(
      weights = weights_init,
      old_sigmasq = sigmasq_init,
      inputlist = inputlist,
      K = K,
      lambda = NULL,
      maxiter = max_inner_iter,
      init = random_init,
      samediff = samediff,
      plot = plot_M_step,
      alpha2 = 1
    )
    lambda_seq = M_step_res0$fit$lambda
    M_step_res = M_step_simu(
      weights = weights_init,
      old_sigmasq = sigmasq_init,
      inputlist = inputlist,
      K = K,
      lambda = lambda_seq,
      maxiter = max_inner_iter,
      init = random_init,
      samediff = samediff,
      plot = plot_M_step,
      alpha2 = alpha2
    )#0.007148256
  }
  #######
  steps = 0
  
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
  M_step_res_in = list()
  
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
    Z_old = apply(weights_old, 1, which.max)
    
    
    coefnorm_old = coefnorm_new
    weights_old = weights_new
    lambdaold = lambdanew
    
    weights_new = E_step_simu(inputlist = inputlist,
                              M_step_res = M_step_res,
                              K = K)
    Z_new = apply(weights_new, 1, which.max)
    M_step_res = M_step_simu(
      weights = weights_new,
      old_sigmasq = M_step_res$sigmasq,
      inputlist = inputlist,
      K = K,
      maxiter = max_inner_iter,
      lambdamax = inflation * M_step_res$fit$lambda[M_step_res$bestone],
      samediff = samediff,
      plot = plot_M_step,
      alpha2 = alpha2,
      lambda = lambda_seq
    )
    coefnorm_new = unlist(M_step_res$alpha)
    lambdanew = M_step_res$fit$lambda[M_step_res$bestone]
    
    Zpath[loopind, ] = Z_new
    
    M_step_res_in[[loopind]] = list(
      M_step_res = M_step_res,
      weights_new = weights_new,
      sigmasq = M_step_res$sigmasq,
      Z_new = Z_new,
      loopind = loopind
    )
    loopind = loopind + 1
  }
  
  steps = loopind
  diverge = FALSE
  if (loopind >= max_iter)
    diverge = TRUE
  
  res_alpha = M_step_res_in[[loopind - 1]]
  res_alpha$diverge = diverge
  
  
  Zind = matrix(0, 3, 1)
  Zind[1] = as.numeric(names(which.max(table(Z_new[1:numbers[1]]))))
  Zind[2] = as.numeric(names(which.max(table(Z_new[(1 + numbers[1]):(numbers[1] +
                                                                       numbers[2])]))))
  Zind[3] = as.numeric(names(which.max(table(Z_new[(1 + numbers[1] + numbers[2]):(numbers[1] +
                                                                                    numbers[2] + numbers[3])]))))
  if (K == 3 & length(table(Z_new)) == 3 & length(table(Zind)) == 3) {
    weightstemp = res_alpha$weights_new
    
    sigmasqtemp = res_alpha$sigmasq
    
    for (j in 1:3) {
      res_alpha$weights_new[, j] = weightstemp[, Zind[j]]
      res_alpha$sigmasq[j] = sigmasqtemp[Zind[j]]
    }
    
    M_step_res = M_step_simu(
      weights = res_alpha$weights_new,
      old_sigmasq = res_alpha$sigmasq,
      samediff = samediff,
      inputlist = inputlist,
      K = K,
      maxiter = max_inner_iter,
      lambdamax = 10,
      plot = plot_M_step,
      alpha2 = alpha2,
      lambda = lambda_seq
    )
    Z_new = apply(res_alpha$weights_new, 1, which.max)
    res_alpha$M_step_res = M_step_res
    res_alpha$Z_new = Z_new
  }
  
  
  BIC = res_alpha$M_step_res$BIC
  AIC = res_alpha$M_step_res$AIC
  beta_est = list()
  eval.beta = eval.basis(seq(0, 1, , 100), inputlist$beta_basis)
  
  for (j in 1:p) {
    beta_est[[j]] = matrix(0, 100, K)
    for (k in 1:K) {
      beta_est[[j]][, k] =
        res_alpha$M_step_res$b[[k]][(2 + (j - 1) * nbeta_basis):(1 + j *
                                                                   nbeta_basis)] %*% t(eval.beta)
    }
  }
  
  res = list(
    beta_est = beta_est,
    updateBIC = BIC,
    updateAIC = AIC,
    K = K,
    intercept = sapply(res_alpha$M_step_res$alpha, function(x)
      x[1]),
    Rsq = res_alpha$M_step_res$Rsq,
    loopind = loopind,
    Zpath = Zpath,
    diverge = diverge,
    res_alpha = res_alpha,
    steps = steps,
    alpha2 = alpha2,
    lRP = lRP
  )
  
  return(res)
}
#########################################################################################
beta_eval_rd = function(nbeta_basis,
                        K,
                        updateres,
                        N,
                        Z_new,
                        beta_basis,
                        p) {
  eval.beta = eval.basis(seq(1, 10, , 100), beta_basis)
  beta_est = list()
  
  Zind = 1:K
  
  par(mfrow = c(3, K))
  for (j in 1:p) {
    beta_est[[j]] = matrix(0, 100, K)
    for (k in 1:K) {
      beta_est[[j]][, Zind[k]] =
        updateres$b[[Zind[k]]][(2 + (j - 1) * nbeta_basis):(1 + j * nbeta_basis)] %*%
        t(eval.beta)
    }
  }
  return(beta_est)
}
