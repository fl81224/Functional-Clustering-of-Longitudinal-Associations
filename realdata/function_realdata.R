alphaupdate_rd_vs = function(Z, inputlist, K, lambda = NULL) {
  N = inputlist$N
  Ts = inputlist$Ts
  Y = inputlist$Y
  X_i = inputlist$X_i
  
  p = inputlist$p
  Q = inputlist$Q
  Linv = inputlist$Linv
  nbeta_basis = inputlist$nbeta_basis
  
  Y_vec = unlist(Y)
  Q_mat = matrix(0, N * 10, p * nbeta_basis * K + K)
  for (i in 1:N) {
    kind = Z[i]
    Q_mat[(1 + 10 * (i - 1)):(10 * i),
          (1 + (p * nbeta_basis + 1) * (kind - 1)):((p * nbeta_basis + 1) *
                                                      (kind))] = c(matrix(1, 10, 1), Q[[i]])
  }
  
  
  
  group = NULL
  for (k in 1:K) {
    group = c(group, ((p + 1) * (k - 1) + 1))
    for (i in 1:p) {
      group = c(group, rep((p + 1) * (0) + i + 1 , nbeta_basis))
    }
  }
  
  if (is.null(lambda)) {
    fit = tryCatch({
      grpreg(
        X = Q_mat,
        y = Y_vec,
        group = group,
        penalty = "grSCAD",
        alpha = 0.5,
        max.iter = 2000
      )
    },
    error = function(e)
      1)
  } else{
    fit = tryCatch({
      grpreg(
        X = Q_mat,
        y = Y_vec,
        group = group,
        penalty = "grSCAD",
        alpha = 0.5,
        max.iter = 5000,
        lambda = lambda
      )
    },
    error = function(e)
      1)
  }
  
  updateBIC = matrix(0, length(fit$lambda), 1)
  Q_bigmat = kronecker(diag(1, K), cbind(matrix(1, N * 10, 1), do.call(rbind, Q)))
  
  
  Ypre = cbind(rep(1, N * Ts * K), Q_bigmat) %*% fit$beta
  
  
  
  for (lam in 1:length(fit$lambda)) {
    alpha = list()
    for (k in 1:K) {
      alpha[[k]] = fit$beta[(2 + (nbeta_basis * p + 1) * (k - 1)):((nbeta_basis *
                                                                      p + 1) * (k) + 1), lam]
      alpha[[k]][1] = alpha[[k]][1] + fit$beta[1, lam]
    }
    
    
    termloss = array(((rep(unlist(
      Y
    ), K) - Ypre[, lam]) ^ 2), c(Ts, N, K))
    termloss_indi = t(apply(termloss, c(2, 3), mean))
    
    
    
    prop = matrix(apply(termloss_indi, 2, function(x)
      x / sum(x)), K, N)
    prop = apply(prop, 1, mean)
    
    new_sigmasq = matrix(0, K, 1)
    for (k in 1:K) {
      new_sigmasq[k] = sum(termloss_indi[k, ] * as.numeric(Z == k)) / sum(as.numeric(Z ==
                                                                                       k))
    }
    
    RSS = (new_sigmasq)
    
    termloss_divide = apply(termloss, c(1, 2), function(x)
      log(sum(prop *
                (
                  exp(-1 * as.vector(x) / (2 * RSS)) /
                    sqrt(RSS)
                ))))
    ll = sum(termloss_divide)
    
    updateBIC[lam] = -2 * ll + sum((fit$df[lam]) * log(10 * N))
  }
  
  
  lastone = which.min(updateBIC)
  alpha = list()
  for (k in 1:K) {
    alpha[[k]] = fit$beta[(2 + (nbeta_basis * p + 1) * (k - 1)):((nbeta_basis *
                                                                    p + 1) * (k) + 1), lastone]
    alpha[[k]][1] = alpha[[k]][1] + fit$beta[1, lastone]
  }
  
  b = alpha
  for (k in 1:K) {
    for (j in 1:p) {
      var_index = (2 + (j - 1) * nbeta_basis):(1 + (j) * nbeta_basis)
      b[[k]][var_index] = Linv %*% alpha[[k]][var_index]
    }
  }
  
  termloss = array(((rep(unlist(
    Y
  ), K) - Ypre[, lastone]) ^ 2), c(Ts, N, K))
  termloss_indi = t(apply(termloss, c(2, 3), mean))
  prop = matrix(apply(termloss_indi, 2, function(x)
    x / sum(x)), K, N)
  prop = apply(prop, 1, mean)
  
  Ypre_all = Ypre[, lastone]
  Rsq = matrix(0, K, K)
  for (k1 in 1:K) {
    kind = Z == k1
    Y_true = matrix(unlist(Y), 10, N)[, kind]
    for (k2 in 1:K) {
      Y_pre = matrix(Ypre_all[(1 + (k2 - 1) * N * 10):(k2 * N * 10)], 10, N)[, kind]
      RSS = mean((Y_true - Y_pre) ^ 2)
      ERR = mean((Y_true - mean(Y_true)) ^ 2)
      Rsq[k1, k2] = 1 - RSS / ERR
    }
  }
  new_sigmasq = matrix(0, K, 1)
  for (k in 1:K) {
    new_sigmasq[k] = sum(termloss_indi[k, ] * as.numeric(Z == k)) / sum(as.numeric(Z ==
                                                                                     k))
  }
  RSS = (new_sigmasq)
  termloss_divide = apply(termloss, c(1, 2), function(x)
    log(sum(prop *
              (
                exp(-1 * as.vector(x) / (2 * RSS)) /
                  sqrt(RSS)
              ))))
  ll = sum(termloss_divide)
  
  print(ll)
  res = list(
    alpha = alpha,
    b = b,
    df = fit$df[lastone],
    RSS = RSS,
    Rsq = Rsq,
    BIC = updateBIC[lastone],
    ll = ll,
    fit = fit,
    bestone = lastone
  )
  return(res)
}
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
clustervs_rd_onetime = function(Rseed, K, inputlist, task, save = FALSE) {
  set.seed(seed = Rseed)
  
  N = inputlist$N
  Ts = inputlist$Ts
  Y = inputlist$Y
  X_i = inputlist$X_i
  
  p = inputlist$p
  Q = inputlist$Q
  Linv = inputlist$Linv
  nbeta_basis = inputlist$nbeta_basis
  lambda = seq(0.005, 0.0001, , 10)
  
  Z_old = rep(2, N)
  Z_new = init_label(Rseed = Rseed, N = N, K = K)
  
  loopind = 1
  max_iter = 100
  Zpath = matrix(0, max_iter, N)
  
  coefnorm_old = 0
  coefnorm_new = 100
  
  aa = Sys.time()
  while (abs(coefnorm_old - coefnorm_new) / (coefnorm_new + coefnorm_old) >
         1e-3 & loopind < max_iter) {
    if (length(table(Z_new)) != K) {
      break
    }
    Z_old = Z_new
    coefnorm_old = coefnorm_new
    
    updateres = alphaupdate_rd_vs(Z = Z_new,
                                  inputlist = inputlist,
                                  K = K)
    updatealpha = updateres$alpha
    coefnorm_new = norm(unlist(updatealpha))
    
    if (is.null(updateres$df)) {
      BIC = 1e10
      break
    }
    
    Z_new = updateZ(
      N = N,
      Y = Y,
      Q = Q,
      K = K,
      updatealpha = updatealpha
    )
    Z_new = Z_new$temp_vec
    Zpath[loopind, ] = Z_new
    loopind = loopind + 1
  }
  
  updateres = alphaupdate_rd_vs(Z = Z_new, inputlist = inputlist, K = K)
  
  updateBIC = 1e10
  
  if (length(table(Z_new)) == K & (is.null(updateres$df) != 1)) {
    BIC_LL = ComputeBIC(
      Z_new = Z_new,
      updateres = updateres,
      N = N,
      K = K,
      Y = Y,
      Q = Q
    )
    updateBIC = BIC_LL$BIC
    
    beta_est = beta_eval_rd(
      nbeta_basis = nbeta_basis,
      K = K,
      updateres = updateres,
      N = N,
      Z_new = Z_new,
      beta_basis = inputlist$beta_basis,
      p = p
    )
    
    res = list(
      beta_est = beta_est,
      updateBIC = updateBIC,
      Z = Z_new,
      K = K,
      intercept = sapply(updateres$alpha, function(x)
        x[1]),
      Rsq = updateres$Rsq,
      loopind = loopind,
      Zpath = Zpath,
      lambda = lambda,
      diverge = FALSE
    )
    
    return(res)
  } else{
    updateBIC = 1e10
    
    res = list(
      updateBIC = updateBIC,
      Z = Z_new,
      K = K,
      loopind = loopind,
      Zpath = Zpath,
      lambda = lambda,
      diverge = is.null(updateres$df)
    )
    return(res)
  }
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
