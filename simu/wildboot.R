#################################################################################
wild_bootstrap_onetime_simu = function(Rseedboot,
                                       inputlist,
                                       inputlist99,
                                       updateres,
                                       K = 3,
                                       task,
                                       save = T) {
  set.seed(Rseedboot)
  
  N = inputlist$N
  Ts = inputlist$Ts
  Y = inputlist$Y
  X_i = inputlist$X_i
  
  p = inputlist$p
  Q = inputlist$Q
  Linv = inputlist$Linv
  nbeta_basis = inputlist$nbeta_basis
  
  
  
  Z = updateres$res_alpha_select$bestres$res_alpha$Z_new
  alpha = updateres$res_alpha_select$bestres$res_alpha$M_step_res$alpha
  prop = updateres$res_alpha_select$bestres$res_alpha$M_step_res$prop
  
  
  alpha = updateres$res_alpha_select$bestres$res_alpha$M_step_res$alpha
  
  sigmasq = updateres$res_alpha_select$bestres$res_alpha$M_step_res$sigmasq
  
  prop = updateres$res_alpha_select$bestres$res_alpha$M_step_res$prop
  
  
  termloss_array = array(0, c(N, K, length(Ts)))
  fitted_array = array(0, c(N, K, length(Ts)))
  weights = matrix(0, N, K)
  
  for (i in 1:N) {
    termloss = matrix(0, K, length(Ts))
    for (k in 1:K) {
      termloss[k, ] = ((Y[[i]] - Q[[i]] %*%
                          alpha[[k]][-1] - alpha[[k]][1]) ^ 2)
      termloss_array[i, k, ] = as.vector((Y[[i]] - Q[[i]] %*%
                                            alpha[[k]][-1] - alpha[[k]][1]))
      fitted_array[i, k, ] = as.vector((Q[[i]] %*%
                                          alpha[[k]][-1] + alpha[[k]][1]))
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
  
  
  termloss_array_boot = termloss_array
  
  
  ind = NULL
  for (k in 1:3) {
    ind_of_k = which(Z == k)
    ind = c(ind, sample(ind_of_k, length(ind_of_k), replace = T))
    
  }
  
  
  ind = sort(ind)
  
  Y_boot = fitted_array
  
  for (k in 1:K) {
    Y_boot[, k, ] =  Y_boot[, k, ] + termloss_array_boot[ind, k, ]
  }
  
  
  
  k_sample = NULL
  for (i in 1:N) {
    k_sample = c(k_sample, sample(1:K, 1, prob = weights[i, ]))
  }
  
  Y_boot_list = list()
  for (i in 1:N) {
    Y_boot_list[[i]] = Y_boot[i, k_sample[i], ]
  }
  
  
  Y_boot = Y_boot_list
  
  
  inputlist_boot = inputlist
  inputlist_boot$Y = Y_boot
  
  weights_boot = updateres$res_alpha_select$bestres$res_alpha$weights_new
  
  for (i in 1:N) {
    ind = which.max(weights_boot[i, ])
    weights_boot[i, ] = 0
    weights_boot[i, ind] = 1
  }
  
  
  alpha2 = updateres$res_alpha_select$bestres$alpha2
  
  M_step_res = M_step_simu(
    weights = weights_boot,
    old_sigmasq = updateres$res_alpha_select$bestres$res_alpha$sigmasq,
    inputlist = inputlist_boot,
    K = K,
    lambda = NULL,
    lambdamax = 100,
    samediff = "same",
    plot = F,
    alpha2 = alpha2,
    maxiter = 3000
  )
  
  
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
  
  res = list(beta_est = beta_est)
  return(res)
}

##################################################################################################

#####################################################################################################################################

betatrue = list()
beta_reso = 100
Ts = seq(0, 1, , beta_reso)
for (i in 1:3)
  betatrue[[i]] = matrix(0, beta_reso, 3)

betatrue[[1]][, 2] = sin(0.5 * pi * Ts + 1.5 * pi) + (Ts)
betatrue[[1]][, 1] = sin(0.5 * pi * Ts + pi) - Ts
betatrue[[1]][, 3] = (Ts - 0.5) ^ 2

betatrue[[2]][, 2] = sin(2 * pi * Ts) - (Ts)
betatrue[[2]][, 1] = (cos(2 * pi * Ts) - 1) ^ 2
betatrue[[2]][, 3] = -1 * sin(2 * pi * Ts) + (-1)


betatrue[[3]][, 1] = -1 * sin(0.5 * pi * Ts + 1.5 * pi) - (Ts - 0.5) - 1
betatrue[[3]][, 2] = -1 * (cos(2 * pi * Ts) - 1) ^ 2
betatrue[[3]][, 3] = -1 * (-1 * sin(0.5 * pi * Ts + 1.5 * pi)) + 1

betanorm = 0
for (i in 1:3) {
  for (j in 1:3) {
    x = betatrue[[i]][, j]
    betatrue[[i]][, j] = (x) / sqrt(sum((x[1 + 11 * (0:9)]) ^ 2)) / 10
    betanorm = betanorm + sum(betatrue[[i]][, j] ^ 2)
  }
}
