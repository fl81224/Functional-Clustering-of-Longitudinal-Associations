simu_onetime_vsrp = function(task,
                             Rseed,
                             noise,
                             N,
                             numbers,
                             p,
                             samediff,
                             inflation,
                             save = TRUE,
                             depTX,
                             depFX = 0,
                             max_outer_iter) {
  inputlist = data_generation(
    Rseed = Rseed,
    noise = noise,
    N = N,
    numbers = numbers,
    p = p,
    time_reso = 10,
    nbeta_basis = 10,
    lambda_RP = c(0.5, 0.1, 0.001, 0.00001),
    depTX = depTX
  )
  
  inputlist4 = data_generation(
    Rseed = Rseed,
    noise = noise,
    N = N,
    numbers = numbers,
    p = p,
    time_reso = 10,
    nbeta_basis = 10,
    lambda_RP = c(0.5, 0.1, 0.001, 0.00001),
    depTX = depTX
  )
  
  max_inner_iter = 750
  if (p == 180) {
    max_inner_iter = 500
  }
  
  
  weights_init = sigmasq_init = NULL
  
  temp_alpha = list()
  bestres = list()
  
  k = 1
  ind = 1
  for (alpha2 in c(1, 0.1, 0.01, 0.001)) {
    temp_RP = list()
    
    inflat = 100
    
    for (lRP_ind in 1:3) {
      if (lRP_ind != 1) {
        weights_init = temp_RP[[lRP_ind - 1]]$res_alpha$weights_new
        sigmasq_init = temp_RP[[lRP_ind - 1]]$res_alpha$sigmasq
        max_outer_iter0 = 10
      }
      temp_RP[[lRP_ind]] = EM_simu(
        weights_init = weights_init,
        sigmasq_init = sigmasq_init,
        inputlist = inputlist[[lRP_ind]],
        Rseed = Rseed,
        K = k,
        task = Rseed,
        inflation = inflat,
        samediff = samediff,
        numbers = numbers,
        plot_M_step = FALSE,
        max_outer_iter = max_outer_iter,
        max_inner_iter = max_inner_iter,
        save = FALSE,
        alpha2 = alpha2
      )
    }
    temp_alpha[[ind]] = temp_RP[[which.min(sapply(temp_RP, function(x)
      x$updateAIC))]]
    ind = ind + 1
  }
  
  bestres[[k]] = temp_alpha[[which.min(sapply(temp_alpha, function(x)
    x$updateAIC))]]
  
  vs = which(sapply(bestres[[1]]$beta_est, function(x)
    sum(x ^ 2)) != 0)
  
  if (identical(vs, integer(0))) {
    var_important = matrix(0, length(bestres[[1]]$res_alpha$M_step_res$fit$lambda) , p)
    for (lam in 1:length(bestres[[1]]$res_alpha$M_step_res$fit$lambda)) {
      for (j in 1:p) {
        if (sum(bestres[[1]]$res_alpha$M_step_res$fit$beta[(3 + 10 * (j - 1)):(2 +
                                                                               10 * j) , lam]) ^ 2 != 0) {
          var_important[lam, j] = 1
        }
      }
    }
    ind = which(apply(var_important, 1, sum) != 0)[1]
    vs = which(var_important[ind, ] != 0)
  }
  
  bestres[[1]]$res_alpha$M_step_res$fit$beta
  
  inputlist = data_generation(
    Rseed = Rseed,
    noise = noise,
    N = N,
    numbers = numbers,
    p = p,
    time_reso = 10,
    nbeta_basis = 10,
    lambda_RP = c(0.999),
    depTX = depTX,
    vs = vs
  )
  
  
  ##################################################
  
  temp = list()
  for (k in 1:4) {
    temp[[k]] = EM_novs(
      inputlist = inputlist[[1]],
      Rseed = Rseed,
      K = k,
      task = Rseed,
      numbers = numbers,
      plot_M_step = FALSE,
      max_outer_iter = 20,
      inflation = 2,
      save = FALSE
    )
  }
  
  bestK = which.min(sapply(temp, function(x)
    x$updateBIC))
  
  res = list()
  res$bestres = temp[[bestK]]
  res$vs = vs
  res$res3 = temp[[3]]
  
  return(res)
}
