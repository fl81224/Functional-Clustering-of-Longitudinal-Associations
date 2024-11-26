simu_onetime = function(task,
                        Rseed,
                        noise,
                        N,
                        numbers,
                        p,
                        samediff,
                        inflation,
                        save = TRUE,
                        depTX,
                        depFX,
                        max_outer_iter = 20,
                        k_seq = 1:4) {
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
    nbeta_basis = 4,
    lambda_RP = c(0.1),
    depTX = depTX
  )
  
  
  max_inner_iter = 1000
  text = list()
  bestres = list()
  bestres_alpha1 = list()
  res3 = list()
  
  max_outer_iter0 = max_outer_iter
  
  for (k in k_seq) {
    ind = 1
    temp_alpha = list()
    for (alpha2 in c(1, 0.1, 0.01)) {
      temp_RP = list()
      
      inflat = 100
      
      for (lRP_ind in 1:2) {
        cat("K=",k," alpha2=", alpha2, " LRP=", lRP_ind, '\n')
        
        weights_init = sigmasq_init = NULL
        
        if (k != 1 & alpha2 == 1 & lRP_ind == 1) {
          inflat = inflation
          max_outer_iter0 = max_outer_iter
          temp_start = EM_simu(
            inputlist = inputlist4[[1]],
            Rseed = Rseed,
            K = k,
            task = Rseed,
            inflation = inflat,
            samediff = samediff,
            numbers = numbers,
            plot_M_step = FALSE,
            max_outer_iter = max_outer_iter0,
            max_inner_iter = max_inner_iter,
            save = FALSE,
            alpha2 = alpha2
          )
          weights_init = temp_start$res_alpha$weights_new
          sigmasq_init = temp_start$res_alpha$sigma_sq
        }
        
        #Path Fitting
        if (lRP_ind == 1 & alpha2 != 1) {
          weights_init = temp_alpha[[ind - 1]]$res_alpha$weights_new
          sigmasq_init = temp_alpha[[ind - 1]]$res_alpha$sigmasq
          max_outer_iter0 = 10
        }
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
          max_outer_iter = max_outer_iter0,
          max_inner_iter = max_inner_iter,
          save = FALSE,
          alpha2 = alpha2
        )
        
        
        if (sum(apply(temp_RP[[lRP_ind]]$res_alpha$weights_new, 2, mean) <
                0.05) > 0) {
          temp_RP = list()
          temp_RP[[1]] = list()
          temp_RP[[1]]$updateAIC = 1e10
          temp_RP[[1]]$updateBIC = 1e10
          break
        }
        
        
      }
      temp_alpha[[ind]] = temp_RP[[which.min(sapply(temp_RP, function(x)
        x$updateAIC))]]
      
      
      if (temp_alpha[[ind]]$updateAIC == 1e10) {
        break
      }
      ind = ind + 1
    }
    bestres[[k]] = temp_alpha[[which.min(sapply(temp_alpha, function(x)
      x$updateAIC))]]
    bestres_alpha1[[k]] = temp_alpha[[1]]
  }
  
  res = list()
  bestK = which.min(sapply(bestres, function(x)
    x$updateBIC))
  res$res_alpha_select$bestres = bestres[[bestK]]
  res$res_alpha_select$res3 = bestres[[3]]
  res$res_alpha_1$bestres = bestres_alpha1[[bestK]]
  res$res_alpha_1$res3 = bestres_alpha1[[3]]
  
  return(res)
}
