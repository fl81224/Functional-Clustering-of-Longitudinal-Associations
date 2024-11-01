cl = makeCluster(100)
registerDoParallel(cl)

for (N in c(180)) {
  for (p in c(90)) {
    for (noise in c(1 / 12)) {
      for (depTX in c(0.4)) {
        inflation = 2
        numbers = rep(N / 3, 3)
        
        Rseed_ind = which(sapply(res, function(x)
          x$res_alpha_select$bestres$K) == 3)
        
        
        bootres_beta_list = list()
        for (Rseed in Rseed_ind) {
          cat(Rseed)
          
          inputlist = data_generation(
            Rseed = simuind[Rseed_ind],
            noise = rep(noise, 3),
            N = N,
            numbers = numbers,
            p = p,
            time_reso = 10,
            nbeta_basis = 10,
            lambda_RP = c(0.5, 0.1, 0.001, 0.00001,0.9),
            depTX = depTX
          )
          
          inputlist99 = data_generation(
            Rseed = Rseed,
            noise = rep(noise, 3),
            N = N,
            numbers = numbers,
            p = p,
            time_reso = 10,
            nbeta_basis = 10,
            lambda_RP = c(0.999),
            depTX = depTX
          )
          
          inputlist99 = inputlist99[[1]]
          
          updateres = res[[Rseed]]
          
          lRP_ind = which(c(0.5, 0.1, 0.001, 0.00001) == updateres$res_alpha_select$bestres$lRP)
 
          
          
          bootres = foreach(
            k = c(1:300),
            .packages = c(
              "fda",
              "aricode",
              "grpreg",
              "glmnet",
              "flexmix",
              "Hmisc",
              "dplyr",
              "fdapace"
            )
          ) %dopar%
            
            tryCatch(
              {
                
                wild_bootstrap_onetime_simu(
                  Rseedboot = k,
                  inputlist0=inputlist[[lRP_ind]],
                  inputlist = inputlist[[lRP_ind]],
                  inputlist99 = inputlist99,
                  updateres = updateres,
                  task = k
                )
              },
              error=function(e )1
            )
        
            
          
          bootres_beta = lapply(bootres, function(x)
            x[[1]])
          
          bootres_beta_list[[Rseed]] = bootres_beta
          
        }
      }
    }
  }
  
}
stopImplicitCluster()
stopCluster(cl)












        