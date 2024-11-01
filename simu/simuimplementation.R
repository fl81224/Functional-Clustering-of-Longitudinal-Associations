rm(list = ls())
library(fda)
library(Matrix)
library(foreach)
library(doParallel)
library(aricode)
library(grpreg)
library(glmnet)
library(corpcor)
library(flexmix)

source("./grpreg2/function_grpreg_GXY.R")
source("./grpreg2/function_grpreg_ortho.R")
source("./grpreg2/grpreg2.R")

source("./simu/emsimunovs.R")
source("./simu/emsimuvs.R")
source("./simu/function_datageneration.R")
source("./simu/simu_onetime.R")
source("./simu/simu_onetime_novs.R")
source("./simu/simu_onetime_vsrp.R")
source("./simu/othermethod.R")
aa = Sys.time()

set.seed(418)
simuind = sample(1:1e5, 100)
simuind = 1:100

clustertimes = 100



cl = makeCluster(1)
registerDoParallel(cl)
for (N in c(180)) {
  for (p in c(240)) {
    for (noise in c(1 / 12)) {
      for (depTX in c(0.8)) {
        inflation = 2
        numbers = rep(N / 3, 3)
        
        
        samediff = "same"
        res = foreach(
          k = c(1:100),
          .packages = c("fda", "aricode", "grpreg", "glmnet", "flexmix")
        ) %dopar%
          simu_onetime(
            task = k,
            Rseed = simuind[k],
            noise = rep(noise, 3),
            N = N,
            numbers = numbers,
            p = p,
            samediff = samediff,
            inflation = 2,
            depTX = depTX,
            max_outer_iter = 30,
            k_seq = c(1:4),
            save = TRUE
          )
        
        
        samediff = "diff"
        res = foreach(
          k = c(1:100),
          .packages = c("fda", "aricode", "grpreg", "glmnet", "flexmix")
        ) %dopar%
          simu_onetime(
            task = k,
            Rseed = simuind[k],
            noise = rep(noise, 3),
            N = N,
            numbers = numbers,
            p = p,
            samediff = samediff,
            inflation = 2,
            depTX = depTX,
            max_outer_iter = 30,
            k_seq = c(1:4),
            save = TRUE
          )
        
        
        res = foreach(
          k = c(1:100),
          .packages = c("fda", "aricode", "grpreg", "glmnet", "flexmix")
        ) %dopar%
          simu_onetime_novs(
            task = k,
            Rseed = simuind[k],
            noise = rep(noise, 3),
            N = N,
            numbers = numbers,
            p = p,
            samediff = samediff,
            inflation = 2,
            depTX = depTX,
            max_outer_iter = 30,
            save = TRUE
          )
        
        
        res = foreach(
          k = c(1:100),
          .packages = c("fda", "aricode", "grpreg", "glmnet", "flexmix")
        ) %dopar%
          simu_onetime_vsrp(
            task = k,
            Rseed = simuind[k],
            noise = rep(noise, 3),
            N = N,
            numbers = numbers,
            p = p,
            samediff = samediff,
            inflation = 2,
            depTX = depTX,
            max_outer_iter = 30
          )
        
        
        res = foreach(
          k = c(1:100),
          .packages = c("fda", "aricode", "grpreg", "glmnet", "flexmix")
        ) %dopar%
          flexmix_simu(
            Rseed = simuind[k],
            noise = rep(noise, 3),
            N = N,
            numbers = numbers,
            p = p,
            depTX = depTX,
            depFX = 0
          )
        
        
      }
    }
  }
  
}
stopImplicitCluster()
stopCluster(cl)
