library(fda)
library(Matrix)
library(foreach)
library(doParallel)
library(aricode)
library(grpreg)
library(glmnet)
library(usmap)
library(ggplot2)


source("./grpreg2/function_grpreg_GXY.R")
source("./grpreg2/function_grpreg_ortho.R")
source("./grpreg2/grpreg3.R")

source("./simu/emsimuvs.R")
source("./simu/function_datageneration.R")

load("./realdata/newimputeXY.Rdata")
source("./realdata/emrealdata.R")
source("./realdata/realdata_reg_onetime.R")

logY_stroke=lapply(Y_stroke, function(x)log(x+0.1))

for(i in 1:3226)XData[[i]]=XData[[i]][-c(2:10,132)]

X_i=lapply(XData,function(x)do.call(cbind,x))
for(i in 1:3226){colnames(X_i[[i]])=names(XData[[1]])}
for(i in 1:3226){X_i[[i]]=cbind(X_i[[i]],rep(1,10));colnames(X_i[[i]])[125]="time"}


X_var=list()
for(i in 1:125){X_var[[i]]=sapply(X_i,function(x)x[,i])}
X_var_sd=lapply(X_var,function(x) (x-mean(x))/sd(x))
X_var_sd[[125]]=matrix(1,10,3226)
X_i_sd=X_i
for(i in 1:3226)for(j in 1:125)X_i_sd[[i]][,j]=X_var_sd[[j]][,i]

#####
###########################################################################################################
#############################################################################################################
time_reso=10; nbeta_basis=10; N=3226; Ts=10; p=125; norder=4

ind=1
inputlist=list()
for(lambda_RP in c(0.5,0.1,0.01,0.001)){
  temp=HQ_generation(time_reso = time_reso,nbeta_basis = nbeta_basis,X_i=X_i_sd,lambda_RP = lambda_RP,N=N,norder = norder)
  H=temp$H; Q=temp$Q; beta_basis=temp$beta_basis; L=temp$L; Linv=temp$Linv
  inputlist[[ind]]=list(N=N,Ts=Ts,Y=logY_stroke,X_i=X_i_sd,p=p,
                        H=H,Q=Q,beta_basis=beta_basis,L=L,Linv=Linv,
                        nbeta_basis=nbeta_basis,lRP=lambda_RP)
  ind=ind+1
}


cl=makeCluster(3)
registerDoParallel(cl)
res_fullsample_init=foreach(k=c(1:50)
                            ,.packages = c("fda","MFPCA","aricode","grpreg","glmnet") ) %dopar% 
  simu_onetime_realdata_subsample(inputlist=inputlist,
                                  Rseed=k,N=N,p=p,samediff="same",
                                  inflation=100,save=TRUE,max_outer_iter=30
  )
  
stopImplicitCluster()
stopCluster(cl)

best_init=which.min(sapply(res_fullsample_init,function(x)x[[2]]$updateBIC))
# 
res=simu_onetime_realdata(weights_init = res_fullsample_init[[best_init]][[2]]$res_alpha$weights_new,
                          sigmasq_init = res_fullsample_init[[best_init]][[2]]$res_alpha$sigmasq,
                          inputlist=inputlist,Rseed=1,N=N,p=p,samediff="same",
                          inflation=100,save=TRUE,max_outer_iter=20)





########################################################################  
