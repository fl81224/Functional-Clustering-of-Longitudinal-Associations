simu_onetime_realdata=function(weights_init=NULL,sigmasq_init=NULL,
                               inputlist,Rseed,N,p,samediff,inflation,save=TRUE,
                               max_outer_iter){

  max_inner_iter=250; inflation=100
  text=list();  bestres=list(); bestres_alpha1=list(); res=list()
  
  weights_init0=weights_init
  sigmasq_init0=sigmasq_init
  
  ind=1
  for(k in 2){
    {max_inner_iter=750}
    
    temp_alpha=list()
    for(alpha2 in c(0.01)){
      
      max_inner_iter=3000
      
      temp_RP=list(); inflat=100
      
      for(lRP_ind in 1:4){
        if(lRP_ind==1 ){
          weights_init = weights_init0
          sigmasq_init = sigmasq_init0
        }
        
        if(lRP_ind !=1){
          weights_init = temp_RP[[lRP_ind-1]]$res_alpha$weights_new
          sigmasq_init = temp_RP[[lRP_ind-1]]$res_alpha$sigmasq
        }
        temp_RP[[lRP_ind]]=EM_realdata(weights_init = weights_init,
                                       sigmasq_init = sigmasq_init,
                                       inputlist=inputlist[[lRP_ind]],Rseed=Rseed,K=k,task=Rseed,
                                       inflation=inflat,samediff = samediff,
                                       plot_M_step = FALSE,max_outer_iter =max_outer_iter,
                                       max_inner_iter=max_inner_iter,
                                       save=FALSE,alpha2 = alpha2)
      }
      
      temp_alpha[[ind]]=temp_RP[[which.min(sapply(temp_RP,function(x)x$updateAIC))]]
      ind=ind+1
      
      
    }
    bestres[[k]]=temp_alpha[[which.min(sapply(temp_alpha,function(x)x$updateAIC))]]
    bestres_alpha1[[k]]=temp_alpha[[1]]
    
  }
  
  res=list()
  res$res_alpha_select=bestres
  res$res_alpha_1=bestres_alpha1
  res$res_alpha=temp_alpha
  return(res)
}
####################################################################################################
simu_onetime_realdata_subsample=function(weights_init=NULL,sigmasq_init=NULL,
                                         inputlist,Rseed,N,p,samediff,inflation,save=TRUE,
                               max_outer_iter){

  max_inner_iter=250; inflation=100
  text=list();  bestres=list(); bestres_alpha1=list(); res=list()
  
  ind=1
  for(k in 2){
    {max_inner_iter=1500}
    
    temp_alpha=list()
    
    for(alpha2 in c(1)){
      
      temp_RP=list(); inflat=100
      
      for(lRP_ind in 1){
        
        
        if(lRP_ind !=1){
          weights_init = temp_RP[[lRP_ind-1]]$res_alpha$weights_new
          sigmasq_init = temp_RP[[lRP_ind-1]]$res_alpha$sigmasq
        }
        temp_RP[[lRP_ind]]=EM_realdata(weights_init = weights_init,
                                       sigmasq_init = sigmasq_init,
                                       inputlist=inputlist[[lRP_ind]],Rseed=Rseed,K=k,task=Rseed,
                                       inflation=inflat,samediff = samediff,
                                       plot_M_step = FALSE,max_outer_iter =max_outer_iter,
                                       max_inner_iter=max_inner_iter,
                                       save=FALSE,alpha2 = alpha2)
      }
      
      temp_alpha[[ind]]=temp_RP[[which.min(sapply(temp_RP,function(x)x$updateAIC))]]
      ind=ind+1
    }
    bestres[[k]]=temp_alpha[[which.min(sapply(temp_alpha,function(x)x$updateAIC))]]

  }
  return(bestres)
}
