#LE importance sampling
LEisamp=function(modelresults, nreps=1000) {

  #get E, tau, theta values from best model
  E=modelresults$modelstats$E
  tau=modelresults$modelstats$tau
  theta=modelresults$modelstats$theta
  hyperpars=data.frame(bestE=E, bestTau=tau, bestTheta=theta, 
                       Emax = modelresults$modelstats$Emax, taumax = modelresults$modelstats$taumax)
  #get best model form
  modelform=modelresults$form
  modnum=case_when(modelform =="ut-ut" ~ 1,
                   modelform =="log-log" ~ 2,
                   modelform =="fd-ut" ~ 3,
                   modelform =="gr-ut" ~ 4,
                   modelform =="gr-log" ~ 5)

  #rerun to get covariance matrices
  modelresults2=smap_model_options(data=modelresults$resultsdf, 
                                   hypars=hyperpars, model=modnum, y="Obs_abund")
  #covariance matrices
  covm=modelresults2$smap_results$smap_coefficient_covariances[[1]]
  resultsdf=modelresults2$resultsdf
  ss_opt=sum((resultsdf$obs-resultsdf$pred)^2,na.rm = T)
  n=length(which(!is.na(resultsdf$pred)))
  predictors=cbind.data.frame(modelresults2$ser_lags[-1],int=1)

  #resampling
  coefcols=c(paste0("c_",1:E),"c_0")
  coefs_orig=resultsdf[,coefcols]
  resultsdfcopy=resultsdf
  modelresultscopy=modelresults
  ss=numeric(length = nreps)
  JLE=numeric(length = nreps)
  
  resamples=array(NA, dim = c(dim(coefs_orig),nreps))
  for(i in 1:nrow(resultsdfcopy)) {
    if(!is.na(coefs_orig[i,1])) {
      resamples[i,,]=t(MASS::mvrnorm(n=nreps,mu=as.numeric(coefs_orig[i,]),Sigma=covm[[i]]))
    }
  }
  for(j in 1:nreps) {
    resultsdfcopy[,coefcols]=resamples[,,j]
    coefs_new=resultsdfcopy[,coefcols]
    resultsdfcopy$pred=rowSums(predictors*coefs_new)
    ss[j]=sum((resultsdfcopy$obs-resultsdfcopy$pred)^2,na.rm = T)
    modelresultscopy$resultsdf=resultsdfcopy

    jacobians=getJacobians(modelresultscopy)
    JLE[j]=LEglobal(modelresultscopy, jacobians)
  }
  # plot(JLE)
  # for(j in 1:nreps) {
  #   for(i in 1:nrow(resultsdfcopy)) {
  #     if(!is.na(coefs_orig[i,1])) {
  #       resultsdfcopy[i,c(paste0("c_",1:E),"c_0")]=MASS::mvrnorm(mu=as.numeric(coefs_orig[i,]),Sigma=covm[[i]])
  #     }
  #   }
  #   coefs_new=resultsdfcopy[,c(paste0("c_",1:E),"c_0")]
  #   resultsdfcopy$pred=rowSums(predictors*coefs_new)
  #   ss[j]=sum((resultsdfcopy$obs-resultsdfcopy$pred)^2,na.rm = T)
  #   modelresultscopy$resultsdf=resultsdfcopy
  # 
  #   jacobians=getJacobians(modelresultscopy)
  #   JLE[j]=LEglobal(modelresultscopy, jacobians)
  # }
  
  output=data.frame(JLE=JLE)
  output$W=exp(-n/2*log(ss/ss_opt))
  output$W=output$W/sum(output$W)
  output=output[order(output$JLE),]
  output$Wcum=cumsum(output$W)
  rownames(output)=NULL
  
  return(output)
}

# mvnormsamp=function(n,mu,Sigma) {
#   A=chol(Sigma)
#   r=numeric(n*length(mu)) 
#   z=matrix(rnorm(n=n*length(mu),mean = 0,sd = 1),ncol = n, nrow = length(mu))
#   f=t(A)%*%z+mu
#   return(f)
# }

#Residual bootstrap method
# LEbootrep=function(modelresults, nreps=100) {
#   LEvec=numeric(length = nreps)
#   for(i in 1:nreps) {
#     LEvec[i]=LEboot(modelresults)
#   }
#   return(LEvec)
# } 

LEboot=function(modelresults, nreps=100) {
  
  LEvec=numeric(length = nreps)
  
  n=nrow(modelresults$resultsdf)
  resid=modelresults$resultsdf$obs-modelresults$resultsdf$pred
  
  yboot=matrix(nrow = n, ncol = nreps)
  for(i in 1:nreps) {
    yboot[,i]=modelresults$resultsdf$pred+sample(x=na.omit(resid),size=n,replace=T)
  }
  
  #get E, tau, theta values from best model
  E=modelresults$modelstats$E
  tau=modelresults$modelstats$tau
  theta=modelresults$modelstats$theta
  hyperpars=data.frame(bestE=E, bestTau=tau, bestTheta=theta, 
                       Emax = modelresults$modelstats$Emax, taumax = modelresults$modelstats$taumax)
  #get best model form
  modelform=modelresults$form
  modnum=case_when(modelform =="ut-ut" ~ 1,
                   modelform =="log-log" ~ 2,
                   modelform =="fd-ut" ~ 3,
                   modelform =="gr-ut" ~ 4,
                   modelform =="gr-log" ~ 5)
  
  for(i in 1:nreps) {
    modelresultsboot=smap_model_options(data=modelresults$resultsdf, model=modnum, y="Obs_abund", 
                                        hypars=hyperpars, yboot=yboot[,i])
    jacobians=getJacobians(modelresultsboot)
    LEvec[i]=LEglobal(modelresultsboot, jacobians)
  }
  return(LEvec)
}

#This computes global LE
#The global LE is calculated for the whole time series.
LEglobal=function(modelresults, jacobians) {
  tau=modelresults$modelstats$tau
  len=dim(jacobians)[3]
  ndim=dim(jacobians)[1]
  
  #global le
  runlen=rle(!is.na(jacobians[1,1,]))
  serlen=max(runlen$lengths[runlen$values==TRUE])
  if(ndim==1) {
    #average over all timepoints
    gle=mean(log(abs(jacobians[1,1,])), na.rm=T)/tau
  } else {
    #compute over longest run of non-missing values
    starti=sum(runlen$lengths[1:(which.max(runlen$lengths)-1)])+1
    endi=sum(runlen$lengths[1:(which.max(runlen$lengths))])
    LEtemp=numeric(tau)
    for(i in 1:tau){
      indices=seq(from=starti+i-1, to=endi, by=tau)
      Jacs=jacobians[,,indices]
      nk=dim(Jacs)[3]
      Jk=Jacs[,,1]
      QR=qr(Jk)
      R=qr.R(QR)
      Q=qr.Q(QR)
      Rcum=R
      for(k in 2:nk) {
        Jk=Jacs[,,k]
        QR=qr(Jk%*%Q)
        R=qr.R(QR)
        Q=qr.Q(QR)
        Rcum=R%*%Rcum
        #print(diag(Rcum))
      }
      LEtemp[i]=log(max(abs(diag(Rcum))))/tau/nk
    }
    gle=mean(LEtemp)
  }
  return(gle)
}

#Asymptotic SE for LE
LEasymp=function(modelresults, jacobians) {
  
  tau=modelresults$modelstats$tau
  runlen=rle(!is.na(jacobians[1,1,]))

  #compute over longest run of non-missing values
  starti=sum(runlen$lengths[1:(which.max(runlen$lengths)-1)])+1
  endi=sum(runlen$lengths[1:(which.max(runlen$lengths))])
  LEtemp=numeric(tau)
  L=list()
  for(i in 1:tau){
    indices=seq(from=starti+i-1, to=endi, by=tau)
    Jacs=jacobians[,,indices,drop=F]
    nk=dim(Jacs)[3]
    Jk=Jacs[,,1]
    QR=qr(Jk)
    R=qr.R(QR)
    Q=qr.Q(QR)
    Rcum=R
    L[[i]]=numeric(nk)
    L[[i]][1]=log(max(abs(diag(Rcum))))
    for(k in 2:nk) {
      Jk=Jacs[,,k]
      QR=qr(Jk%*%Q)
      R=qr.R(QR)
      Q=qr.Q(QR)
      Rcum=R%*%Rcum
      L[[i]][k]=log(max(abs(diag(Rcum)))) 
    }
    LEtemp[i]=log(max(abs(diag(Rcum))))/tau/nk
  }
  gle=mean(LEtemp)
  
  #truncate to length of shortest series
  n=min(sapply(L,length))
  lags=(-n+1):(n-1)
  
  #Kernel to interpolate ACF - supposedly improves estimate of SE
  K=function(x) (1-abs(x))*(abs(x)<1)
  rho=0.9 #guessing autocorrelation like AR(1)
  a1=4*rho^2/((1-rho)^2*(1+rho)^2)
  Kt=1.14*(a1*n)^(1/3)
  w=K(lags/Kt) #weights to use
  
  V=numeric(tau)
  for(i in 1:length(L)) {
    Li=L[[i]][1:n]
    Nu=c(Li[1],diff(Li))-gle
    Q=outer(Nu,Nu) #outer products
    C=numeric(length(lags))
    for(j in lags) {
      C[n+j]=sum(mgcv::sdiag(Q,j)) #sum at each lag - to approximate ACF
    }
    C=C/n #approximate ACF
    V[i]=w%*%C/n
  }

  V=mean(V) #variance, one for each tau
  SE=sqrt(V)/tau
  LElower=gle-1.64*SE
  
  return(list(gle=gle, SE=SE, LElower=LElower))
}
