# Functions for GPDD stability analysis

#get max E, theta, R2, and s-map coeffs for each series

#true R2 function
getR2=function(resultsdf, obse="obs", prede="pred") {
  d=na.omit(select(resultsdf, obse, prede))
  R2=1-sum((d[,obse]-d[,prede])^2)/sum((d[,obse]-mean(d[,obse]))^2)
}

besthyper=function(data=NULL, y, pgr=NULL, Efix=NULL, taufix=NULL) { 
  #get best E, tau, theta
  #y is the time series
  #if providing dataframe containing y, give y as character (column name)
  if(!is.null(data)) {
    ser=data[,y]
  } else {
    ser=y
  }
  #set max E based on effective ts length (longest string of non-NA values)
  len=rle(!is.na(ser))
  serlen=max(len$lengths[len$values==TRUE])
  #Emax=min(12,round(sqrt(serlen),0))
  Emax=min(12,floor(sqrt(serlen))) #round down
  
  if(!is.null(taufix)) {
    tautry=taufix #set tau to fixed value
  } else {
    tautry=1:Emax
  }
  if(!is.null(Efix)) {
    Etry=Efix #set E to fixed value
  } else {
    Etry=1:Emax
  }
  simplex_results=expand.grid(Etry, tautry)
  colnames(simplex_results)=c("E","tau")
  simplex_results=filter(simplex_results, E*tau<=Emax)
  simplex_results$rmse=NA
  simplex_results$npred=NA
  
  #generate lags
  ser_lags=make_block(ser, max_lag = Emax+1, tau=1)[,-1]
  #if using growth rate as response, substitute it
  if(!is.null(pgr)) {
    ser_lags$col1=pgr
  }
  #ser_lags=na.omit(ser_lags) #standardize targets
  #get best E and tau using simplex
  for(i in 1:nrow(simplex_results)) {
    simptemp=block_lnlp(ser_lags, tp = 0, method = "simplex", 
                        columns = paste0("col1_",simplex_results$tau[i]*(1:simplex_results$E[i])), 
                        target_column = "col1",
                        first_column_time = F, silent = T, stats_only = F)
    simplex_results$rmse[i]=simptemp$rmse
    simplex_results$npred[i]=simptemp$num_pred
  }
  simplex_results$error=with(simplex_results, ifelse(npred-E^2<0, NA, npred*rmse^2/(npred-E^2)))
  bestE=simplex_results$E[which.min(simplex_results$error)]
  bestTau=simplex_results$tau[which.min(simplex_results$error)]
  #get best theta
  thetatest = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  ser_lags=make_block(ser, max_lag = bestE+1, tau=bestTau)[,-1]
  if(!is.null(pgr)) {
    ser_lags$col1=pgr
  }
  smap_results=block_lnlp(ser_lags, tp = 0, method = "s-map", 
                          columns = paste0("col1_",bestTau*(1:bestE)),
                          target_column = "col1", theta=thetatest,
                          first_column_time = F, silent = T)
  bestTheta=smap_results$theta[which.min(smap_results$rmse)]
  
  return(data.frame(bestE=bestE, bestTau=bestTau, bestTheta=bestTheta, Emax=Emax))
}

smap_model=function(data, y, ylog, gr=NULL, fd=NULL, Efix=NULL, taufix=NULL) {
  #can supply either gr (growth rate) or fd (first difference) or none
  #ylog indicates whether y is log transformed
  #if fd, must use ylog=F
  ser=data[,y]
  pgr=NULL
  if(!is.null(gr)) pgr=data[,gr]
  if(!is.null(fd)) pgr=data[,fd]
  #get hyper parameters
  hyperpars=besthyper(y=ser, pgr=pgr, Efix=Efix, taufix=taufix)
  #store smap output
  ser_lags=make_block(ser, max_lag = hyperpars$bestE+1, tau=hyperpars$bestTau)[,-1]
  if(!is.null(pgr)) {
    ser_lags$col1=pgr
  }
  smap_results=block_lnlp(ser_lags, tp = 0, method = "s-map", 
                          columns = paste0("col1_",hyperpars$bestTau*(1:hyperpars$bestE)), 
                          target_column = "col1", theta=hyperpars$bestTheta,
                          first_column_time = F, silent = T,
                          stats_only = F, save_smap_coefficients = T)
  resultsdf=cbind(smap_results$model_output[[1]], smap_results$smap_coefficients[[1]])
  #untransformed abundance
  resultsdf$Obs_abund=data$PopRescale
  #get untranformed abundance predictions based on model type
  if(is.null(fd) & is.null(gr) & ylog==F) {
    resultsdf$Pred_abund=resultsdf$pred
    form="ut-ut"
  }  
  if(is.null(fd) & is.null(gr) & ylog==T) {
    resultsdf$Pred_abund=exp(resultsdf$pred)
    form="log-log"
  }
  if(!is.null(fd)) {
    resultsdf$Pred_abund=lag(resultsdf$Obs_abund)+resultsdf$pred
    form="fd-ut"
  }
  if(!is.null(gr)) {
    resultsdf$Pred_abund=lag(resultsdf$Obs_abund)*exp(resultsdf$pred)
    #ylog is not used to get predictions, but used for jacobian construction
    if(ylog==T) form="gr-log"
    if(ylog==F) form="gr-ut"
  }
  #R2 for model
  modelR2=getR2(resultsdf) 
  #R2 for untransformed abundance
  modelR2_abund=getR2(resultsdf, obse="Obs_abund", prede="Pred_abund")
  modelstats=data.frame(E=hyperpars$bestE, tau=hyperpars$bestTau, theta=hyperpars$bestTheta, 
                        Emax=hyperpars$Emax, num_pred=smap_results$num_pred,
                        R2model=modelR2, R2abund=modelR2_abund, rho=smap_results$rho)
  return(list(modelstats=modelstats, resultsdf=resultsdf, form=form))
}

#select among multiple models
smap_model_options=function(data, model, Efix=NULL, taufix=NULL) {
  if(model==1) {
    modelresults=smap_model(data, y="PopRescale", ylog=F, Efix=Efix, taufix=taufix)
  }
  if(model==2) {
    modelresults=smap_model(data, y="PopRescale_log", ylog=T, Efix=Efix, taufix=taufix)
  }
  if(model==3) {
    modelresults=smap_model(data, fd="PopRescale_fd", y="PopRescale", ylog=F, Efix=Efix, taufix=taufix)
  }
  if(model==4) {
    modelresults=smap_model(data, fd="PopRescale_gr", y="PopRescale", ylog=F, Efix=Efix, taufix=taufix)
  }
  if(model==5) {
    modelresults=smap_model(data, fd="PopRescale_gr", y="PopRescale_log", ylog=T, Efix=Efix, taufix=taufix)
  }
  return(modelresults)
}

#construct Jacobians
getJacobians=function(modelresults){
  #this will vary depending on the model, data transform
  form=modelresults$form
  
  ndim=modelresults$modelstats$E #dimension
  len=nrow(modelresults$resultsdf)
  coefs=modelresults$resultsdf[,paste0("c_",1:ndim)]
  jacobians=array(dim = c(ndim,ndim,len))
  
  if(form=="ut-ut") {
    if(ndim==1) {
      jacobians[,,]=coefs
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1=matrix(nrow = ndim, ncol = ndim, 0)
          J1[1,]=as.numeric(coefs[i,])
          J1[2:ndim,1:(ndim-1)]=diag(ndim-1)
          jacobians[,,i]=J1
        }
      }
    }
  }
  
  if(form=="fd-ut") {
    if(ndim==1) {
      jacobians[,,]=coefs+1
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1=matrix(nrow = ndim, ncol = ndim, 0)
          J1[1,]=as.numeric(coefs[i,])+c(1,rep(0,ndim-1))
          J1[2:ndim,1:(ndim-1)]=diag(ndim-1)
          jacobians[,,i]=J1
        }
      }
    }
  }
  
  if(form=="log-log") {
    x_obs=modelresults$resultsdf$Obs_abund
    r_x=exp(modelresults$resultsdf$pred)
    if(ndim==1) {
      jacobians[,,]=(coefs)*r_x/lag(x_obs)
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1<-J2<-J3<-matrix(nrow = ndim, ncol = ndim, 0)
          diag(J1)<-c(r_x[i], x_obs[(i-1):(i-ndim+1)])
          #c(r_x[i], lag(x_obs)[i], lag(x_obs,2)[i], lag(x_obs,3)[i], lag(x_obs,4)[i])
          J2[1,]=as.numeric(coefs[i,])
          J2[2:ndim,1:(ndim-1)]=diag(ndim-1)
          diag(J3)<-1/x_obs[(i-1):(i-ndim)] 
          #c(1/(lag(x_obs)[i]), 1/(lag(x_obs,2)[i]), 1/(lag(x_obs,3)[i]), 1/(lag(x_obs,4)[i]), 1/(lag(x_obs,5)[i]))
          jacobians[,,i]=J1%*%J2%*%J3
        }
      }
    }
  }
  
  if(form=="gr-ut") {
    x_obs=modelresults$resultsdf$Obs_abund
    r_x=exp(modelresults$resultsdf$pred)
    if(ndim==1) {
      jacobians[,,]=(coefs*lag(x_obs)+1)*r_x
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1<-J2<-matrix(nrow = ndim, ncol = ndim, 0)
          diag(J1)<-c(r_x[i], rep(1,ndim-1))
          #c(r_x[i], 1, 1, 1)
          J2[1,]=as.numeric(coefs[i,])*x_obs[i-1]+c(1,rep(0,ndim-1))
          J2[2:ndim,1:(ndim-1)]=diag(ndim-1)
          jacobians[,,i]=J1%*%J2
        }
      }
    }
  }
  
  if(form=="gr-log") {
    x_obs=modelresults$resultsdf$Obs_abund
    r_x=exp(modelresults$resultsdf$pred)
    if(ndim==1) {
      jacobians[,,]=(coefs+1)*r_x
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1<-J2<-J3<-matrix(nrow = ndim, ncol = ndim, 0)
          diag(J1)<-c(r_x[i]*x_obs[i-1], x_obs[(i-1):(i-ndim+1)])
          #c(r_x[i]*lag(x_obs)[i], lag(x_obs)[i], lag(x_obs,2)[i], lag(x_obs,3)[i], lag(x_obs,4)[i])
          J2[1,]=as.numeric(coefs[i,])+c(1,rep(0,ndim-1))
          J2[2:ndim,1:(ndim-1)]=diag(ndim-1)
          diag(J3)<-1/x_obs[(i-1):(i-ndim)] 
          #c(1/(lag(x_obs)[i]), 1/(lag(x_obs,2)[i]), 1/(lag(x_obs,3)[i]), 1/(lag(x_obs,4)[i]), 1/(lag(x_obs,5)[i]))
          jacobians[,,i]=J1%*%J2%*%J3
        }
      }
    }
  }
  
  return(jacobians)
}

#get local and global LEs
getStability=function(jacobians) {
  len=dim(jacobians)[3]
  ndim=dim(jacobians)[1]
  
  if(ndim==1) {
    lle=log(abs(jacobians[1,1,]))
  } else {
    lle=numeric(length = len)
    for(i in 1:len) {
      Jac=jacobians[,,i]
      if(any(is.na(Jac))) {
        lle[i]=NA
      } else {
        lle[i]=log(max(abs(eigen(Jac, only.values = T)$values))) 
      }
    }
  }
  
  #proportion positive
  lle_pp=length(which(lle>0))/length(which(!is.na(lle)))
  
  #naive avg of local LEs
  lle_avg=mean(lle, na.rm=T)
  
  #global le
  runlen=rle(!is.na(jacobians[1,1,]))
  serlen=max(runlen$lengths[runlen$values==TRUE])
  if(serlen<(len/2.5) & serlen<25) {
    #do not compute if the longest run of non-missing values is less than 40% of ts length or less than 25
    gle=NA 
  } else {
    if(ndim==1) {
      #average over all timepoints
      gle=mean(log(abs(jacobians[1,1,])), na.rm=T)
    } else {
      #compute over longest run of non-missing values
      starti=sum(runlen$lengths[1:(which.max(runlen$lengths)-1)])+1
      endi=sum(runlen$lengths[1:(which.max(runlen$lengths))])
      Jacs=jacobians[,,starti:endi]
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
      }
      gle=1/nk*log(max(abs(diag(Rcum))))
    }
  }
  
  return(list(lle=lle, lle_pp=lle_pp, lle_avg=lle_avg, gle=gle))
}

#1d LE from best model
LE1d=function(data, bestmodel) {
  modelresults=smap_model_options(data, bestmodel, Efix=1)
  jacobians=getJacobians(modelresults)
  stability=getStability(jacobians)
  return(stability)
}

# Regression-based LE (Rosenstein method)
regLE=function(data, y) {
  #recompute bestE with simplex
  ser=data[,y]
  len=rle(!is.na(ser))
  serlen=max(len$lengths[len$values==TRUE])
  Emax=min(10,round(sqrt(serlen),0))
  simplex_results=simplex(ser, tau=1, tp=1, E=1:Emax, silent = T)
  bestE=simplex_results$E[which.min(simplex_results$rmse)]
  
  xi=make_block(ser, max_lag = bestE, tau=1)[-(1:(bestE-1)),-1] #generate E lags
  if(bestE==1) {ni=length(xi)} else {ni=nrow(xi)}
  navals=which(!complete.cases(xi))
  rownames(xi)=NULL
  if(bestE==1) {xi[navals]=NA} else {xi[navals,]=NA}
  
  
  #distij=as.matrix(dist(xi))+10000*diag(ni) 
  distij=as.matrix(dist(xi)) #distance matrix
  exclusion=1 #exclusion radius
  delta <- row(distij) - col(distij)
  distij[delta <= exclusion & delta >= -1*exclusion] <- distij[delta <= exclusion & delta >= -1*exclusion]+10000
  
  md=apply(distij,1,min,na.rm=T) #min distance vector
  md[navals]=NA
  minind=apply(distij,1,which.min) #nearest point index vector
  minind[navals]=NA
  smax=4 #max steps into future
  avg_dist=numeric(smax+1)
  md[md==0]<-NA #remove distances of 0 (otherwise fails when taking logs)
  avg_dist[1]=mean(log(md), na.rm=T) 
  for(s in 1:smax) {
    ind1=(1:ni)+s
    ind2=minind+s
    rows=which(ind1<=ni & ind2<=ni)
    if(bestE==1){
      ds=abs(xi[ind1[rows]]-xi[ind2[rows]])
    } else {
      ds=sqrt(rowSums((xi[ind1[rows],]-xi[ind2[rows],])^2))
    }
    ds[ds==0]<-NA #remove distances of 0 (otherwise fails when taking logs)
    avg_dist[s+1]=mean(log(ds), na.rm=T)
  }
  
  #regression
  xr=matrix(ncol=2, nrow=smax+1, data = c(rep(1,smax+1), 0:smax))
  # #drop initial distance
  # xr=matrix(ncol=2, nrow=smax, data = c(rep(1,smax), 1:smax))
  # avg_dist=avg_dist[2:(smax+1)]
  bb=solve(t(xr)%*%xr)%*%t(xr)%*%avg_dist
  verr=t(avg_dist-xr%*%bb)%*%(avg_dist-xr%*%bb)/(smax-1)
  vcoef=verr[1,1]*solve(t(xr)%*%xr)
  # reg=lm(avg_dist~c(0:smax))
  # summary(reg)
  # plot(0:smax, avg_dist)
  
  LEreg=bb[2]
  LEreg_se=sqrt(vcoef[2,2])
  
  return(data.frame(LEreg=LEreg, LEreg_se=LEreg_se))
}

#converts timesteps to months
timescale_mo=function(SamplingInterval, y) {
  ifelse(SamplingInterval=="annual", y*12,
         ifelse(SamplingInterval=="seasonal", y*6,
                ifelse(SamplingInterval=="bimonthly", y*2,
                       ifelse(SamplingInterval=="monthly", y,
                              ifelse(SamplingInterval=="4-week", y/1.07,
                                     ifelse(SamplingInterval=="weekly", y/4.286,NA))))))
  
}
