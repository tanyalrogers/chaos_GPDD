# Functions for GPDD stability analysis

#get max E, theta, R2, and s-map coeffs for each series

#true R2 function
getR2=function(resultsdf, obse="obs", prede="pred") {
  d=na.omit(select(resultsdf, obse, prede))
  R2=1-sum((d[,obse]-d[,prede])^2)/sum((d[,obse]-mean(d[,obse]))^2)
}

#predictors same as response variable
smap_model=function(data, y, ylog, Efix=NULL) {
  ser=data[,y]
  #set max E based on effective ts length (longest string of non-NA values)
  len=rle(!is.na(ser))
  serlen=max(len$lengths[len$values==TRUE])
  Emax=min(10,round(sqrt(serlen),0))
  if(!is.null(Efix)) {
    bestE=Efix #set E to fixed value
  } else {
    #get best E using simplex
    simplex_results=simplex(ser, tau=1, tp=1, E=1:Emax, silent = T)
    bestE=simplex_results$E[which.min(simplex_results$rmse)]
  }
  #get best theta
  smap_results=s_map(ser, tau=1, tp=1, E=bestE, silent = T)
  bestTheta=smap_results$theta[which.min(smap_results$rmse)]
  #store smap output
  smap_results=s_map(ser, tau=1, tp=1, E=bestE, theta=bestTheta, 
                     stats_only = F, save_smap_coefficients = T, silent = T)
  resultsdf=cbind(smap_results$model_output[[1]], smap_results$smap_coefficients[[1]])
  #shift df forward so will fit timepoint 1
  resultsdf[2:nrow(resultsdf),]=resultsdf[1:(nrow(resultsdf)-1),]
  resultsdf[1,]=NA
  resultsdf[1,1]=1
  #get true loo R2
  modelR2=getR2(resultsdf)
  if(ylog==T) {
    resultsdf$Obs_abund=exp(ser) #exp(resultsdf$obs)
    resultsdf$Pred_abund=exp(resultsdf$pred)
    modelR2_abund=getR2(resultsdf, obse="Obs_abund", prede="Pred_abund")
    form="log-log"
  }
  if(ylog==F) {
    modelR2_abund=modelR2
    form="ut-ut"
  }
  modelstats=data.frame(E=bestE, Emax=Emax, theta=bestTheta, num_pred=smap_results$num_pred,
                        R2model=modelR2, R2abund=modelR2_abund, rho=smap_results$rho)
  return(list(modelstats=modelstats, resultsdf=resultsdf, form=form))
}

#growth rate or first diff as response variable
smap_model_gr_fd=function(data, gr=NULL, fd=NULL, y, ylog, Efix=NULL) {
  #supply either gr (growth rate) or fd (first difference)
  #if fd, must use ylog=F
  ser=data[,y]
  if(!is.null(gr)) pgr=data[,gr]
  if(!is.null(fd)) pgr=data[,fd]
  #set max E based on effective ts length (longest string of non-NA values)
  len=rle(!is.na(ser))
  serlen=max(len$lengths[len$values==TRUE])
  Emax=min(10,round(sqrt(serlen),0))
  #generate lags
  ser_lags=make_block(ser, max_lag = Emax+1, tau=1)[,-1]
  ser_lags=cbind(ser_lags, pgr)
  if(!is.null(Efix)) {
    bestE=Efix #set E to fixed value
  } else {
    #get best E using simplex
    simplex_results=numeric(length = Emax)
    for(i in 1:Emax) {
      simplex_results[i]=block_lnlp(ser_lags, tp = 0, method = "simplex", 
                                    columns = paste0("col1_",1:i), 
                                    target_column = "pgr", 
                                    first_column_time = F, silent = T)$rmse
    }
    bestE=which.min(simplex_results)
  }

  #get best theta
  thetatest = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  smap_results=block_lnlp(ser_lags, tp = 0, method = "s-map", 
                          columns = paste0("col1_",1:bestE), 
                          target_column = "pgr", theta=thetatest,
                          first_column_time = F, silent = T)
  bestTheta=smap_results$theta[which.min(smap_results$rmse)]
  #store smap output
  smap_results=block_lnlp(ser_lags, tp = 0, method = "s-map", 
                          columns = paste0("col1_",1:bestE), 
                          target_column = "pgr", theta=bestTheta,
                          first_column_time = F, silent = T,
                          stats_only = F, save_smap_coefficients = T)
  resultsdf=cbind(smap_results$model_output[[1]], smap_results$smap_coefficients[[1]])
  if(!is.null(gr)) {
    resultsdf$Obs_abund=data$PopRescale #exp(ser)
    resultsdf$Pred_abund=lag(resultsdf$Obs_abund)*exp(resultsdf$pred)
    #ylog not used to get predictions, but used for jacobian
    if(ylog==T) form="gr-log"
    if(ylog==F) form="gr-ut"
  }
  if(!is.null(fd)) {
    resultsdf$Obs_abund=data$PopRescale #exp(ser)
    resultsdf$Pred_abund=lag(resultsdf$Obs_abund)+resultsdf$pred
    form="fd-ut"
  }
  modelR2=getR2(resultsdf) #get true loo R2
  modelR2_abund=getR2(resultsdf, obse="Obs_abund", prede="Pred_abund")
  modelstats=data.frame(E=bestE, Emax=Emax, theta=bestTheta, num_pred=smap_results$num_pred,
                        R2model=modelR2, R2abund=modelR2_abund, rho=smap_results$rho)
  return(list(modelstats=modelstats, resultsdf=resultsdf, form=form))
}

#select among multiple models
smap_model_options=function(data, model, Efix=NULL) {
  if(model==1) {
    modelresults=smap_model(data, y="PopRescale", ylog=F, Efix=Efix)
  }
  if(model==2) {
    modelresults=smap_model(data, y="PopRescale_log", ylog=T, Efix=Efix)
  }
  if(model==3) {
    modelresults=smap_model_gr_fd(data, fd="PopRescale_fd", y="PopRescale", ylog=F, Efix=Efix)
  }
  if(model==4) {
    modelresults=smap_model_gr_fd(data, fd="PopRescale_gr", y="PopRescale", ylog=F, Efix=Efix)
  }
  if(model==5) {
    modelresults=smap_model_gr_fd(data, fd="PopRescale_gr", y="PopRescale_log", ylog=T, Efix=Efix)
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
