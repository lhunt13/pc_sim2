# functions


# nsim = number of simulations used to perform g-comp
# fits = the model fits you're simulating from
# K = number of time points in study
gcomp = function(nsim,fits,K){
  D_lag = matrix(NA,nrow=nsim,ncol=K)
  Y_lag = matrix(NA,nrow=nsim,ncol=K)
  B_lag = matrix(NA,nrow=nsim,ncol=K)
  Dsim = matrix(NA,nrow=nsim,ncol=K)
  Ysim = matrix(NA,nrow=nsim,ncol=K)
  Bsim = matrix(NA,nrow=nsim,ncol=K)
  
  D_lag[,1] = rep(0,nsim)
  Y_lag[,1] = rep(0,nsim)
  B_lag[,1] = rep(0,nsim)
  
  for(k in 1:K){
    intervention = data.frame(D=1,Y=1,B=1,k=k,B_lag=B_lag[,k],A_lag=1)
    
    # simulate D_k
    X = model.matrix(fits$formD, data=intervention)
    eta = apply(X,1,function(x){sum(x*fits$coefD)})
    probs = 1/(1 + exp(-eta))
    sim = rbinom(nsim,1,probs)
    Dsim[,k] = ifelse(D_lag[,k]==1,1,ifelse(Y_lag[,k]==1,0,sim))
    
    # simulate Y_k
    X = model.matrix(fits$formY, data=intervention)
    eta = apply(X,1,function(x){sum(x*fits$coefY)})
    probs = 1/(1 + exp(-eta))
    sim = rbinom(nsim,1,probs)
    Ysim[,k] = ifelse(Y_lag[,k]==1,1,ifelse(Dsim[,k]==1,0,sim))
    
    # simulate B_k
    X = model.matrix(fits$formB, data=intervention)
    eta = apply(X,1,function(x){sum(x*fits$coefB)})
    probs = 1/(1 + exp(-eta))
    sim = rbinom(nsim,1,probs)
    Bsim[,k] = ifelse(B_lag[,k]==0,sim,2)
    
    if(k < K){
      D_lag[,k+1] = Dsim[,k]
      Y_lag[,k+1] = Ysim[,k]
      B_lag[,k+1] = Bsim[,k]
    }
  }
  
  py = colMeans(Ysim)
  pd = colMeans(Dsim)
  names(py) = paste0("day",1:K)
  names(pd) = paste0("day",1:K)
  return(list(risk = py,competing_risk = pd))
}

# do same but treating competing risk as censoring event
gcomp2 = function(nsim,fits,K){
  Y_lag2 = matrix(NA,nrow=nsim,ncol=K)
  B_lag = matrix(NA,nrow=nsim,ncol=K)
  Ysim = matrix(NA,nrow=nsim,ncol=K)
  Bsim = matrix(NA,nrow=nsim,ncol=K)
  
  Y_lag2[,1] = rep(0,nsim)
  B_lag[,1] = rep(0,nsim)
  
  for(k in 1:K){
    intervention = data.frame(Y2=1,B=1,k=k,B_lag=B_lag[,k],A_lag=1)
    
    # simulate Y_k
    X = model.matrix(fits$formY, data=intervention)
    eta = apply(X,1,function(x){sum(x*fits$coefY)})
    probs = 1/(1 + exp(-eta))
    sim = rbinom(nsim,1,probs)
    Ysim[,k] = ifelse(Y_lag2[,k]==1,1,sim)
    
    # simulate B_k
    X = model.matrix(fits$formB, data=intervention)
    eta = apply(X,1,function(x){sum(x*fits$coefB)})
    probs = 1/(1 + exp(-eta))
    sim = rbinom(nsim,1,probs)
    Bsim[,k] = ifelse(B_lag[,k]==0,sim,2)
    
    if(k < K){
      Y_lag2[,k+1] = Ysim[,k]
      B_lag[,k+1] = Bsim[,k]
    }
  }
  
  py = colMeans(Ysim)
  names(py) = paste0("day",1:K)
  return(py)
}

# perform analysis with a dataset for A and L
do_analysis = function(data_L,data_A,nsim,K){
  # fit Lovenox models
  modD_L = speedglm(D ~ k + A_lag + I(B_lag > 0) + I(B_lag > 0):A_lag, 
                    family = binomial(), 
                    data = subset(data_L,D_lag==0 & Y_lag==0 & C_lag==0))
  modY_L = speedglm(Y ~ k + A_lag + I(B_lag > 0) + I(B_lag > 0):A_lag, 
                    family = binomial(), 
                    data = subset(data_L, D==0 & Y_lag==0 & C_lag==0))
  modB_L = speedglm(B ~ k + A_lag, 
                    family = binomial(), 
                    data = subset(data_L, B_lag==0 & C_lag==0 & Y==0 & D==0))
  modA_L = speedglm(A ~ k + A_lag + I(B > 0), 
                    family = binomial(), 
                    data = subset(data_L,B_lag==0 & C_lag==0 & Y==0 & D==0))
  modC_L = speedglm(C ~ k + I(B > 0), 
                    family = binomial(), 
                    data = subset(data_L, C_lag==0 & D==0 & Y==0))
  fits_L = list(coefD = coef(modD_L), coefY = coef(modY_L), coefB = coef(modB_L), 
                coefA = coef(modA_L), coefC = coef(modC_L),
                formD = formula(modD_L), formY = formula(modY_L), formB = formula(modB_L),
                formA = formula(modA_L), formC = formula(modC_L))
  
  # fit Aspirin models
  modD_A = speedglm(D ~ k + A_lag + I(B_lag > 0) + I(B_lag > 0):A_lag, 
                    family = binomial(), 
                    data = subset(data_A,D_lag==0 & Y_lag==0 & C_lag==0))
  modY_A = speedglm(Y ~ k + A_lag + I(B_lag > 0) + I(B_lag > 0):A_lag, 
                    family = binomial(), 
                    data = subset(data_A, D==0 & Y_lag==0 & C_lag==0))
  modB_A = speedglm(B ~ k + A_lag, 
                    family = binomial(), 
                    data = subset(data_A, B_lag==0 & C_lag==0 & Y==0 & D==0))
  modA_A = speedglm(A ~ k + A_lag + I(B > 0), 
                    family = binomial(), 
                    data = subset(data_A,B_lag==0 & C_lag==0 & Y==0 & D==0))
  modC_A = speedglm(C ~ k + I(B > 0), 
                    family = binomial(), 
                    data = subset(data_A, C_lag==0 & D==0 & Y==0))
  fits_A = list(coefD = coef(modD_A), coefY = coef(modY_A), coefB = coef(modB_A), 
                coefA = coef(modA_A), coefC = coef(modC_A),
                formD = formula(modD_A), formY = formula(modY_A), formB = formula(modB_A),
                formA = formula(modA_A), formC = formula(modC_A))
  
  
  # g-comp & ITT with competing risks
  result_good_L = gcomp(nsim,fits_L,90)
  result_good_A = gcomp(nsim,fits_A,90)
  setDT(data_L)
  Y_L = data.table::dcast(data_L, i ~ k, value.var = c('Y'))
  Y_L[Y_L==-1] = NA
  D_L = data.table::dcast(data_L, i ~ k, value.var = c('D'))
  D_L[D_L==-1] = NA
  setDT(data_A)
  Y_A = data.table::dcast(data_A, i ~ k, value.var = c('Y'))
  Y_A[Y_A==-1] = NA
  D_A = data.table::dcast(data_A, i ~ k, value.var = c('D'))
  D_A[D_A==-1] = NA
  
  gcomp_risk_area = sum(result_good_L$risk) - sum(result_good_A$risk)
  gcomp_comp_area = sum(result_good_L$competing_risk) - sum(result_good_A$competing_risk)
  gcomp_risk_last = result_good_L$risk[90] - result_good_A$risk[90]
  gcomp_comp_last = result_good_L$competing_risk[90] - result_good_A$competing_risk[90]
  itt_risk_area = sum(colMeans(Y_L,na.rm=T)[-1]) - sum(colMeans(Y_A,na.rm=T)[-1])
  itt_comp_area = sum(colMeans(D_L,na.rm=T)[-1]) - sum(colMeans(D_A,na.rm=T)[-1])
  itt_risk_last = colMeans(Y_L,na.rm=T)[91] - colMeans(Y_A,na.rm=T)[91]
  itt_comp_last = colMeans(D_L,na.rm=T)[91] - colMeans(D_A,na.rm=T)[91]
  
  ### g-comp & ITT where death is treated as censoring
  # fit Lovenox models
  data_L$C_lag2 = data_L$C_lag + 1*(data_L$D==1)
  data_L$C2 = data_L$C + 1*(data_L$D==1)
  data_L$Y2 = ifelse(data_L$D==1,-1,data_L$Y)
  data_L$Y_lag2 = ifelse(data_L$D_lag==1,-1,data_L$Y_lag)
  modY_Lb = speedglm(Y2 ~ k + A_lag + I(B_lag > 0) + I(B_lag > 0):A_lag, 
                     family = binomial(), 
                     data = subset(data_L, Y_lag2==0 & C_lag2==0))
  modB_Lb = speedglm(B ~ k + A_lag, 
                     family = binomial(), 
                     data = subset(data_L, B_lag==0 & C_lag2==0 & Y2==0))
  modA_Lb = speedglm(A ~ k + A_lag + I(B > 0), 
                     family = binomial(), 
                     data = subset(data_L,B_lag==0 & C_lag2==0 & Y2==0 ))
  modC_Lb = speedglm(C2 ~ k + I(B > 0), 
                     family = binomial(), 
                     data = subset(data_L, C_lag2==0 & Y2==0))
  fits_Lb = list(coefY = coef(modY_Lb), coefB = coef(modB_Lb), 
                 coefA = coef(modA_Lb), coefC = coef(modC_Lb), 
                 formY = formula(modY_Lb), formB = formula(modB_Lb),
                 formA = formula(modA_Lb), formC = formula(modC_Lb))
  
  # fit Aspirin models
  data_A$C_lag2 = data_A$C_lag + 1*(data_A$D==1)
  data_A$C2 = data_A$C + 1*(data_A$D==1)
  data_A$Y2 = ifelse(data_A$D==1,-1,data_A$Y)
  data_A$Y_lag2 = ifelse(data_A$D_lag==1,-1,data_A$Y_lag)
  modY_Ab = speedglm(Y2 ~ k + A_lag + I(B_lag > 0) + I(B_lag > 0):A_lag, 
                     family = binomial(), 
                     data = subset(data_A, Y_lag2==0 & C_lag2==0))
  modB_Ab = speedglm(B ~ k + A_lag, 
                     family = binomial(), 
                     data = subset(data_A, B_lag==0 & C_lag2==0 & Y2==0))
  modA_Ab = speedglm(A ~ k + A_lag + I(B > 0), 
                     family = binomial(), 
                     data = subset(data_A,B_lag==0 & C_lag2==0 & Y2==0 ))
  modC_Ab = speedglm(C2 ~ k + I(B > 0), 
                     family = binomial(), 
                     data = subset(data_A, C_lag2==0 & Y2==0))
  fits_Ab = list(coefY = coef(modY_Ab), coefB = coef(modB_Ab), 
                 coefA = coef(modA_Ab), coefC = coef(modC_Ab), 
                 formY = formula(modY_Ab), formB = formula(modB_Ab),
                 formA = formula(modA_Ab), formC = formula(modC_Ab))
  
  # get estimates
  result_bad_L = gcomp2(nsim,fits_Lb,90)
  result_bad_A = gcomp2(nsim,fits_Ab,90)
  Y_L2 = data.table::dcast(data_L, i ~ k, value.var = c('Y2'))
  Y_L2[Y_L2==-1] = NA
  Y_A2 = data.table::dcast(data_A, i ~ k, value.var = c('Y2'))
  Y_A2[Y_A2==-1] = NA
  
  gcomp_area2 = sum(result_bad_L) - sum(result_bad_A)
  gcomp_last2 = result_bad_L[90] - result_bad_A[90]
  itt_area2 = sum(colMeans(Y_L2,na.rm=T)[-1]) - sum(colMeans(Y_A2,na.rm=T)[-1])
  itt_last2 = colMeans(Y_L2,na.rm=T)[91] - colMeans(Y_A2,na.rm=T)[91]
  
  
  estimates = c(gcomp_risk_area,gcomp_comp_area,gcomp_risk_last,gcomp_comp_last,
                itt_risk_area,itt_comp_area,itt_risk_last,itt_comp_last,
                gcomp_area2,gcomp_last2,itt_area2,itt_last2)
  names(estimates) = c('gcomp_risk_area','gcomp_comp_area','gcomp_risk_last','gcomp_comp_last',
                       'itt_risk_area','itt_comp_area','itt_risk_last','itt_comp_last',
                       'gcomp_area2','gcomp_last2','itt_area2','itt_last2')
  
  return(list(estimates = estimates,fits_L=fits_L,fits_A=fits_A))
}


# n = sample size
# truth = distribution you're sampling from
# K = number of time points in study
sim_data = function(n,truth,K){
  Dsim = matrix(NA,nrow=n,ncol=K)
  Ysim = matrix(NA,nrow=n,ncol=K)
  Bsim = matrix(NA,nrow=n,ncol=K)
  Asim = matrix(NA,nrow=n,ncol=K)
  Csim = matrix(NA,nrow=n,ncol=K)
  D_lag = matrix(NA,nrow=n,ncol=K)
  Y_lag = matrix(NA,nrow=n,ncol=K)
  B_lag = matrix(NA,nrow=n,ncol=K)
  A_lag = matrix(NA,nrow=n,ncol=K)
  C_lag = matrix(NA,nrow=n,ncol=K)
  D_lag[,1] = 0
  Y_lag[,1] = 0
  B_lag[,1] = 0
  A_lag[,1] = 1
  C_lag[,1] = 0
  
  for(k in 1:K){
    newdata = data.frame(A_lag = A_lag[,k],
                         B_lag = B_lag[,k],
                         D = 1, Y = 1, B = 1, k=k)
    X = model.matrix(truth$formD,data=newdata)
    eta = apply(X,1,function(x){sum(x*truth$coefD)})
    probs = 1/(1 + exp(-eta))
    Dsim[,k] = ifelse(D_lag[,k]==1,1,
                      ifelse(Y_lag[,k]==1,0,
                             ifelse(C_lag[,k]==1,-1,
                                    rbinom(n,1,probs))))
    
    X = model.matrix(truth$formY,data=newdata)
    eta = apply(X,1,function(x){sum(x*truth$coefY)})
    probs = 1/(1 + exp(-eta))
    Ysim[,k] = ifelse(Y_lag[,k]==1,1,
                      ifelse(Dsim[,k]==1,0,
                             ifelse(C_lag[,k]==1,-1,
                                    rbinom(n,1,probs))))
    
    
    X = model.matrix(truth$formB,data=newdata)
    eta = apply(X,1,function(x){sum(x*truth$coefB)})
    probs = 1/(1 + exp(-eta))
    Bsim[,k] = ifelse(Dsim[,k] + Ysim[,k]==1,-1,
                      ifelse(C_lag[,k]==1,-1,
                             ifelse(B_lag[,k]==0,rbinom(n,1,probs),2)))
    
    
    newdata = data.frame(A_lag = A_lag[,k],
                         B_lag = B_lag[,k],
                         B = Bsim[,k],
                         A = 1, k=k)
    X = model.matrix(truth$formA,data=newdata)
    eta = apply(X,1,function(x){sum(x*truth$coefA)})
    probs = 1/(1 + exp(-eta))
    Asim[,k] = ifelse(Dsim[,k] + Ysim[,k]==1,-1,
                      ifelse(C_lag[,k]==1,-1,
                             ifelse(B_lag[,k]==1,A_lag[,k],rbinom(n,1,probs))))
    
    newdata = data.frame(A_lag = A_lag[,k],
                         B_lag = B_lag[,k],
                         B = Bsim[,k],
                         A = Asim[,k],
                         C = 1, k=k)
    X = model.matrix(truth$formC,data=newdata)
    eta = apply(X,1,function(x){sum(x*truth$coefC)})
    probs = 1/(1 + exp(-eta))
    Csim[,k] = ifelse(Dsim[,k] + Ysim[,k]==1,0,
                      ifelse(C_lag[,k]==1,1,rbinom(n,1,probs)))
    
    
    if(k < K){
      D_lag[,k+1] = Dsim[,k]
      Y_lag[,k+1] = Ysim[,k]
      B_lag[,k+1] = Bsim[,k]
      A_lag[,k+1] = Asim[,k]
      C_lag[,k+1] = Csim[,k]
    }
  }
  
  data = data.frame(i=rep(1:n,each=K),
                    k=rep(1:K,n),
                    D=c(t(Dsim)),
                    Y=c(t(Ysim)),
                    B=c(t(Bsim)),
                    A=c(t(Asim)),
                    C=c(t(Csim)),
                    D_lag=c(t(D_lag)),
                    Y_lag=c(t(Y_lag)),
                    B_lag=c(t(B_lag)),
                    A_lag=c(t(A_lag)),
                    C_lag=c(t(C_lag)))
  
  return(data)
}





