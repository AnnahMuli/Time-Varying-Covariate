############################################################################################
#  Joint Model: Weibull baseline hazard and gamma frailty                                  #
#  h0 = lambda*rho*t^(rh0-1)   and   H0 = lambda*t^(rho)                                   #
############################################################################################
#                                                                                          #
#  The function arguments are:                                                             #
#   - data     : a data.frame containing all the variables;                                #
#   - cluster  : the name of the variable in data containing cluster IDs;                  #
#   - id       : the name of the variable in data containing individual IDs;               #
#   - Cov      : the longitudinal covariate measurements                                   #
#   - CovTime  : time when the covariate measurements were taken                           #
#   - status   : censoring indicator                                                       #
#   - SurvTime : Survival time                                                             #
#   - A0       : Age at entry                                                              #
#   - inip     : initial values for the parameters (Lambda,Rho,alpha,theta,sigma_u,sigma_e)#
#   - hess     : hessian matrix to generate standard error                                 #
#   - randslope: include id level random slope in the longitudinal model                   #
############################################################################################
#                                                                                          #
#   Last modification : March  27, 2025                                                    #
#                                                                                          #
############################################################################################

fit_JM<- function(data, cluster, id, inip, Cov, CovTime, status, covariate,SurvTime,A0, hess, randslope) 
{
  library(survival); library(xtable); library(dplyr); library(numDeriv);library(lme4) ; library(doBy) ; library(cubature)
  
  Dat <-data
  Dat$id <- Dat[[deparse(substitute(id))]] 
  Dat$start <- Dat[[deparse(substitute(CovTime))]] 
  Dat$zij <-Dat[[deparse(substitute(covariate))]] 
  Dat$status <- Dat[[deparse(substitute(status))]] 
  Dat$cluster <- Dat[[deparse(substitute(cluster))]]
  Dat$ATF <- Dat[[deparse(substitute(SurvTime))]]
  Dat$A0 <- Dat[[deparse(substitute(A0))]]
  
  Dat<- Dat[Dat$A0<Dat$ATF,] # Only include individual who survived beyond the age at entry
  
  if (randslope == F){
    model = lmer(zij ~ start+(1|id)+(1|cluster), data=Dat)
    model
    #-------------------------------------------------------
    #        PARAMETERS ESTIMATES
    #-------------------------------------------------------
    Beta0<-coef(summary(model) )[1,1]
    Beta1<-coef(summary(model) )[2,1]
    
    #--------------------------------------------------------
    bij<-ranef(model)$id
    B_ij<-cbind(id =rownames(bij),bij)
    names(B_ij)[2] <- "bij"
    Merged_data<-merge(Dat, B_ij, by="id")
    head(Merged_data)
    
    
    ui<-ranef(model)$cluster
    Ui_cluster<-cbind(cluster=rownames(ui),ui)
    names(Ui_cluster)[2] <- "ui"
    Merged_data1<-merge(Merged_data,Ui_cluster, by="cluster")
    head(Merged_data1)
    dim(Merged_data1)
    LongData11<-Merged_data1[Merged_data1$ATF>Merged_data1$start,]
    LongData<-subset(LongData11, ave(id, cluster, FUN = function(x) length(unique(x))) > 1)
    
    dim(LongData)
    
    SurvData1=subset(LongData, !duplicated(id) )
    SurvData<-subset(SurvData1, ave(id, cluster, FUN = function(x) length(unique(x))) > 1)
    dim(SurvData)
    ## ---------------------------------------------------------------------
    ## ---------------------------------------------------------------------
    ## Weibull baseline hazard
    ## h0 = lambda*rho*t^(rh0-1)   and   H0 = lambda*t^(rho)
    ## ---------------------------------------------------------------------
    likelihood.Joint<-function(p)
    { 
      lambda       <-  exp(p[4] )
      rho      <- exp( p[6] )
      theta <- exp(p[2] ) # s=log(var.eps) <=> exp(s)=var.eps
      alpha  <- p[3]
      sigma_u <- exp( p[1]   )
      sigma_e <-  exp( p[5]  )
      n <- length(unique(LongData$cluster))
      L.ni <-  as.vector(table(LongData$cluster))
      yit<- LongData$zij-Beta0-LongData$bij-Beta1*LongData$start 
      L.yit<-split(yit, LongData$cluster )
      S.status<- aggregate( SurvData$status  , by =  list(SurvData$cluster)  ,FUN=sum )[2] 
      hit<- ( (lambda*rho*SurvData$ATF^(rho-1) )* exp(  alpha*(Beta0+SurvData$bij+Beta1*SurvData$ATF )  ) )^SurvData$status
      S.hit<-split(hit, SurvData$cluster )
      S.ni<-as.vector(table(SurvData$cluster))
      Hit<- unlist (sapply(1:nrow(SurvData),   function(i)  adaptIntegrate( function(s)
        exp(alpha*(Beta0 + SurvData$bij[i] ) )* s^(rho-1)* exp(alpha*Beta1*s) ,lowerLimit= c(0) , upperLimit= c( SurvData$ATF[i]) )[1]),
        use.names=FALSE)
      S.Hit<-split(Hit, SurvData$cluster )
      H0it<- unlist (sapply(1:nrow(SurvData),   function(i)  adaptIntegrate( function(s)
        exp(alpha*(Beta0 + SurvData$bij[i] ) )* s^(rho-1)* exp(alpha*Beta1*s) ,lowerLimit= c(0) , upperLimit= c( SurvData$A0[i]) )[1]),
        use.names=FALSE)
      S0.Hit<-split(H0it, SurvData$cluster )
      
      lik<-sum(log( unlist( 
        sapply(
          1:n, function(i)
            adaptIntegrate( function(uv)
              as.numeric(      ( (2*pi*(sigma_e^2 ) )^(-L.ni[i]/2) * exp(   -1/(2*(sigma_e^2 )  ) *  sum( (L.yit[[i]] - uv[1])^2  )  ) )*
                                 ( prod( S.hit[[i]] )*( exp(alpha*uv[1]) )^S.status[i,]    )*
                                 (( ( 1/theta+ sum(S0.Hit[[i]])*lambda*rho*exp(alpha*uv[1]) )^(  1/theta  ) )/ gamma(1/theta) )*
                                 ( gamma( S.status[i,] + 1/theta ) /( ( 1/theta+ sum(S.Hit[[i]])*lambda*rho*exp(alpha*uv[1]) )^( S.status[i,] + 1/theta  ) ) ) )* 
                dnorm(uv[1],mean=0,sd=sigma_u ) , 
              lowerLimit= c(-10 ) , upperLimit= c(10) )[1]
        )                                        
      ) ) )
      -lik
    }
    
    
    system.time(t2<- nlminb(inip, likelihood.Joint,  hessian = hess,
                            control = list(trace=T,  rel.tol = 1e-8, step.min= 1e-3)) )
    
    output <- c(sigma_u_JM = exp(t2$par[1]), Theta_JM=exp(t2$par[2]),
                alpha_JM= t2$par[3], Lambda_JM=exp(t2$par[4]),
                sigma_e_JM=exp(t2$par[5]),Rho_JM= exp(t2$par[6]) )
    
    if (hess == T){
      std3 <- sqrt(diag(solve(hessian(likelihood.Joint, t2$par))))
      
      Output2<-round( matrix( c( t2$par[3], std3[3], t2$par[2], std3[2] ,
                                 exp(t2$par[4]), exp(t2$par[4])*std3[4],
                                 exp(t2$par[6]), exp(t2$par[6])*std3[6],
                                 exp(t2$par[1]), exp(t2$par[1])*std3[1] ,
                                 exp(t2$par[5]),  exp(t2$par[5])*std3[5] ),
                              nrow=6, byrow=T),3)
      colnames(Output2) <- c("Estimate","SE")
      rownames(Output2) <- c("alpha","theta", "lambda","rho","Sigma_u","Sigma_e")
      return(Output2)
    }  else { output  } 
    
    }
  
  else {
    if (randslope == T) {
      
      model<-  lmer(zij ~ start+(1|id)+(0+start|id)+(1|cluster), data=Dat,REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead") )
     
      #-------------------------------------------------------
      #        PARAMETERS ESTIMATES
      #-------------------------------------------------------
      Beta0<-coef(summary(model) )[1,1] 
      Beta1<-coef(summary(model) )[2,1]
      
      #--------------------------------------------------------
      bij<-ranef(model)$id
      bij0<- bij[[1]] 
      bij1<-bij[[2]] 
     
      B_ij<-   cbind(id = rownames(bij), bij0, bij1)  
      names(B_ij)[2] <- "bij0"
      names(B_ij)[3] <- "bij1"
      
      Merged_data<-merge(Dat, B_ij, by="id")
      
      
      ui<-ranef(model)$cluster
      Ui_cluster<-cbind(cluster=rownames(ui),ui)
      names(Ui_cluster)[2] <- "ui"
      Merged_data1<-merge(Merged_data,Ui_cluster, by="cluster")
  
      
      LongData11<-Merged_data1[Merged_data1$ATF>Merged_data1$start,]
      LongData<-subset(LongData11, ave(id, cluster, FUN = function(x) length(unique(x))) > 1) 
      LongData$bij0<-as.numeric(LongData$bij0)
      LongData$bij1<-as.numeric(LongData$bij1)
   
      SurvData1=subset(LongData, !duplicated(id) )
      SurvData<-subset(SurvData1, ave(id, cluster, FUN = function(x) length(unique(x))) > 1 )
      
      ## ---------------------------------------------------------------------
      ## ---------------------------------------------------------------------
      ## Weibull baseline hazard
      ## h0 = lambda*rho*t^(rh0-1)   and   H0 = lambda*t^(rho)
      ## ---------------------------------------------------------------------
      
      
      likelihood.Joint<-function(p)
      { 
        lambda       <-  exp(p[4] )
        rho      <- exp( p[6] )
        theta <- exp(p[2]  ) # s=log(var.eps) <=> exp(s)=var.eps
        alpha  <- p[3]
        sigma_u <- exp( p[1]    )
        sigma_e <-  exp( p[5]   )
        n <- length(unique(LongData$cluster))
        L.ni <-  as.vector(table(LongData$cluster))
        yit<- LongData$zij-Beta0-LongData$bij0-LongData$bij1*LongData$start-Beta1*LongData$start 
        L.yit<-split(yit, LongData$cluster )
        #uuii <-1 # uv[1] random effect
        #vi_exp_ui<- 1 #uv[2]*exp(alpha*uv[1]) #__i.e_____vi*exp(alpha*ui)
        S.status<- aggregate( SurvData$status  , by =  list(SurvData$cluster)  ,FUN=sum )[2] 
        hit<- ( (lambda*rho*SurvData$ATF^(rho-1) )* exp(  alpha*(Beta0+SurvData$bij0+SurvData$bij1*SurvData$ATF+Beta1*SurvData$ATF )  ) )^SurvData$status
        S.hit<-split(hit, SurvData$cluster )
        S.ni<-as.vector(table(SurvData$cluster))
        Hit<- unlist (sapply(1:nrow(SurvData),   function(i)  adaptIntegrate( function(s)
          exp(alpha*(Beta0 + SurvData$bij0[i] ) )* s^(rho-1)* exp(alpha *(Beta1+ bij1[i] )*s     ) ,lowerLimit= c(0) , upperLimit= c( SurvData$ATF[i]) )[1]),
          use.names=FALSE)
        S.Hit<-split(Hit, SurvData$cluster )
        H0it<- unlist (sapply(1:nrow(SurvData),   function(i)  adaptIntegrate( function(s)
          exp(alpha*(Beta0 + SurvData$bij0[i] ) )* s^(rho-1)* exp(alpha *(Beta1+ bij1[i] )*s     ) ,lowerLimit= c(0) , upperLimit= c( SurvData$A0[i]) )[1]),
          use.names=FALSE)
        S0.Hit<-split(H0it, SurvData$cluster )
        
        c<- 1e-6 # To help the likelihood not end up as infinity in some cases
        
        lik<-sum(log( c + unlist( 
          sapply(
            1:n, function(i)
              adaptIntegrate( function(uv)
                as.numeric(      ( (2*pi*(sigma_e^2 ) )^(-L.ni[i]/2) * exp(   -1/(2*(sigma_e^2 )  ) *  sum( (L.yit[[i]] - uv[1])^2  )  ) )*
                                   ( prod( S.hit[[i]] )*( exp(alpha*uv[1]) )^S.status[i,]    )*
                                   (( ( 1/theta+ sum(S0.Hit[[i]])*lambda*rho*exp(alpha*uv[1]) )^(  1/theta  ) )/ gamma(1/theta) )*
                                   ( gamma( S.status[i,] + 1/theta ) /( ( 1/theta+ sum(S.Hit[[i]])*lambda*rho*exp(alpha*uv[1]) )^( S.status[i,] + 1/theta  ) ) ) )* 
                  dnorm(uv[1],mean=0,sd=sigma_u ) , 
                lowerLimit= c(-10 ) , upperLimit= c(10) )[1]
          )                                        
        ) ) )
        -lik
      }

      system.time(t2<- nlminb(inip, likelihood.Joint,  hessian = hess,
                              control = list(trace=T,  rel.tol = 1e-6, step.min= 1e-3)) )
      
      output <- c(sigma_u_JM = exp(t2$par[1]), Theta_JM=exp(t2$par[2]),
                  alpha_JM= t2$par[3], Lambda_JM=exp(t2$par[4]),
                  sigma_e_JM=exp(t2$par[5]),Rho_JM= exp(t2$par[6]) )
      
      if (hess == T){
        std3 <- sqrt(diag(solve(hessian(likelihood.Joint, t2$par))))
        
        Output2<-round( matrix( c( t2$par[3], std3[3], t2$par[2], std3[2] ,
                                   exp(t2$par[4]), exp(t2$par[4])*std3[4],
                                   exp(t2$par[6]), exp(t2$par[6])*std3[6],
                                   exp(t2$par[1]), exp(t2$par[1])*std3[1] ,
                                   exp(t2$par[5]),  exp(t2$par[5])*std3[5] ),
                                nrow=6, byrow=T),3)
        colnames(Output2) <- c("Estimate","SE")
        rownames(Output2) <- c("alpha","theta", "lambda","rho","Sigma_u","Sigma_e")
        return(Output2)
      }  else { output  }
    }       }
  
}


DatH<-read.csv( file="C:/Users/staammu/Desktop/SimulData.csv" )
attach(DatH)
fit_JM( data=DatH, cluster = pair, id=id, inip=c(  log(0.1), log(0.5), 2, log(0.05), log(0.05), log(2) ) , Cov=zij,  CovTime=start, status=status,SurvTime=ATF, A0=A0, covariate=zij, hess=F, randslope=F)

