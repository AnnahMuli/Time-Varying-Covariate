
#################################################################################
#  RiskSet regression calibration: Weibull baseline hazard and gamma frailty    #
#  h0 = lambda*rho*t^(rh0-1)   and   H0 = lambda*t^(rho)                        #
#################################################################################
#                                                                               #
#  The function arguments are:                                                  #
#   - data     : a data.frame containing all the variables;                     #
#   - cluster  : the name of the variable in data containing cluster IDs;       #
#   - id       : the name of the variable in data containing individual IDs;    #
#   - Cov      : the longitudinal covariate measurements
#   - CovTime  : time when the covariate measurements were taken                #
#   - status   : censoring indicator                                            #
#   - SurvTime : Survival time                                                  #
#   - A0       : Age at entry                                                   #
#   - inip     : initial values for the parameters (Lambda,Rho,alpha,theta)     #
#   - hess     : hessian matrix to generate standard error                      #
#   - randslope: include id level random slope in the longitudinal model        #
#################################################################################
#                                                                               #
#   Last modification : March  27, 2025                                         #
#                                                                               #
#################################################################################

fit_RRC <- function(data, cluster, id, inip, Cov, CovTime, status, covariate,SurvTime,A0, hess, randslope) 
{
  library(survival); library(xtable); library(dplyr); library(numDeriv);library(lme4) ; library(doBy)
  
  Dat <-data
  Dat$id <- Dat[[deparse(substitute(id))]] 
  Dat$start <- Dat[[deparse(substitute(CovTime))]] 
  Dat$zij <-Dat[[deparse(substitute(covariate))]] 
  Dat$status <- Dat[[deparse(substitute(status))]] 
  Dat$cluster <- Dat[[deparse(substitute(cluster))]]
  Dat$ATF <- Dat[[deparse(substitute(SurvTime))]]
  Dat$A0 <- Dat[[deparse(substitute(A0))]]
  
  
  Dat<- Dat[Dat$A0<Dat$ATF,] # Only include individual who survived beyond the age at entry
  
  Dat<- Dat[order(Dat$ATF), ]
  tt <- sort(unique( c(Dat$A0, Dat$ATF[status==1] ) )) # Unique event times or age at entry.
  nt <- length(tt)

  if (randslope == F){
    predicted_Cov <- function (atf) {
      ds_subset<-  Dat[Dat$ATF >= atf & Dat$A0<=atf & Dat$start<=atf, ]
      fit =  lmer(zij ~ start + (1|id) + (1|cluster), data=ds_subset) 
      newSubset <- data.frame(id = unique(ds_subset$id), cluster =subset( ds_subset,!duplicated(id) )$cluster,
                           
                              A0 =subset(ds_subset, !duplicated(id) )$A0,
                              ATF =subset(ds_subset, !duplicated(id) )$ATF,
                              status =subset(ds_subset, !duplicated(id) )$status, start = atf )
      newSubset$predicted_Cov <- predict(fit, newdata = newSubset)
      newSubset
      
    }
    
    LMM_Output<-tt[1:(nt-1)]|>purrr::map_dfr(predicted_Cov) 
    
    LMM_Output_Merged<- LMM_Output
    Dat2<- LMM_Output_Merged[LMM_Output_Merged$ATF>LMM_Output_Merged$start,]
    #Arrange data in the start stop format. The idea is to create start-stop intervals where the stop of an interval is the start of the next interval.
    #The last stop is equal to when an individual experiences the event/is censored. The first start is has BMD measurements predicted at A0 or measurements taken at entry for individuals without a past.
    Dat3<-Dat2 %>% group_by(id) %>%
      mutate(stop = lead(start, default = 1e6),
             stop = ifelse(stop > ATF, ATF, stop), .after = start)%>%
      filter(stop<= ATF) %>%
      ungroup()
    Dat3<- orderBy(~id, data=Dat3)
    Dat3$status[Dat3$stop<Dat3$ATF] <-0 # censor all intervals with exception of the last interval for an individual. The status of the last interval will be either censored or observed.
    dat11 <- Dat3[Dat3$stop>Dat3$start,] 
    #This next step is for obtaining only fully observed clusters (i.e only twin pairs without singletons). For the output below it will be based on fully observed clusters.
    # Fully observed clusters.
    dat1<-subset(dat11, ave(id, cluster, FUN = function(x) length(unique(x))) > 1)  
    
    #The following will be required in computation of the likelihood.
    Di<-aggregate(dat1$status, by=list ( dat1$cluster ) ,FUN=sum) #Total  number of NOT censored observations per cluster
    e<-colSums(Di[ 2 ] )  #e gives the sum of uncensored observation in the whole dataset.
    di<-Di [ 2 ] #1st column of Di is pair identifier so needs to be omitted. 
    GG<-length(unique( dat1$cluster) ) # Number of clusters
    Output3<-round(   matrix( c( e, GG  ) , nrow=2, byrow=T),5)
    
    #Likelihood estimation. This is the implementation of the proposed method.
    Lik_RRC<-function(p)
    {
      ###exp(p[1])<-lambda ; exp( p[2] )<-rho ;  p[3]<-theta ; p[4]<-alpha
      cumhazStart<-  as.numeric( aggregate(out ~ cluster, transform( dat1,
                                                                     out = exp(p[1])*(start^exp(p[2]))*exp( predicted_Cov*p[4] ) ), FUN = sum)[,2] )
      cumhazEntry<-as.numeric( aggregate(out ~ cluster, transform(subset(dat1, !duplicated(id)), out = exp(p[1])*(start^exp(p[2]) )*exp( predicted_Cov*p[4]) ), FUN = sum)[,2] )
      cumhazStop<- as.numeric( aggregate(out ~ cluster, transform( dat1,
                                                                   out =exp(p[1])*(stop^exp(p[2]))*exp(predicted_Cov*p[4] ) ), FUN = sum)[,2] )
      lnhaz<- as.numeric( aggregate(out ~ cluster, transform( dat1,
                                                              out = status*( predicted_Cov*p[4] + log(exp(p[1])*exp(p[2])*stop^(exp(p[2])-1) ) ) ), FUN = sum)[,2] )
      ## Likelihood Equation  ###
      lik<-e*log( (p[3] ) )- #***1st term***
        GG*log( gamma(1/(p[3] )) )+ #***2nd term***
        sum(log(gamma(di +1/(p[3] ))) )+#***3rd term***
        sum(lnhaz)+ ##***4th term***
        sum( (1/(p[3]) )*log(1+cumhazEntry*(p[3] )))- # *** 5 term***
        sum( ( di   +1/(p[3] ) )*log(1+  cumhazEntry*(p[3] ) + (cumhazStop- cumhazStart)*(p[3] ))  ) # ***6th term***  
      -lik
    } 
    
    t2<-nlminb(inip,Lik_RRC,hessian=hess, control = list(rel.tol = 1e-6) )
    
    
    output <- c(alpha=t2$par[4], theta = t2$par[3],
                lambda = exp(t2$par[1]),
                rho =exp(t2$par[2])  )
    if (hess == T){
      
      std3 <- sqrt(diag(solve(hessian( Lik_RRC, t2$par))))
      
      Output2<-round( matrix( c( t2$par[4], std3[4], t2$par[3], std3[3] ,  exp(t2$par[1]), std3[1]*exp(t2$par[1]), exp( t2$par[2]) , std3[2]*exp(t2$par[2])  ), nrow=4, byrow=T),3)
      colnames(Output2) <- c("Estimate","SE")
      rownames(Output2) <- c("alpha","theta","Lambda","Rho" )
      return(Output2)
    }  else { output  }
    
  }
  
  
  
  
  else {
    if (randslope == T) {
      #We use purrr::map_dfr() as it allows perform a linear mixed model, perform prediction for every individual in the risk set at the specific event times then combines all the new subsets with predicted covariate into a single data.frame.
      predicted_Cov <- function(atf) {
        
        # Subsetting the dataset based on atf
        ds_subset <- Dat[Dat$ATF >= atf & Dat$A0 <= atf & Dat$start <= atf, ]
        
        
        fit <- lmer(zij ~ start+(1|id)+(0+start|id)+(1|cluster), data=ds_subset ,REML = FALSE, 
                    control = lmerControl(optimizer ="Nelder_Mead") )
        
        # Creating new data for prediction
        newSubset <- data.frame(id = unique(ds_subset$id), 
                                cluster = subset(ds_subset, !duplicated(id))$cluster,
                                A0 = subset(ds_subset, !duplicated(id))$A0,
                                ATF = subset(ds_subset, !duplicated(id))$ATF,
                                status = subset(ds_subset, !duplicated(id))$status, 
                                start = atf)
        
        
        # Making predictions 
        newSubset$predicted_Cov <- predict(fit, newdata = newSubset)
        
        return(newSubset)
      }
      
      LMM_Output<-tt[1:(nt-1)]|>purrr::map_dfr(predicted_Cov) 
      head(LMM_Output,10) 
      
      LMM_Output_Merged<- LMM_Output
      Dat2<- LMM_Output_Merged[LMM_Output_Merged$ATF>LMM_Output_Merged$start,]
      #Arrange data in the start stop format. The idea is to create start-stop intervals where the stop of an interval is the start of the next interval.
      #The last stop is equal to when an individual experiences the event/is censored. The first start is has BMD measurements predicted at A0 or measurements taken at entry for individuals without a past.
      Dat3<-Dat2 %>% group_by(id) %>%
        mutate(stop = lead(start, default = 1e6),
               stop = ifelse(stop > ATF, ATF, stop), .after = start)%>%
        filter(stop<= ATF) %>%
        ungroup()
      Dat3<- orderBy(~id, data=Dat3)
      Dat3$status[Dat3$stop<Dat3$ATF] <-0 # censor all intervals with exception of the last interval for an individual. The status of the last interval will be either censored or observed.
      dat11 <- Dat3[Dat3$stop>Dat3$start,] 
      #This next step is for obtaining only fully observed clusters (i.e only twin pairs without singletons). For the output below it will be based on fully observed clusters.
      # Fully observed clusters.
      dat1<-subset(dat11, ave(id, cluster, FUN = function(x) length(unique(x))) > 1)  
      
      #The following will be required in computation of the likelihood.
      Di<-aggregate(dat1$status, by=list ( dat1$cluster ) ,FUN=sum) #Total  number of NOT censored observations per cluster
      e<-colSums(Di[ 2 ] )  #e gives the sum of uncensored observation in the whole dataset.
      di<-Di [ 2 ] #1st column of Di is pair identifier so needs to be omitted. 
      GG<-length(unique( dat1$cluster) ) # Number of clusters
      Output3<-round(   matrix( c( e, GG  ) , nrow=2, byrow=T),5)
      
      #Likelihood estimation. This is the implementation of the proposed method.
      Lik_RRC<-function(p)
      {
        ###exp(p[1])<-lambda ; exp( p[2] )<-rho ;  p[3]<-theta ; p[4]<-alpha
        cumhazStart<-  as.numeric( aggregate(out ~ cluster, transform( dat1,
                                                                       out = exp(p[1])*(start^exp(p[2]))*exp( predicted_Cov*p[4] ) ), FUN = sum)[,2] )
        cumhazEntry<-as.numeric( aggregate(out ~ cluster, transform(subset(dat1, !duplicated(id)), out = exp(p[1])*(start^exp(p[2]) )*exp( predicted_Cov*p[4]) ), FUN = sum)[,2] )
        cumhazStop<- as.numeric( aggregate(out ~ cluster, transform( dat1,
                                                                     out =exp(p[1])*(stop^exp(p[2]))*exp(predicted_Cov*p[4] ) ), FUN = sum)[,2] )
        lnhaz<- as.numeric( aggregate(out ~ cluster, transform( dat1,
                                                                out = status*( predicted_Cov*p[4] + log(exp(p[1])*exp(p[2])*stop^(exp(p[2])-1) ) ) ), FUN = sum)[,2] )
        ## Likelihood Equation  ###
        lik<-e*log( (p[3] ) )- #***1st term***
          GG*log( gamma(1/(p[3] )) )+ #***2nd term***
          sum(log(gamma(di +1/(p[3] ))) )+#***3rd term***
          sum(lnhaz)+ ##***4th term***
          sum( (1/(p[3]) )*log(1+cumhazEntry*(p[3] )))- # *** 5 term***
          sum( ( di   +1/(p[3] ) )*log(1+  cumhazEntry*(p[3] ) + (cumhazStop- cumhazStart)*(p[3] ))  ) # ***6th term***  
        -lik
      }
      #Obtaining the MLE's
      
      t2<-nlminb(inip,Lik_RRC,hessian=hess, control = list(rel.tol = 1e-8) )
      
      
      output <- c(alpha=t2$par[4], theta = t2$par[3],
                  lambda = exp(t2$par[1]),
                  rho =exp(t2$par[2])  )
      if (hess == T){
        
        std3 <- sqrt(diag(solve(hessian( Lik_RRC, t2$par))))
        
        Output2<-round( matrix( c( t2$par[4], std3[4], t2$par[3], std3[3] ,  exp(t2$par[1]), std3[1]*exp(t2$par[1]), exp( t2$par[2]) , std3[2]*exp(t2$par[2])  ), nrow=4, byrow=T),3)
        colnames(Output2) <- c("Estimate","SE")
        rownames(Output2) <- c("alpha","theta","Lambda","Rho" )
        return(Output2)
      }  else { output  }
      }
    
  }
  
}

# Data Application step;
DatH<-read.csv( file="C:/Users/staammu/Desktop/SimulData.csv" )
attach(DatH)

fit_RRC( data=DatH, cluster = pair, id=id, inip=c(  log(0.01), log(2),0.5, 2 ) , Cov=zij,  CovTime=start, status=status,SurvTime=ATF, A0=A0, covariate=zij, hess=T, randslope=F)

