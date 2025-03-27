
#################################################################################
#  Last observation carried forward: Weibull baseline hazard and gamma frailty  #
#################################################################################
#                                                                               #
#  Its parameters are                                                           #
#   - data     : a data.frame containing all the variables;                     #
#   - cluster  : the name of the variable in data containing cluster IDs;       #
#   - id       : the name of the variable in data containing individual IDs;    #
#   - CovTime  : time when the covariate measurements were taken                #
#   - status   : censoring indicator                                            #
#   - SurvTime : Survival time                                                  #
#   - A0       : Age at entry                                                   #
#   - inip     : initial values for the parameters                              #
#################################################################################
#                                                                               #
#   Last modification : March  27, 2025                                         #
#                                                                               #
#################################################################################

fit_LOCF <- function(data, cluster, id, inip, CovTime, status, covariate,SurvTime,A0) 
{
  library(numDeriv); library(dplyr)
  
  Dat <-data
  Dat$id <- Dat[[deparse(substitute(id))]] 
  Dat$start <- Dat[[deparse(substitute(CovTime))]] 
  Dat$zij <-Dat[[deparse(substitute(covariate))]] 
  Dat$status <- Dat[[deparse(substitute(status))]] 
  Dat$cluster <- Dat[[deparse(substitute(cluster))]]
  Dat$ATF <- Dat[[deparse(substitute(SurvTime))]]
  Dat$A0 <- Dat[[deparse(substitute(A0))]]
  
  Dat3<-Dat %>% group_by(id) %>%
    mutate(stop = lead(start, default = 1e6),
           stop = ifelse(stop > ATF, ATF, stop), .after = start)%>%
    filter(stop<= ATF) %>%
    ungroup()
  ## The last stop is equal to when an individual experiences the event/is censored. 
  Dat4<-Dat3[Dat3$A0<Dat3$stop,]
  
  Dat4$start<-pmax(Dat4$A0,Dat4$start) # The first start is when at entry
  
  Dat4$status[Dat4$stop<Dat4$ATF] <-0 # censor all intervals with exception of the last interval for an individual. The status of the last interval will be either censored or observed.
  dat11 <- Dat4[Dat4$stop>Dat4$start,] 
  
  dat1<-subset(dat11, ave(id, cluster, FUN = function(x) length(unique(x))) > 1)
  
  
  #The following will be required in computation of the likelihood.
  Di<-aggregate(dat1$status, by=list ( dat1$cluster ) ,FUN=sum) #Total  number of NOT censored observations per cluster
  e<-colSums(Di[ 2 ] )  #e gives the sum of uncensored observation in the whole dataset.
  
  di<-Di [ 2 ] #1st column of Di is cluster identifier so needs to be omitted. 
  GG<-length(unique( dat1$cluster) ) # Number of clusters
  
  Lik_LOCF <- function(p){
    
    ###exp(p[1])<-lambda ; exp( p[2] )<-rho ;  p[3]<-theta ; p[4]<-alpha
    cumhazStart<-  as.numeric( aggregate(out ~ cluster, transform( dat1,
                                                                   out = exp(p[1])*(start^exp(p[2]))*exp( zij*p[4] ) ), FUN = sum)[,2] )
    cumhazEntry<-as.numeric( aggregate(out ~ cluster, transform(subset(dat1, !duplicated(id)), out = exp(p[1])*(start^exp(p[2]) )*exp( zij*p[4]) ), FUN = sum)[,2] )
    cumhazStop<- as.numeric( aggregate(out ~ cluster, transform( dat1,
                                                                 out =exp(p[1])*(stop^exp(p[2]))*exp(zij*p[4] ) ), FUN = sum)[,2] )
    lnhaz<- as.numeric( aggregate(out ~ cluster, transform( dat1,
                                                            out = status*( zij*p[4] + log(exp(p[1])*exp(p[2])*stop^(exp(p[2])-1) ) ) ), FUN = sum)[,2] )
    ## Likelihood Equation  ###
    lik<-e*log( (p[3] ) )- #***1st term***
      GG*log( gamma(1/(p[3] )) )+ #***2nd term***
      sum(log(gamma(di +1/(p[3] ))) )+#***3rd term***
      sum(lnhaz)+ ##***4th term***
      sum( (1/(p[3]) )*log(1+cumhazEntry*(p[3] )))- # *** 5 term***
      sum( ( di   +1/(p[3] ) )*log(1+  cumhazEntry*(p[3] ) + (cumhazStop- cumhazStart)*(p[3] ))  ) # ***6th term***
    -lik
  }
  
  
  t2<-nlminb(inip, Lik_LOCF,hessian=hess, control = list(rel.tol = 1e-8) )
  std3 <- sqrt(diag(solve(hessian( Lik_LOCF, t2$par))))
  
  Output<-round( matrix( c( t2$par[4], std3[4], t2$par[3], std3[3] ,  exp(t2$par[1]), std3[1]*exp(t2$par[1]), exp( t2$par[2]) , std3[2]*exp(t2$par[2])  ), nrow=4, byrow=T),3)
  colnames(Output) <- c("Estimate","SE")
  rownames(Output) <- c("alpha","theta","Lambda","Rho" )
  Output
  
  
}




#================================================================================================;
#Fitting to the data
DatH<-read.csv( file="C:\\Users\\mmamm\\Desktop\\Data Manipulation\\result12611.csv" )
attach(DatH)
fit_LOCF ( data=DatH, cluster = pair, id=id, inip=c(  log(0.1), log(1.5),0.3,-3.5), CovTime=start, status=status,SurvTime=ATF, A0=A0, covariate=zij)

