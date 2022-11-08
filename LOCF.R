library(numDeriv)
Dat<-read.csv( file="C:\\Users\\mmamm\\Desktop\\Data Manipulation\\result12611.csv" )
attach(Dat)
head(Dat,10) 
# Now arrange in start stop format
Dat3<-Dat2 %>% group_by(id) %>%
  mutate(stop = lead(start, default = Inf),
         stop = ifelse(stop > ATF, ATF, stop), .after = start)%>%
  filter(stop<= ATF) %>%
  ungroup()
## The last stop is equal to when an individual experiences the event/is censored. 
Dat4<-Dat3[Dat3$A0<Dat3$stop,]
head(Dat3,10)
Dat4$start<-pmax(Dat4$A0,Dat4$start) # The first start is when at entry
 head( Dat4, 10 )
length(unique(Dat4$id))

Dat4$status[Dat4$stop<Dat4$ATF] <-0 # censor all intervals with exception of the last interval for an individual. The status of the last interval will be either censored or observed.
dat11 <- Dat4[Dat4$stop>Dat4$start,] 

dat1<-subset(dat11, ave(id, pair, FUN = function(x) length(unique(x))) > 1)
sum(dat1$status)

#The following will be required in computation of the likelihood.
Di<-aggregate(dat1$status, by=list ( dat1$pair ) ,FUN=sum) #Total  number of NOT censored observations per cluster
e<-colSums(Di[ 2 ] )  #e gives the sum of uncensored observation in the whole dataset.
di<-Di [ 2 ] #1st column of Di is pair identifier so needs to be omitted. 
GG<-length(unique( dat1$pair) ) # Number of clusters
Output3<-round(   matrix( c( e, GG  ) , nrow=2, byrow=T),5)
rownames(Output3) <- c("Observed events","Clusters")
colnames(Output3)<-c("Number")
Output3


#============================================#
likelihood.Weibull2<-function(p)
{
  ###exp(p[1])<-lambda ; exp( p[2] )<-rho ;  p[3]<-theta ; p[4]<-alpha
  cumhazStart<-  as.numeric( aggregate(out ~ pair, transform( dat1,
                                                              out = exp(p[1])*(start^exp(p[2]))*exp( zij*p[4] ) ), FUN = sum)[,2] )
  cumhazEntry<-as.numeric( aggregate(out ~ pair, transform(subset(dat1, !duplicated(id)), out = exp(p[1])*(start^exp(p[2]) )*exp( zij*p[4]) ), FUN = sum)[,2] )
  cumhazStop<- as.numeric( aggregate(out ~ pair, transform( dat1,
                                                            out =exp(p[1])*(stop^exp(p[2]))*exp(zij*p[4] ) ), FUN = sum)[,2] )
  lnhaz<- as.numeric( aggregate(out ~ pair, transform( dat1,
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
#Obtaining the MLE's
initial2<-c(  log(0.003), log(1.5), 0.3,-3.5)
t2<-nlminb(initial2,likelihood.Weibull2,hessian=T, control = list(rel.tol = 1e-8) )
std3 <- sqrt(diag(solve(hessian(likelihood.Weibull2, t2$par))))
#alpha.SE<-data.frame(t2$par[4] , std3[4]) # alpha estimate+standard error
#theta.se<-data.frame(t2$par[3] , std3[3]) #Theta estimate+standard error
Output<-round( matrix( c( t2$par[4], std3[4], t2$par[3], std3[3]), nrow=2, byrow=T),3)
colnames(Output) <- c("Estimate","SE")
rownames(Output) <- c("alpha","theta")
#print(xtable(Output,digits=3),floating=FALSE,latex.environments=NULL) # Output to Latex table format.
Output


Output2<-round(   matrix( c( exp(t2$par[1]), exp( t2$par[2])  ) , nrow=2, byrow=T),5)
rownames(Output2) <- c("Lambda","Rho")
colnames(Output2)<-c("Estimate")
Output2
