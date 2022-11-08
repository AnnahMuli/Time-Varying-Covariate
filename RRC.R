library(survival); library(xtable); library(dplyr); library(numDeriv);library(lme4)
 library(doBy)
Dat<-read.csv( file="C:\\Users\\mmamm\\Desktop\\Data Manipulation\\result12611.csv" )
attach(Dat)
#Individuals with first BMD measurement at entry i.e they do not have past BMD values at entry.
Dat0<-Dat[order(start),]
Dat22<-subset(Dat0, !duplicated(id) )
Dat33<-Dat22[ Dat22$start==Dat22$A0, ]
Datt<-Dat33 %>% rename(predicted_BMD=zij )
Subset1<-Datt[c(1,2,5,6,7,3,4 )] # To help during merging and arranging data in start stop format.
head(Subset1)

#Ordering the unique event times and age at entry:
tt <- sort(unique( c(Dat$A0, Dat$ATF[status==1] ) )) # Unique event times or age at entry.
nt <- length(tt)

#Because at age 50, no individual has a past we don't do predictions at 50.
#At age 51 and 52 individuals at risk do not have repeated measurements, we fit a linear mixed model with intercept at twin grouping level. For the other event times, we fit a linear mixed model with random intercept at both individual and cluster level. Therefore, I analyze these two risk sets separately from the other risk sets.
#We use covariate measurements taken before time tt for the individuals in the risk set and predict covariate values for the individuals at the respective ages.
#We use purrr::map_dfr() as it allows perform a linear mixed model for each risk set, perform prediction for every individual in the risk set at the specific event times then combines all the new subsets with predicted covariate into a single data.frame.
predicted_BMD1 <- function (atf) {
  ds_subset<-  Dat[Dat$ATF >= atf & Dat$A0<=atf & Dat$start<atf, ]
  fit = lmer(zij ~ start + (1|pair), data = ds_subset)
  newSubset <- data.frame(id = unique(ds_subset$id), pair =subset( ds_subset,!duplicated(id) )$pair,
        #PublicID =subset( ds_subset,!duplicated(id) )$PublicID ,
                          A0 =subset(ds_subset, !duplicated(id) )$A0,
                          ATF =subset(ds_subset, !duplicated(id) )$ATF,
                          status =subset(ds_subset, !duplicated(id) )$status, start = atf )
  newSubset$predicted_BMD <- predict(fit, newdata = newSubset)
  newSubset
}

LMM_Output1<-tt[2:3] |> purrr::map_dfr(predicted_BMD1)
head(LMM_Output1,10)
#Next is to fit multiple linear mixed models for each of the remaining unique event/age at entry times and perform prediction too. 

predicted_BMD <- function (atf) {
  ds_subset<-  Dat[Dat$ATF >= atf & Dat$A0<=atf & Dat$start<atf, ]
fit = lmer(zij ~ start  + (1 | id) + (1|pair), data = ds_subset) # Model fit
newSubset <- data.frame(id = unique(ds_subset$id), pair =subset( ds_subset,!duplicated(id) )$pair,
        # PublicID =subset( ds_subset,!duplicated(id) )$PublicID ,
                              A0 =subset(ds_subset, !duplicated(id) )$A0,
                              ATF =subset(ds_subset, !duplicated(id) )$ATF,
                              status =subset(ds_subset, !duplicated(id) )$status, start = atf )
  newSubset$predicted_BMD <- predict(fit, newdata = newSubset) # Prediction of subject specific covariate at the event time
  newSubset
}

LMM_Output<-tt[4:nt]|>purrr::map_dfr(predicted_BMD)
head(LMM_Output,10) 

LMM_Output_Merged<-rbind( Subset1, LMM_Output1, LMM_Output)
Dat2<- LMM_Output_Merged[LMM_Output_Merged$ATF>LMM_Output_Merged$start,]

#Next is to arrange data in the start stop format. The idea is to create start-stop intervals where the stop of an interval is the start of the next interval.
#The last stop is equal to when an individual experiences the event/is censored. The first start is has BMD measurements predicted at A0 or measurements taken at entry for individuals without a past.
Dat3<-Dat2 %>% group_by(id) %>%
  mutate(stop = lead(start, default = Inf),
         stop = ifelse(stop > ATF, ATF, stop), .after = start)%>%
  filter(stop<= ATF) %>%
  ungroup()
Dat3<- orderBy(~id, data=Dat3)
Dat3$status[Dat3$stop<Dat3$ATF] <-0 # censor all intervals with exception of the last interval for an individual. The status of the last interval will be either censored or observed.
dat11 <- Dat3[Dat3$stop>Dat3$start,] 

#This next step is for obtaining only fully observed clusters (i.e only twin pairs without singletons). For the output below it will be based on fully observed clusters.
dat1<-subset(dat11, ave(id, pair, FUN = function(x) length(unique(x))) > 1)  # Fully observed clusters.
#The following will be required in computation of the likelihood.
Di<-aggregate(dat1$status, by=list ( dat1$pair ) ,FUN=sum) #Total  number of NOT censored observations per cluster
e<-colSums(Di[ 2 ] )  #e gives the sum of uncensored observation in the whole dataset.
di<-Di [ 2 ] #1st column of Di is pair identifier so needs to be omitted. 
GG<-length(unique( dat1$pair) ) # Number of clusters
Output3<-round(   matrix( c( e, GG  ) , nrow=2, byrow=T),5)
rownames(Output3) <- c("Observed events","Clusters")
colnames(Output3)<-c("Number")
Output3

#Likelihood estimation based. This is the implementation of the proposed method.
#============================================#
likelihood.Weibull2<-function(p)
{
  ###exp(p[1])<-lambda ; exp( p[2] )<-rho ;  p[3]<-theta ; p[4]<-alpha
cumhazStart<-  as.numeric( aggregate(out ~ pair, transform( dat1,
     out = exp(p[1])*(start^exp(p[2]))*exp( predicted_BMD*p[4] ) ), FUN = sum)[,2] )
cumhazEntry<-as.numeric( aggregate(out ~ pair, transform(subset(dat1, !duplicated(id)), out = exp(p[1])*(start^exp(p[2]) )*exp( predicted_BMD*p[4]) ), FUN = sum)[,2] )
cumhazStop<- as.numeric( aggregate(out ~ pair, transform( dat1,
       out =exp(p[1])*(stop^exp(p[2]))*exp(predicted_BMD*p[4] ) ), FUN = sum)[,2] )
lnhaz<- as.numeric( aggregate(out ~ pair, transform( dat1,
        out = status*( predicted_BMD*p[4] + log(exp(p[1])*exp(p[2])*stop^(exp(p[2])-1) ) ) ), FUN = sum)[,2] )
  ## Likelihood Equation ###
   lik<-e*log( (p[3] ) )- #***1st term***
    GG*log( gamma(1/(p[3] )) )+ #***2nd term***
   sum(log(gamma(di +1/(p[3] ))) )+#***3rd term***
  sum(lnhaz)+ ##***4th term***
 sum( (1/(p[3]) )*log(1+cumhazEntry*(p[3] )))- # *** 5 term***
sum( ( di   +1/(p[3] ) )*log(1+  cumhazEntry*(p[3] ) + (cumhazStop- cumhazStart)*(p[3] ))  ) # ***6th term***  
  -lik
}
#Obtaining the MLE's
initial2<-c(  log(0.001), log(1.5),0.3,-3.5)
t2<-nlminb(initial2,likelihood.Weibull2,hessian=T, control = list(rel.tol = 1e-8) )
std3 <- sqrt(diag(solve(hessian(likelihood.Weibull2, t2$par))))
#alpha.SE<-data.frame(t2$par[4] , std3[4]) # alpha estimate+standard error
#theta.se<-data.frame(t2$par[3] , std3[3]) #Theta estimate+standard error
Output<-round( matrix( c( t2$par[4], std3[4], t2$par[3], std3[3]), nrow=2, byrow=T),3)
colnames(Output) <- c("Estimate","SE")
rownames(Output) <- c("alpha","theta")
#print(xtable(Output,digits=3),floating=FALSE,latex.environments=NULL) # Output to Latex table format.
Output

#We can also obtain the Weibull baseline hazard parameters as follows:
Output2<-round(   matrix( c( exp(t2$par[1]), exp( t2$par[2])  ) , nrow=2, byrow=T),5)
rownames(Output2) <- c("Lambda","Rho")
colnames(Output2)<-c("Estimate")
Output2
