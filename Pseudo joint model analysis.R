library(survival); library(parfm); library(frailtypack); library(xtable); library(dplyr); library(numDeriv);library(lme4)
library(doBy); library(cubature)

Dat<-read.csv( file="C:\\Users\\mmamm\\Desktop\\Data Manipulation\\result12611.csv" )

Dat<-Dat[Dat$ATF>Dat$start,]
head(Dat,10)
dim(Dat)

model = lmer(zij ~ start+(1|id)+(1|pair), data=Dat)
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


ui<-ranef(model)$pair
Ui_pair<-cbind(pair=rownames(ui),ui)
names(Ui_pair)[2] <- "ui"
Merged_data1<-merge(Merged_data,Ui_pair, by="pair")
head(Merged_data1)
dim(Merged_data1)
#Merged_data1<-orderBy(~id, data=Merged_data1)
LongData11<-Merged_data1[Merged_data1$ATF>Merged_data1$start,]
LongData1<-LongData11[1:100,]
dim(LongData1)
LongData<-subset(LongData1, ave(id, pair, FUN = function(x) length(unique(x))) > 1)

dim(LongData)

SurvData1=subset(LongData, !duplicated(id) )
SurvData<-subset(SurvData1, ave(id, pair, FUN = function(x) length(unique(x))) > 1)
dim(SurvData)
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
## Weibull baseline hazard
## h0 = lambda*rho*t^(rh0-1)   and   H0 = lambda*t^(rho)
## ---------------------------------------------------------------------
likelihood.Joint<-function(p)
  { 
  ids <- SurvData$id
  n_i <- length(unique(ids))
  G<-length(unique(LongData$pair)) # #G: Total number of clusters
  nij <- as.vector(tapply(ids, ids, length)) # same as as.vector(table(long$id))  # #nij: Number of measurements per individual
  nij1j2 <-as.vector(table(LongData$pair))
  H_it<- rep(NA, n_i)
  likk<- rep(NA, G)
lambda       <-  exp(p[1] )
rho      <-exp( p[2] )
theta <- exp(p[3]) # 
alpha  <- p[4]
sigma_u <-  exp(p[5]  )
sigma_e<-  exp( p[6]  )
Hit<-rep(NA, n_i)
for(i in 1:G) ## s<-3
{
  likk[i]<-adaptIntegrate( function(uv)  as.numeric( sapply( split( ( ( lambda*rho*SurvData$ATF^(rho-1)* uv[2]*exp(  alpha*(SurvData$bij+uv[1]+Beta1*SurvData$ATF )  ) )^SurvData$status )*
 (for(id in 1:n_i){ Hit[id]<-
 integrate(function(s)
 lambda*rho*s^(rho-1)*uv[2]* exp( alpha*(Beta0+SurvData$bij[id]+uv[1]+Beta1*s) ),
  0, SurvData$ATF[id] )$value
}) , SurvData$pair),prod) *((2*pi*(sigma_e^2) )^(-nij1j2/2) * exp(   -1/(2*(sigma_e^2)) * 
 aggregate((LongData$zij-Beta0-LongData$bij-uv[1]-Beta1*LongData$start)^2 , by=list ( LongData$pair ) ,FUN=sum)[,2]  )) ) [i]*
  dgamma(x=uv[2],shape=1/theta, scale =theta )* dnorm(uv[1],mean=0,sd=sigma_u )  , 
 lowerLimit= c(-20,0 ) , upperLimit= c(20, 100) )[1]
}
lik<-sum(log( unlist( likk ) ) )
-lik
}
initial2<-c(log(0.01), log(2), log(0.5), -3.5,log(0.09),log(0.04))
t2<-nlminb(initial2,likelihood.Joint, control=list(trace=TRUE))
t2