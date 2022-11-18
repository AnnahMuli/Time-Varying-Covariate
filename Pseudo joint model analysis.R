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
    lambda       <-  exp(p[4] )
    rho      <- exp( p[6] )
    theta <- exp(p[2] ) # s=log(var.eps) <=> exp(s)=var.eps
    alpha  <- p[3]
    sigma_u <- exp( p[1]   )
    sigma_e <-  exp( p[5]  )
    n <- length(unique(LongData$pair))
    L.ni <-  as.vector(table(LongData$pair))
    yit<- LongData$zij-Beta0-LongData$bij-Beta1*LongData$start 
    L.yit<-split(yit, LongData$pair )
    S.status<- aggregate( SurvData$status  , by =  list(SurvData$pair)  ,FUN=sum )[2] 
    hit<- ( (lambda*rho*SurvData$ATF^(rho-1) )* exp(  alpha*(Beta0+SurvData$bij+Beta1*SurvData$ATF )  ) )^SurvData$status
    S.hit<-split(hit, SurvData$pair )
    S.ni<-as.vector(table(SurvData$pair))
    Hit<- unlist (sapply(1:nrow(SurvData),   function(i)  adaptIntegrate( function(s)
      exp(alpha*(Beta0 + SurvData$bij[i] ) )* s^(rho-1)* exp(alpha*Beta1*s) ,lowerLimit= c(0) , upperLimit= c( SurvData$ATF[i]) )[1]),
      use.names=FALSE)
    S.Hit<-split(Hit, SurvData$pair )
    H0it<- unlist (sapply(1:nrow(SurvData),   function(i)  adaptIntegrate( function(s)
      exp(alpha*(Beta0 + SurvData$bij[i] ) )* s^(rho-1)* exp(alpha*Beta1*s) ,lowerLimit= c(0) , upperLimit= c( SurvData$A0[i]) )[1]),
      use.names=FALSE)
    S0.Hit<-split(H0it, SurvData$pair )
    
    lik<-sum(log( unlist( 
      sapply(
        1:n, function(i)
          adaptIntegrate( function(uv)
            as.numeric(      ( (2*pi*(sigma_e^2 ) )^(-L.ni[i]/2) * exp(   -1/(2*(sigma_e^2 )  ) *  sum( (L.yit[[i]] - uv[1])^2  )  ) )*
                               ( prod( S.hit[[i]] )*( exp(alpha*uv[1]) )^S.status[i,]    )*
                               (( ( 1/theta+ sum(S0.Hit[[i]])*lambda*rho*exp(alpha*uv[1]) )^(  1/theta  ) )/ gamma(1/theta) )*
                               ( gamma( S.status[i,] + 1/theta ) /( ( 1/theta+ sum(S.Hit[[i]])*lambda*rho*exp(alpha*uv[1]) )^( S.status[i,] + 1/theta  ) ) ) )* 
              #   ( 1/(sqrt(2*pi* (sigma_u^2) ) ) )*exp( -uv[1]^2/(2*(sigma_u^2) )),
              dnorm(uv[1],mean=0,sd=sigma_u ) , 
            lowerLimit= c(-10 ) , upperLimit= c(10) )[1]
      )                                        
    ) ) )
    -lik
  }
  initial2<-c( log(0.1), log(0.5), 2,log(0.001) , log(0.6), log(2) )
  system.time(t2<- nlminb(initial2, likelihood.Joint, gradient = NULL, hessian = NULL,
                          control = list(trace=T,  rel.tol = 1e-6, step.min= 1e-3)) )
  cbind(sigma_u_JM = exp(t2$par[1]), Theta_JM=exp(t2$par[2]),
        alpha_JM= t2$par[3], Lambda_JM=exp(t2$par[4]),
        sigma_e_JM=exp(t2$par[5]),Rho_JM= exp(t2$par[6]) )
