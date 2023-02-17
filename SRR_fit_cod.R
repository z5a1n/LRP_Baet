#R version 4.1.2
library(ggplot2)
library(rstan)
library(rstanarm)
library(shinystan)

source("Functions.R")

DF <- read.csv("Data/cod_RS_data.csv")
R<-DF$Rec
S<-DF$SSB
maxS <- 200
maxR <- 250

P <- ggplot() + geom_point(data=DF, aes(y=Rec,x=SSB,colour=Year)) + 
  theme_classic() +
  labs(x='SSB', y="Recruitment") 

P + scale_x_continuous(limits=c(0,maxS), expand = c(0, 0)) +
   scale_y_continuous(limits=c(0,maxR), expand = c(0, 0)) 

##########################################################################################
# Bayesian fits using priors consistent with Perala et al. (2022)
# Final run use warmup of 10 000 and iter of 100 000 following Perala et al. (2022)
##########################################################################################
nw <- 10000
ni <- 100000

sBH_data <- list(Y=length(S),S=S,R=DF$R,RinfLow=0,RinfUp=max(R)*2,SqLow=0,SqUp=max(S)*2,cLow=0,cUp=5,sigmaLow=0,sigmaUp=1)
f_sBH <- stan(file="Stan/sBH_c5.stan",data=sBH_data,chains=4,warmup=nw,iter=ni,cores=1,refresh=0)

launch_shinystan(f_sBH)

SL_data <- list(Y=length(S),S=S,R=R,kLow=0,kUp=max(R)*2,SkLow=0,SkUp=max(S)*2,cLow=0,cUp=5,sigmaLow=0,sigmaUp=1)
f_SL <- stan(file="Stan/SL_c5.stan",data=SL_data,chains=4,warmup=nw,iter=ni,cores=1,refresh=0)

launch_shinystan(f_SL)

stan_dens(f_sBH, pars="c")
stan_dens(f_SL, pars="c")

# Extract parameter estimates
p_sBH <- sBH_pars(f_sBH)
p_SL <- SL_pars(f_SL)

P + geom_function(fun=function(x) (p_sBH$Rinf/(1+(p_sBH$Sq/x)^p_sBH$c)),colour="black",linetype=2) +
  geom_function(fun=function(x) (p_SL$k*(x/p_SL$Sk)^p_SL$c * exp(p_SL$c*(1-x/p_SL$Sk))),colour="blue",linetype=2) + 
  scale_x_continuous(limits=c(0,800), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,500), expand = c(0, 0)) 

