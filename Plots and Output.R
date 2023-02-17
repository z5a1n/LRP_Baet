#R version 4.1.2
library(ggplot2)
library(rootSolve)
source("Functions.R")
#output <- list(); saveRDS(output,"output.rda")
output <- readRDS("output.rda")

stock <- "BS"
#stock <- "GoR" 
#stock <- "NSAS" 
#stock <- "WBSS" 

D <- read.csv("Data/herring_RS_data.csv")
D$SSB <- D$SSB/1000 # put in kt
D$Rec <- D$Rec/1000000 # put in billions

DF <- D[D$Stock==stock,]

if(stock=="BS") {start = 500;pi <- 1; ageR <- 1; maxS <- 2000; maxR <- 20; x1<-1600; y1<-25}  # Used for plotting and root finding in ref pt calculations
if(stock=="GoR") {start = 50;pi <- 2; ageR <- 0; maxS <- 600; maxR <- 8; x1<-600; y1<-30} 
if(stock=="NSAS") {start = 1000;pi <- 3; ageR <- 0; maxS <- 10000; maxR <- 150; x1<-12000; y1<-1000} 
if(stock=="WBSS") {start = 200;pi <- 4; ageR <- 1; maxS <- 1000; maxR <- 10; x1<-1250; y1<-150} 

R<-DF$Rec; S<-DF$SSB

################################################################################
#Biological Data and Selectivity for Ref Pt Calcs 
################################################################################
BD <- read.csv("Data/herring_bio_data.csv")
BD <- BD[BD$Stock==stock,1:13] #exclude units and reference columns

if(stock %in% c("GoR")){cols <- 5:12}
if(stock %in% c("BS","WBSS","NSAS")){cols <- 5:13}
M <- colSums(BD[BD$Parameter=="M",cols])/nrow(BD[BD$Parameter=="M",cols]) # 1 to 8+ or 0 to 8+
waa <- colSums(BD[BD$Parameter=="waa",cols])/nrow(BD[BD$Parameter=="waa",cols])
mat <- colSums(BD[BD$Parameter=="mat",cols])/nrow(BD[BD$Parameter=="mat",cols])
sel <- colSums(BD[BD$Parameter=="vul",cols])/max(colSums(BD[BD$Parameter=="vul",cols]))

phi0 <- sum(survivorship_F(M=M,n_ages = length(M))*waa*mat) #g per recruit

##########################################################################################
# Parameter estimates from SRR_fit_herring.R
##########################################################################################

#extract parameters from final run
if(stock=="BS"){
  p_BH <- c(12.6159213,424.7663383)
  p_sBH <- c(11.4495992,335.532772,1.8018155)
  p_rick <- c(7.475654,638.7172959)
  p_SL <- c(8.6044496,602.5070237,2.1428783)
}
if(stock=="GoR"){
  p_BH <- c(5.6772294,81.7341825)
  p_sBH <- c(4.977087,66.5420464,2.2981945)
  p_rick <- c(3.1319455,105.7343989)
  p_SL <- c(3.6206014,103.2636184,3.2732624)
}
if(stock=="NSAS"){
  p_BH <- c(93.4441266,1192.918968)
  p_sBH <- c(60.2338751,448.036722,2.1544575)
  p_rick <- c(61.7492847,2409.936957)
  p_SL <- c(61.8293036,2607.964515,0.9405727)
}
if(stock=="WBSS"){
  p_BH <- c(3.9015715,59.0164263)
  p_sBH <- c(3.7837373,86.7006506,3.0031096)
  p_rick <- c(3.5217852,277.2197464)
  p_SL <- c(3.7632599,256.539358,1.7022433)
}
##########################################################################################
# HS fit
##########################################################################################

mle_hs <- optim(fn=HSnll,par=log(start),lower=log(min(S)), upper=log(max(S)), hessian = T, method="L-BFGS-B")
if(mle_hs$convergence!=0){print("Did not conv")}
Sad <- S
Sad[S>exp(mle_hs$par)] <- exp(mle_hs$par)
p_hs <- c(exp(mle_hs$par),exp(mean(log(R/Sad))))

####################################################################################
# Ref Pts      
####################################################################################
# Inflection points
if(p_sBH[3] > 1){ Binf_sBH <- p_sBH[2]/((p_sBH[3]+1)/(p_sBH[3]-1))^(1/p_sBH[3])} else {Bae_inf <- -10}
if(p_SL[3] > 1){ Binf_SL <- (1-1/sqrt(p_SL[3]))*p_SL[2]} else {Binf_SL <- -10} # -10 is a placeholder for NA (i.e., no inflection Point)

#Reference Points
RsBH <- RPcalc(M,waa,mat,sel,SRR="sBH",SRRpars=p_sBH,Req=mean(R),Smax=maxS)
RSL <- RPcalc(M,waa,mat,sel,SRR="SL",SRRpars=p_SL,Req=mean(R),Smax=maxS)
RBH <- RPcalc(M,waa,mat,sel,SRR="BH_Rinf",SRRpars=p_BH,Req=mean(R),Smax=maxS)
Rrick <- RPcalc(M,waa,mat,sel,SRR="SL",SRRpars=c(p_rick[1:2],1,p_rick[3]),Req=mean(R),Smax=maxS)

D_sBH <- data.frame(f=RsBH$f, eSSB = RsBH$eq_SSB, Y = RsBH$yield)
D_sBH <- D_sBH[D_sBH$eSSB>0.00001,]
D_sBH <- D_sBH[order(D_sBH$eSSB),]

D_BH <- data.frame(f=RBH$f, eSSB = RBH$eq_SSB, Y = RBH$yield)
D_BH <- D_BH[D_BH$eSSB>0.00001,]

D_SL <- data.frame(f=RSL$f, eSSB = RSL$eq_SSB, Y = RSL$yield)
D_SL <- D_SL[D_SL$eSSB>0.001 & D_SL$eSSB <= RSL$SSB0,]
D_SL <- D_SL[order(D_SL$eSSB),]

D_rick <- data.frame(f=Rrick$f, eSSB = Rrick$eq_SSB, Y = Rrick$yield)
D_rick <- D_rick[D_rick$eSSB>0.001 & D_rick$eSSB <= Rrick$SSB0,]
D_rick <- D_rick[order(D_rick$eSSB),]

# Allee effect threshold
Baet_sBH <- optimize(function(x){log(p_sBH[1]/(1+(p_sBH[2]/x)^p_sBH[3])/x)}, interval=c(0, maxS), maximum=TRUE)$maximum
Baet_SL <- optimize(function(x){log(p_SL[1] * (x/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-x/p_SL[2]))/x)}, interval=c(0, maxS), maximum=TRUE)$maximum
sBH_y_at_zero <- log(p_sBH[1]/(1+(p_sBH[2]/RsBH$SSB0)^p_sBH[3])/RsBH$SSB0)
SL_y_at_zero <- log(p_SL[1] * (RSL$SSB0/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-RSL$SSB0/p_SL[2]))/RSL$SSB0)

# Allee threshold
Bat_sBH <- uniroot.all(function(x) {log(p_sBH[1]/(1+(p_sBH[2]/x)^p_sBH[3])/x)-sBH_y_at_zero},c(0,RsBH$SSB0/2),tol=0.001)
Bat_SL <- uniroot.all(function(x) {log(p_SL[1] * (x/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-x/p_SL[2]))/x)-SL_y_at_zero},c(0,RSL$SSB0/2),tol=0.001)

if(length(Bat_SL)==0){Bat_SL<-0}
if(length(Bat_sBH)==0){Bat_sBH<-0}

##########################################################################################
# Plot SRR with unfished R/S
##########################################################################################

if(stock=="BS"){xmax=1600}
if(stock=="GoR"){xmax=200}
if(stock=="NSAS"){xmax=6000}
if(stock=="WBSS"){xmax=500}

fs=14

ggplot() +  theme_classic() + labs(x='SSB (kt)', y='Recruitment (billions)') + ggtitle(paste0(stock, " herring")) +
  scale_x_continuous(limits=c(0,xmax), expand = c(0, 0)) + scale_y_continuous(limits=c(0,maxR), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
        title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) +
  geom_function(fun=function(x) (p_sBH[1]/(1+(p_sBH[2]/x)^p_sBH[3])),colour="black",linetype=2,size=1.2) +
  geom_function(fun=function(x) (p_BH[1] / (1+p_BH[2]/x)),colour="black",size=1.2) +
  geom_function(fun=function(x) (p_SL[1] * (x/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-x/p_SL[2]))),colour="blue",linetype=2,size=1.2) +
  geom_function(fun=function(x) (p_rick[1] * (x/p_rick[2]) * exp((1-x/p_rick[2]))),colour="blue",size=1.2) +
  geom_segment(aes(x = 0, y = 0, xend = p_hs[1], yend =  p_hs[1]*p_hs[2]), colour = "green",size=1.2) +
  geom_segment(aes(x = p_hs[1], y = p_hs[1]*p_hs[2], xend = xmax, yend =  p_hs[1]*p_hs[2]), colour = "green",size=1.2) +
  geom_function(fun=function(x) (1/phi0*x),colour="grey",size=1)+ geom_point(data=DF, aes(y=Rec,x=SSB),colour="red4",size=2)

ggsave(filename=paste0("Figs/fign",pi,".png"),plot = last_plot(), width=7,height=5,units="in")

##########################################################################################
# Plot log(R/S)
##########################################################################################

ggplot() + theme_classic() + labs(x='SSB (kt)', y='ln(Recruitment/SSB)') + ggtitle(paste0(stock, " herring")) +
  scale_x_continuous(limits=c(0,xmax), expand = c(0, 0)) + #scale_y_continuous(limits=c(0,0.04), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
        title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) +
  geom_function(fun=function(x) (log(p_sBH[1]/(1+(p_sBH[2]/x)^p_sBH[3])/x)),colour="black",linetype=2,size=1.2) +
  geom_function(fun=function(x) (log(p_BH[1] / (1+p_BH[2]/x)/x)),colour="black",size=1.2) +
  geom_function(fun=function(x) (log(p_SL[1] * (x/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-x/p_SL[2]))/x)),colour="blue",linetype=2,size=1.2) +
  geom_function(fun=function(x) (log(p_rick[1] * (x/p_rick[2]) * exp((1-x/p_rick[2]))/x)),colour="blue",size=1.2) +
  geom_vline(xintercept = Binf_sBH, colour="red",linetype=1,size=1) +
  geom_vline(xintercept = Bat_sBH, colour="red",linetype=1,size=1) +
  geom_vline(xintercept = Baet_sBH, colour="red",linetype=1,size=1) +
  geom_vline(xintercept = Binf_SL, colour="orange",linetype=1,size=1) +
  geom_vline(xintercept = Bat_SL, colour="orange",linetype=1,size=1) +
  geom_vline(xintercept = Baet_SL, colour="orange",linetype=1,size=1) +
  geom_hline(yintercept = SL_y_at_zero, colour="grey",linetype=1,size=0.5) + 
  geom_point(data=DF, aes(y=log(Rec/SSB),x=SSB),colour="red4",size=2) 
  
ggsave(filename=paste0("Figs/fign",pi+4,".png"),plot = last_plot(), width=7,height=5,units="in")

##########################################################################################
# Plot Relative Yield vs Relative SSB
##########################################################################################

ggplot() +  theme_classic() + ggtitle(paste0(stock, " herring")) +
  geom_path(data=D_BH, aes(y=Y/max(Y),x=eSSB/max(eSSB)),size=1.2) +
  geom_path(data=D_sBH, aes(y=Y/max(Y),x=eSSB/max(eSSB)),linetype=2,size=1.2) +
  geom_path(data=D_rick, aes(y=Y/max(Y),x=eSSB/max(eSSB)),colour="blue",size=1.2) +
  geom_path(data=D_SL, aes(y=Y/max(Y),x=eSSB/max(eSSB)),colour="blue",linetype=2,size=1.2) +
  scale_x_continuous(limits=c(0,1), expand = c(0, 0)) + 
  scale_y_continuous(limits=c(0,1.01), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs),
        axis.text.y = element_text(size = fs),
        axis.title.x = element_text(size = fs),
        axis.title.y = element_text(size = fs),
        title = element_text(size = fs),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 20, 10, 10)) +
  labs(x='Relative SSB', y="Relative Yield")

ggsave(filename=paste0("Figs/fig",pi+8,".png"),plot = last_plot(), width=7,height=5,units="in")

##########################################################################################
# Plot Relative Yield vs F
##########################################################################################

if(stock=="BS"){fmax=0.6}
if(stock=="GoR"){fmax=1.5}
if(stock=="NSAS"){fmax=4.3}
if(stock=="WBSS"){fmax=1.75}

ggplot() +  theme_classic() + ggtitle(paste0(stock, " herring")) +
  geom_path(data=D_BH, aes(y=Y/max(Y),x=f),size=1.2) +
  geom_path(data=D_sBH, aes(y=Y/max(Y),x=f),linetype=2,size=1.2) +
  geom_path(data=D_rick, aes(y=Y/max(Y),x=f),colour="blue",size=1.2) +
  geom_path(data=D_SL, aes(y=Y/max(Y),x=f),linetype=2,colour="blue",size=1.2) +
  scale_x_continuous(limits=c(0,fmax), expand = c(0, 0.001)) + 
  scale_y_continuous(limits=c(0,1.01), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs),
        axis.text.y = element_text(size = fs),
        axis.title.x = element_text(size = fs),
        axis.title.y = element_text(size = fs),
        title = element_text(size = fs),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 20, 10, 10)) +
  labs(x='F', y="Relative Yield")

ggsave(filename=paste0("Figs/fig",pi+12,".png"),plot = last_plot(), width=7,height=5,units="in")

##########################################################################################
# Plot Relative Yield vs Relative SSB
##########################################################################################

if(stock=="BS"){ym=100}
if(stock=="GoR"){ym=50}
if(stock=="NSAS"){ym=1200}
if(stock=="WBSS"){ym=120}

#Relative Yield vs Relative SSB
ggplot() +  theme_classic() + ggtitle(paste0(stock, " herring")) +
  geom_path(data=D_BH, aes(y=Y,x=eSSB),size=1.2) +
  geom_path(data=D_rick, aes(y=Y,x=eSSB),colour="blue",size=1.2) +
  scale_x_continuous(limits=c(0,x1), expand = c(0, 0)) + 
  scale_y_continuous(limits=c(0,ym), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs),
        axis.text.y = element_text(size = fs),
        axis.title.x = element_text(size = fs),
        axis.title.y = element_text(size = fs),
        title = element_text(size = fs),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 20, 10, 10)) +
  labs(x='SSB (kt)', y="Yield (kt)") +
  geom_vline(xintercept = Baet_sBH, colour="red",linetype=1,size=1) +
  geom_vline(xintercept = Baet_SL, colour="orange",linetype=1,size=1) +
  geom_vline(xintercept = 0.4*RBH$SSBmsy, colour="black",linetype=4,size=1) +
  geom_vline(xintercept = 0.2*RBH$SSB0, colour="black",linetype=3,size=1) +
  geom_vline(xintercept = 0.4*Rrick$SSBmsy, colour="blue",linetype=4,size=1) +
  geom_vline(xintercept = 0.2*Rrick$SSB0, colour="blue",linetype=3,size=1) +
  geom_vline(xintercept = p_hs[1], colour="green",linetype=1,size=1) 

ggsave(filename=paste0("Figs/fig",pi+16,".png"),plot = last_plot(), width=7,height=5,units="in")

ggplot() +  theme_classic() + ggtitle(paste0(stock, " herring")) +
  geom_path(data=D_sBH, aes(y=Y,x=eSSB),linetype=2,size=1.2) +
  geom_path(data=D_SL, aes(y=Y,x=eSSB),colour="blue",linetype=2,size=1.2) +
  scale_x_continuous(limits=c(0,x1), expand = c(0, 0)) + 
  scale_y_continuous(limits=c(0,ym), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs),
        axis.text.y = element_text(size = fs),
        axis.title.x = element_text(size = fs),
        axis.title.y = element_text(size = fs),
        title = element_text(size = fs),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 20, 10, 10)) +
  labs(x='SSB (kt)', y="Yield (kt)") +
  geom_vline(xintercept = Baet_sBH, colour="red",linetype=1,size=1) +
  geom_vline(xintercept = Baet_SL, colour="orange",linetype=1,size=1) +
  geom_vline(xintercept = 0.4*RsBH$SSBmsy, colour="black",linetype=4,size=1) +
  geom_vline(xintercept = 0.2*RsBH$SSB0, colour="black",linetype=3,size=1) +
  geom_vline(xintercept = 0.4*RSL$SSBmsy, colour="blue",linetype=4,size=1) +
  geom_vline(xintercept = 0.2*RSL$SSB0, colour="blue",linetype=3,size=1) +
  geom_vline(xintercept = p_hs[1], colour="green",linetype=1,size=1) 

ggsave(filename=paste0("Figs/fig",pi+116,".png"),plot = last_plot(), width=7,height=5,units="in")

##########################################################################################
# Save Ref Pts
##########################################################################################

refpts <- c(RBH$SSB0,RsBH$SSB0,Rrick$SSB0,RSL$SSB0,
            RBH$SSBmsy,RsBH$SSBmsy,Rrick$SSBmsy,RSL$SSBmsy,
            RBH$Fmsy,RsBH$Fmsy,Rrick$Fmsy,RSL$Fmsy,
            max(D_BH$f[D_BH$Y>0]),max(D_sBH$f[D_sBH$Y>0]),max(D_rick$f[D_rick$Y>0]),max(D_SL$f[D_SL$Y>0]),
            Bat_sBH,Bat_SL,Binf_sBH,Binf_SL,Baet_sBH,Baet_SL,exp(mle_hs$par))

output[[pi]] <- refpts

saveRDS(output,"output.rda")

##########################################################################################
# Extra Plots [run for BS herring only]
##########################################################################################

WHAT <- rbind(cbind(D_sBH,G=rep("sBH",nrow(D_sBH))),
              #cbind(D_BH,G=rep("BH",nrow(D_BH))),
              #cbind(D_rick,G=rep("Ricker",nrow(D_rick))),
              cbind(D_SL,G=rep("SL",nrow(D_SL))))

#for BS fig
xmax=800
maxR=12
ggplot() +  
  theme_classic() + labs(x='SSB (kt)', y='Recruitment (billions)') + ggtitle(paste0(stock, " herring")) +
  scale_x_continuous(limits=c(0,xmax), expand = c(0, 0)) + scale_y_continuous(limits=c(0,maxR), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
        title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) +
  geom_function(fun=function(x) (p_sBH[1]/(1+(p_sBH[2]/x)^p_sBH[3])),colour="black",linetype=1,size=1.2) +
  geom_function(fun=function(x) (1/sum(survivorship_F(f=max(D_sBH$f),M=M,n_ages = length(M),sel=sel)*waa*mat)*x),colour="red",linetype=1,size=1) +
  geom_function(fun=function(x) (1/sum(survivorship_F(f=0.15,M=M,n_ages = length(M),sel=sel)*waa*mat)*x),colour="purple",linetype=1,size=1) +
  geom_function(fun=function(x) (1/sum(survivorship_F(f=0,M=M,n_ages = length(M),sel=sel)*waa*mat)*x),colour="grey",linetype=1,size=1) +
  geom_vline(xintercept = Bat_sBH, colour="grey",linetype=2,size=1) +
  geom_vline(xintercept = RsBH$SSB0, colour="grey",linetype=2,size=1) +
  geom_vline(xintercept = Baet_sBH, colour="red",linetype=2,size=1) +
  geom_vline(xintercept = D_sBH$eSSB[D_sBH$f==0.15], colour="purple",linetype=2,size=1) 
  
ggsave(filename=paste0("Figs/figextra1b.png"),plot = last_plot(), width=7,height=5,units="in")

ggplot() +  theme_classic() + ggtitle(paste0(stock, " herring")) +
  geom_point(data=WHAT, aes(y=eSSB,x=f,color=Y)) +
  scale_colour_continuous(type='viridis') +
  geom_hline(yintercept = Bat_sBH, colour="red",linetype=1,size=1) + 
  geom_hline(yintercept = Baet_sBH, colour="red",linetype=1,size=1) +
  geom_hline(yintercept = Binf_sBH, colour="red",linetype=1,size=1) +
  geom_hline(yintercept = Bat_SL, colour="orange",linetype=1,size=1) + 
  geom_hline(yintercept = Baet_SL, colour="orange",linetype=1,size=1) +
  geom_hline(yintercept = Binf_SL, colour="orange",linetype=1,size=1) +
  scale_x_continuous(limits=c(0,0.4), expand = c(0, 0.001)) + 
  scale_y_continuous(limits=c(0,x1), expand = c(0,0)) +
  labs(x='F', y="SSB (kt)",color="Yield (kt)") +
  geom_vline(xintercept = 0.15, colour="purple",linetype=2,size=1) 


ggsave(filename=paste0("Figs/figextra2.png"),plot = last_plot(), width=7,height=5,units="in")


##########################################################################################
##########################################################################################
##########################################################################################
# 3Ps cod
##########################################################################################
##########################################################################################
##########################################################################################
rm(list=ls())
library(ggplot2)
library(rootSolve)
source("Functions.R")
output <- readRDS("output.rda")

DF <- read.csv("Data/cod_RS_data.csv")
R<-DF$Rec
S<-DF$SSB
maxS <- 200#max(S)
maxR <- 250#max(R)

###############################################################################################################
BD <- read.csv("Data/cod_bio_data.csv")
cols<-3:ncol(BD)
M <- colSums(BD[BD$Parameter=="M",cols])/nrow(BD[BD$Parameter=="M",cols]) 
waa <- colSums(BD[BD$Parameter=="waa",cols])/nrow(BD[BD$Parameter=="waa",cols])
mat <- colSums(BD[BD$Parameter=="mat",cols])/nrow(BD[BD$Parameter=="mat",cols])
sel <- colSums(BD[BD$Parameter=="vul",cols])/max(colSums(BD[BD$Parameter=="vul",cols]))

phi0 <- sum(survivorship_F(M=M,n_ages = length(M))*waa*mat) #g per recruit
###############################################################################################################

start = 200;pi <- 50; ageR <- 1; maxS <- 200; maxR <- 250; x1<-800; y1<-150

##########################################################################################
# Parameter estimates from SRR_cod.R
##########################################################################################

p_sBH <- c(339.9403761,191.3384186,1.7228125)
p_SL <- c(223.9512085,335.2316159,1.8758351)

if(p_sBH[3] > 1){ Binf_sBH <- p_sBH[2]/((p_sBH[3]+1)/(p_sBH[3]-1))^(1/p_sBH[3])} else {Bae_inf <- -10}
if(p_SL[3] > 1){ Binf_SL <- (1-1/sqrt(p_SL[3]))*p_SL[2]} else {Binf_SL <- -10} # -10 is a placeholder for NA (i.e., no inflection Point)

#Reference Points
RsBH <- RPcalc(M,waa,mat,sel,SRR="sBH",SRRpars=p_sBH,Req=mean(R),Smax=maxS*4)
RSL <- RPcalc(M,waa,mat,sel,SRR="SL",SRRpars=p_SL,Req=mean(R),Smax=maxS*4)

D_sBH <- data.frame(f=RsBH$f, eSSB = RsBH$eq_SSB, Y = RsBH$yield)
D_sBH <- D_sBH[D_sBH$eSSB>0.00001,]
D_sBH <- D_sBH[order(D_sBH$eSSB),]

D_SL <- data.frame(f=RSL$f, eSSB = RSL$eq_SSB, Y = RSL$yield)
D_SL <- D_SL[D_SL$eSSB>0.001 & D_SL$eSSB <= RSL$SSB0,]
D_SL <- D_SL[order(D_SL$eSSB),]

# Allee effect threshold
Baet_sBH <- optimize(function(x){log(p_sBH[1]/(1+(p_sBH[2]/x)^p_sBH[3])/x)}, interval=c(0, maxS), maximum=TRUE)$maximum
Baet_SL <- optimize(function(x){log(p_SL[1] * (x/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-x/p_SL[2]))/x)}, interval=c(0, maxS), maximum=TRUE)$maximum
sBH_y_at_zero <- log(p_sBH[1]/(1+(p_sBH[2]/RsBH$SSB0)^p_sBH[3])/RsBH$SSB0)
SL_y_at_zero <- log(p_SL[1] * (RSL$SSB0/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-RSL$SSB0/p_SL[2]))/RSL$SSB0)

# Allee threshold
Bat_sBH <- uniroot.all(function(x) {log(p_sBH[1]/(1+(p_sBH[2]/x)^p_sBH[3])/x)-sBH_y_at_zero},c(0,100),tol=0.001)
Bat_SL <- uniroot.all(function(x) {log(p_SL[1] * (x/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-x/p_SL[2]))/x)-SL_y_at_zero},c(0,100),tol=0.001)

xmax=200;fmax=0.4;ym=50
fs=14

ggplot() + geom_point(data=DF, aes(y=log(Rec/SSB),x=SSB)) + theme_classic() + labs(x='SSB (kt)', y='ln(Recruitment/SSB)') + ggtitle("3Ps cod") +
  scale_x_continuous(limits=c(0,600), expand = c(0, 0)) + #scale_y_continuous(limits=c(0,0.04), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
        title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) +
  geom_function(fun=function(x) (log(p_sBH[1]/(1+(p_sBH[2]/x)^p_sBH[3])/x)),colour="black",linetype=2,size=1.2) +
  geom_function(fun=function(x) (log(p_SL[1] * (x/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-x/p_SL[2]))/x)),colour="blue",linetype=2,size=1.2) +
  geom_vline(xintercept = Binf_sBH, colour="red",linetype=1,size=1) +
  geom_vline(xintercept = Bat_sBH, colour="red",linetype=2,size=1) +
  geom_vline(xintercept = Baet_sBH, colour="red",linetype=2,size=1) +
  geom_vline(xintercept = Binf_SL, colour="orange",linetype=1,size=1) +
  geom_vline(xintercept = Bat_SL, colour="orange",linetype=2,size=1) +
  geom_vline(xintercept = Baet_SL, colour="orange",linetype=2,size=1) +
  geom_hline(yintercept = SL_y_at_zero, colour="grey",linetype=1,size=0.5) 

ggsave(filename=paste0("Figs/fig23.png"),plot = last_plot(), width=7,height=5,units="in")

ggplot() + theme_classic() + labs(x='SSB (kt)', y='Recruitment (millions)') + ggtitle("3Ps cod") +
  scale_x_continuous(limits=c(0,xmax), expand = c(0, 0)) + scale_y_continuous(limits=c(0,maxR), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
        title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) +
  geom_function(fun=function(x) (p_sBH[1]/(1+(p_sBH[2]/x)^p_sBH[3])),colour="black",linetype=2,size=1.2) +
  geom_function(fun=function(x) (p_SL[1] * (x/p_SL[2])^p_SL[3] * exp(p_SL[3]*(1-x/p_SL[2]))),colour="blue",linetype=2,size=1.2) +
  geom_function(fun=function(x) (1/phi0*x),colour="grey",size=1) + 
  geom_point(data=DF, aes(y=Rec,x=SSB),colour="red4",size=2) 

ggsave(filename=paste0("Figs/fign30.png"),plot = last_plot(), width=7,height=5,units="in")


ggplot() +  theme_classic() + ggtitle("3Ps cod") +
  geom_path(data=D_sBH, aes(y=Y/max(Y),x=eSSB/max(eSSB)),linetype=2,size=1.2) +
  geom_path(data=D_SL, aes(y=Y/max(Y),x=eSSB/max(eSSB)),colour="blue",linetype=2,size=1.2) +
  scale_x_continuous(limits=c(0,1), expand = c(0, 0)) + 
  scale_y_continuous(limits=c(0,1.01), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs),
        axis.text.y = element_text(size = fs),
        axis.title.x = element_text(size = fs),
        axis.title.y = element_text(size = fs),
        title = element_text(size = fs),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 20, 10, 10)) +
  labs(x='Relative SSB', y="Relative Yield")

ggsave(filename=paste0("Figs/fig31.png"),plot = last_plot(), width=7,height=5,units="in")


ggplot() +  theme_classic() + ggtitle("3Ps cod") +
  geom_path(data=D_sBH, aes(y=Y/max(Y),x=f),linetype=2,size=1.2) +
  geom_path(data=D_SL, aes(y=Y/max(Y),x=f),linetype=2,colour="blue",size=1.2) +
  scale_x_continuous(limits=c(0,fmax), expand = c(0, 0.001)) + 
  scale_y_continuous(limits=c(0,1.01), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs),
        axis.text.y = element_text(size = fs),
        axis.title.x = element_text(size = fs),
        axis.title.y = element_text(size = fs),
        title = element_text(size = fs),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 20, 10, 10)) +
  labs(x='F', y="Relative Yield")

ggsave(filename=paste0("Figs/fig32.png"),plot = last_plot(), width=7,height=5,units="in")


#Relative Yield vs Relative SSB
ggplot() +  theme_classic() + ggtitle("3Ps cod") +
  geom_path(data=D_sBH, aes(y=Y,x=eSSB),linetype=2,size=1.2) +
  geom_path(data=D_SL, aes(y=Y,x=eSSB),colour="blue",linetype=2,size=1.2) +
  scale_x_continuous(limits=c(0,x1), expand = c(0, 0)) + 
  scale_y_continuous(limits=c(0,ym), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs),
        axis.text.y = element_text(size = fs),
        axis.title.x = element_text(size = fs),
        axis.title.y = element_text(size = fs),
        title = element_text(size = fs),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 20, 10, 10)) +
  labs(x='SSB (kt)', y="Yield (kt)") +
  geom_vline(xintercept = Baet_sBH, colour="red",linetype=1,size=1) +
  geom_vline(xintercept = Baet_SL, colour="orange",linetype=1,size=1) +
  geom_vline(xintercept = 0.4*RsBH$SSBmsy, colour="black",linetype=3,size=1) +
  geom_vline(xintercept = 0.2*RsBH$SSB0, colour="black",linetype=4,size=1) +
  geom_vline(xintercept = 0.4*RSL$SSBmsy, colour="blue",linetype=3,size=1) +
  geom_vline(xintercept = 0.2*RSL$SSB0, colour="blue",linetype=4,size=1) +
  geom_vline(xintercept = max(S), colour="green",linetype=1,size=1) 

ggsave(filename=paste0("Figs/fig33.png"),plot = last_plot(), width=7,height=5,units="in")


refpts <- c(NA,RsBH$SSB0,NA,RSL$SSB0,
            NA,RsBH$SSBmsy,NA,RSL$SSBmsy,
            NA,RsBH$Fmsy,NA,RSL$Fmsy,
            NA,max(D_sBH$f[D_sBH$Y>0]),NA,max(D_SL$f[D_SL$Y>0]),
            Bat_sBH,Bat_SL,Binf_sBH,Binf_SL,Baet_sBH,Baet_SL,NA)

output[[5]] <- refpts

saveRDS(output,"output.rda")