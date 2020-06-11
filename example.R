library(mvtnorm)
library(parallel)
library(evd)

## Set the working directory as you want ##
setwd("~/max_id_project/")

## Load the functions we need ##
source("Tools_Functions.R")

########################################
###  model simulation and inference #### 
########################################

### We only consider one spatial covariates in this toy example ###
### Therefore, we only have lambda.0 and lambda.1 linked to the range parameter in the dependence structure. ###

# parameter space psi={alpha, beta,lambda.0,lambda.1,nu}

# model 1:alpha>0,beta=0,lambda.0,lambda.1=0,nu=0; stationary max-stable 

# model 2:alpha>0,beta=0,lambda.0,lambda.1,nu=0; non-stationary max-stable 

# model 3:alpha>0,beta>0,lambda.0,lambda.1=0,nu=0; stationary simple max-id

# model 4:alpha>0,beta>0,lambda.0,lambda.1,nu=0; non-stationary simple max-id 

# model 5:alpha>0,beta>0,lambda.0,lambda.1=0,nu; stationary general max-id

# model 6:alpha>0,beta>0,lambda.0,lambda.1,nu; non-stationary general max-id

### take the most complicated model, model 6, for example ###
### simulate the dataset ###

alpha <- 1 
beta <- 1
lambda.0 <- -0.5 #intercept for the range parameter
lambda.1 <- 1 ## slope for the spatial covariate. 
nu <- 0.5 

### parameters not used in this example ###
lambda.2 <- 0 # fixed to zero in this model
a <- NULL; #anisotropy parameter if any
theta <- NULL ## rotation parameter if any

D = 5 ## D*D grid 
x <- y <- c(1:D)/(D+1) #x-axis and y-axis on the grid
cutoff <- 3/(D+1) # cutoff distance for six order neighbors
coord = as.matrix(expand.grid(x,y))

### covariates for the intercept lambda.0 and lambda.1 (slope for the temporal trend) ###
reg=cbind(1,2*pnorm(coord[,1],0.5,0.25)-1) ## spatial covariates. 
reg.t = NULL ## temporal covariates if any

n = 20 ### temporal length 
coord = as.matrix(dist(coord)) ## coord here represents distance matrix

## Parameter space separated into two: parR and parGauss
parR = c(alpha,beta)
parGauss = list(lambda=c(lambda.0,lambda.1),lambda.t=lambda.2,a=a,theta=theta,nu=nu,type=4) ##type here is just to denote the format of parameter space, fixed to be 4 in this toy example, please don't change it.###

### sample the max-id processes on the grid (takes a minute)###
ncores=40 ## number of cores available.
Z<-rmaxidspat(n,coord,parR,parGauss,reg=reg,reg.t=reg.t,N=1000,ncores=ncores) ## you can use multiple cores to speed up ##

### transform the original data into pseudo uniform scores. 
U <- matrix(pG(Z,parR),nrow=n)

### fixed the random seed ###
set.seed(98764672)

### Set the initial values to be around the true values for simplification; Noticed that the parR should be in log scale here.
fixed <- c(F,F,F,F,T,F) ## corresponding to c(alpha, beta, lambda.0,lambda.1,lambda.2,nu)
init <- c(parR,parGauss$lambda,parGauss$lambda.t,parGauss$nu)
init[!fixed] = init[!fixed] + runif(sum(!fixed),0.1,0.3)
init[1:2] <- log(init[1:2])

### Fit the model. It will take hours. ###
fit.result <- fit.pw.parallel(init=init,datU = U,coord=coord,reg=reg,reg.t=reg.t,cutoff=cutoff,proppairs= 1,fixed=fixed,optim =T, hessian=F,sandwich=F,eps = 10^(-2), print.par.file=NULL,ncores=ncores,fit.load=F, fit.save=F,fit.file=NULL)

##########################################################
######## Fit using the real data in our application ######
##########################################################

## if we use the real data to fit the marginal model, please load data.Rdata ##
load("data.Rdata")
y = as.vector(data)
id = rep(1:ncol(data),each=nrow(data))
t = rep((1:nrow(data))/nrow(data),times=44)
alt = rep(reg[,2],each=nrow(data))
lat = rep(data.info$LAT,each=nrow(data))
lon = rep(data.info$LON,each=nrow(data))
dat.mar.fit <- data.frame(y=y,t=t,alt=alt,lat=lat,lon=lon,id=id)[!is.na(y),]
# fit the marginal model
marginal.fit<-gam(list(y~ti(lon,lat,alt,k=5)+ti(t),~1,~1),family=gevlss,data=dat.mar.fit)

fitted.parameters <- fitted(marginal.fit)

## Figure 5 in the manuscript ##
pdf("return_level.pdf",width = 3*3,height=3)
par(mfrow=c(1,3),cex.main=2,mgp=c(2,1,0),cex.lab=2,mar=c(2,4,2,0.5))
s.list = c(9,27,32)
for(s in s.list){
  ylab=expression("Temperature ("*~degree*C*")")
  level.y <- (1917+c(1:101))[!is.na(data[,s])]
  level.1000<-(fitted.parameters[dat.mar.fit$id==s,1]+qgev(1-1/1000,scale=exp(2.871183),shape=-0.1986207))/10
  level.100<-(fitted.parameters[dat.mar.fit$id==s,1]+qgev(1-1/100,scale=exp(2.871183),shape=-0.1986207))/10
  level.10<-(fitted.parameters[dat.mar.fit$id==s,1]+qgev(1-1/10,scale=exp(2.871183),shape=-0.1986207))/10
  level.trend<-fitted.parameters[dat.mar.fit$id==s,1]/10
  level.data<-c(data[,s],data.2019[s])/10
  plot(1917+c(1:102),level.data,col=c(rep("black",101),"red"),ylim=c(min(level.data,na.rm = T)-0.5,max(level.1000)+0.5),cex=0.8,main=paste0("Station ",s),xlab="",ylab=ylab)
  lines(level.y,level.trend,col="black",lwd=1,lty=1)
  lines(level.y,level.1000,col=rgb(1,0,0,1),lwd=1)
  lines(level.y,level.100,col=rgb(1,0,0,0.8),lwd=1)
  lines(level.y,level.10,col=rgb(1,0,0,0.5),lwd=1)
  abline(h=data.2019[s]/10,lty=2,col="black")
  abline(v=c(1950,1975,2000,2018),lty=2)
}
dev.off()

### fit the dependence model using the real data ###
fit.result <- fit.pw.parallel(init=init,datU = U,coord=coord,reg=reg,reg.t=reg.t,cutoff=Inf,proppairs= 1,fixed=rep(F,6),optim =T, hessian=F,sandwich=F,eps = 10^(-2), print.par.file=NULL,ncores=ncores,fit.load=F, fit.save=F,fit.file=NULL)

  