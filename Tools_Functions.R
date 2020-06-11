#######################################################################################################
### intensity of the radial variable R, random number generator
#######################################################################################################
#Density function (PDF) of the point process of R
dF <- function(r,parR,log=FALSE){
  alpha <- parR[1]
  beta <- parR[2]
  logval <- c()
  logval[is.na(r)] <- NA
  ind <- r>0 & !is.na(r)
  if(beta==0){
    logval[ind] <- -(alpha+1)*log(r[ind])+log(alpha)
  } else{# intensity function
    logval[ind] <- -alpha/beta*(r[ind]^beta-1)+log(beta*r[ind]^(-beta-1)+alpha*r[ind]^(-1))
  }
  logval[r<=0 & !is.na(r)] <- -Inf
  if(log){
    return( logval )	
  } else{
    return( exp(logval) )
  }
}

#the intensity dF integrated over [x,infty)
upperF <- function(x,parR,log=FALSE){
  alpha <- parR[1]
  beta <- parR[2]
  logval <- c()
  logval[is.na(x)] <- NA
  ind <- x>0 & !is.na(x)
  if(beta==0){
    logval[ind] <- -alpha*log(x[ind])
  } else{
    logval[ind] <- -beta*log(x[ind]) - alpha/beta*(x[ind]^beta-1)
  }
  logval[x<=0 & !is.na(x)] <- Inf
  if(log){
    return( logval )	
  } else{
    return( exp(logval) )
  }
}

#inverse function of the intensity dF integrated over [x,infty)
upperFinv <- function(y,parR,log=FALSE){
  alpha <- parR[1]
  beta <- parR[2]
  logval <- c()
  for(i in 1:length(y)){
    if(!is.na(y[i]) & y[i]>0){
      fun <- function(x){
        return(log(y[i])-upperF(x=exp(x),parR=parR,log=TRUE))
      }
      logval[i] <- uniroot(f=fun,interval=c(-3,3),extendInt='yes')$root	
    }
  }
  if(log){
    return( logval )	
  } else{
    return( exp(logval) )
  }
}

#generate n points from the point process R_1>R_2>R_3>... (in decreasing order) on [eps,infty)
rF <- function(N, parR){
  return(sort( upperFinv(y=N*runif(N),parR=parR,log=FALSE) ,decreasing = TRUE))
}

# dependence function for Gaussian processes
#parGauss = list(lambda,lambda.t,a,theta,nu,type),reg,reg.t
rho.func <- function(h,r=NULL,parGauss,reg=NULL,reg.t=NULL,mat=F){
  type = parGauss$type; lambda = parGauss$lambda; lambda.t = parGauss$lambda.t
  a = parGauss$a; theta = parGauss$theta; nu = parGauss$nu;
  if(type==1){
    rho.temp = h^2/lambda^2
    rho.r = exp(-sqrt(rho.temp))
    if(mat){return(matrix(c(1,rho.r,rho.r,1),nrow=2,byrow=T))}
    else{return(rho.r)}
  }
  if(type==2){
    if(is.null(reg.t)) {lambda.hat = reg %*%  lambda}else{lambda.hat = reg %*%  lambda + c(reg.t %*% lambda.t)}
    rho.temp = sum(lambda.hat)-log((exp(2*lambda.hat[1])+exp(2*lambda.hat[2]))/2) - sqrt(2*h^2/(exp(2*lambda.hat[1])+exp(2*lambda.hat[2])))
    rho.r = exp(rho.temp)
    if(mat){return(matrix(c(1,rho.r,rho.r,1),nrow=2,byrow=T))}
    else{return(rho.r)}
  }
  if(type==3){
    U.inv = matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),ncol=2,byrow=T)
    A.inv =t(U.inv)%*%diag(c(1,1/a))%*%U.inv/(lambda^2)
    Q = t(h)%*%A.inv%*%h
    rho.r = exp(-sqrt(Q))
    if(mat){return(matrix(c(1,rho.r,rho.r,1),nrow=2,byrow=T))}
    else{return(rho.r)}
  }
  if(type==4){
    if(is.null(reg.t)) {lambda.hat = reg %*%  lambda}else{lambda.hat = reg %*%  lambda + c(reg.t %*% lambda.t)}
    rho.temp = sum(lambda.hat)-log((exp(2*lambda.hat[1])+exp(2*lambda.hat[2]))/2) - sqrt( 2*h^2/(exp(2*lambda.hat[1])+exp(2*lambda.hat[2]))*((1+r)^nu))
    rho.r = exp(rho.temp)
    if(mat){return(matrix(c(1,rho.r,rho.r,1),nrow=2,byrow=T))
    }else{return(rho.r)}
  }
  if(type==5){
    if(is.null(reg.t)) {lambda.hat = reg %*%  lambda}else{lambda.hat = reg %*%  lambda + c(reg.t %*% lambda.t)}
    U.inv = matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),ncol=2,byrow=T)
    A.inv =2*t(U.inv)%*%diag(c(1,1/a))%*%U.inv/(exp(2*lambda[1])+exp(2*lambda[2]))*(1+r)^nu
    Q = t(h)%*%A.inv%*%h
    rho.temp = sum(lambda.hat)-log((exp(2*lambda.hat[1])+exp(2*lambda.hat[2]))/2)-sqrt(Q)
    rho.r = exp(rho.temp)
    if(mat){return(matrix(c(1,rho.r,rho.r,1),nrow=2,byrow=T))}
    else{return(rho.r)}
  }
}

#parGauss is a list
fun2int <- function(r,xi,h,parR,parGauss,reg=NULL,reg.t=NULL){
  xi.list <- rep(list(xi),length(r))
  r.list <- as.list(r)
  sigma.list <- mapply(rho.func,r=r,
                       MoreArgs = list(h=h,parGauss=parGauss,reg=reg,reg.t=reg.t,mat=T),SIMPLIFY = F)
  logpmv <- log( pmax(1-mapply(function(xi,r,sigma){ return(pmvnorm(upper=sign(xi)*exp(log(abs(xi))-log(r)),corr=sigma)[1])},xi=xi.list,r=r.list,sigma=sigma.list),0) )
  return( exp( logpmv + dF(r,parR,log=TRUE) ) )
}

#################################################################################################################################################
### function V and its derivatives
#################################################################################################################################################

#density function of the point process (also -V12)
mV12 <- function(x,h,parR,parGauss,reg=NULL,reg.t=NULL,log=FALSE){
  if(!is.matrix(x)){
    x <- matrix(x,nrow=1)
  } # take matrix as argument
  dGi <- function(xi){
    fun <- function(r,h,parR,parGauss,reg=NULL,reg.t=NULL){
      rho <- mapply(rho.func,r=r,
                    MoreArgs = list(parGauss=parGauss,h=h,reg=reg,reg.t=reg.t,mat=F),SIMPLIFY = T)
      v <- exp( -log(2*pi)-0.5*log(1-rho^2)-((sign(xi[1])*exp(log(abs(xi[1]))-log(r)))^2-2*rho*(sign(xi[1])*exp(log(abs(xi[1]))-log(r)))*(sign(xi[2])*exp(log(abs(xi[2]))-log(r)))+(sign(xi[2])*exp(log(abs(xi[2]))-log(r)))^2)/(2*(1-rho^2))-2*log(r)+dF(r,parR,log=TRUE) ) 
      return(v)
    }
    val <- integrate(fun,lower=0,upper=Inf,h=h,reg.t=reg.t,reg=reg,parR=parR,parGauss=parGauss,rel.tol=10^(-4),stop.on.error=FALSE)$value
    return(val)
  }
  I <- apply(is.na(x),1,sum)==0
  val <- c()
  val[I] <- apply(matrix(x[I,],nrow=sum(I)),1,dGi)
  val[!I] <- NA
  if(log){
    return( log(val) )	
  } else{
    return( val )
  }
}

#Partial derivatives of -V with respect to the k=1,2 element
mVk <- function(x,k,h,parR,parGauss,reg=NULL,reg.t=NULL,log=FALSE){
  if(!is.list(x)){
    if(!is.matrix(x)){
      x <- matrix(x,nrow=1)
    }
    x <- as.list(data.frame(t(x)))
  }# take argument as list
  #k is the vector that contains the index for partial derivatives: I will make k a list
  if(!is.list(k)){
    k <- as.list(k)
  }
  
  g <- function(xi,ki){
    #function of r to be integrated (needs to be defined for different values of r (= r is a vector))
    fun <- function(r,h,parR,parGauss,reg,reg.t){
      rho <- mapply(rho.func,r=r,
                    MoreArgs = list(parGauss=parGauss,h=h,reg=reg,reg.t=reg.t,mat=F),SIMPLIFY = T)
      return( exp( pnorm(sign(xi[-ki]-rho*xi[ki])*exp(log(abs(xi[-ki]-rho*xi[ki]))-log(r)),mean=0,sd=sqrt(1-rho^2),log.p=TRUE) + dnorm(sign(xi[ki])*exp(log(abs(xi[ki]))-log(r)),log=TRUE) -log(r) + dF(r,parR,log=TRUE) ) )
    }
    val <- integrate(fun,lower=0,upper=Inf,h=h,reg=reg,reg.t=reg.t,parR=parR,parGauss=parGauss,rel.tol=10^(-8),stop.on.error=FALSE)$value
    return(val)
  }
  val <- c()
  I <- mapply(function(x) sum(is.na(x))==0,x)
  val[I] <- mapply(g,xi=x[I],ki=k)
  val[!I] <- NA
  
  if(log){
    return( log(val) )	
  } else{
    return( val )
  }
}

# V (integral of point process density over outer region)
V <- function(x,h,parR,parGauss,reg=NULL,reg.t=NULL,log=FALSE){
  if(!is.list(x)){
    if(!is.matrix(x)){
      x <- matrix(x,nrow=1)
    }
    x <- as.list(data.frame(t(x)))
  }
  
  g <- function(xi){
    if (any(xi<=0)){return(Inf)}
    #function of r to be integrated (needs to be defined for different values of r (= r is a vector))
    fun <- function(r,h,parR,parGauss,reg=NULL,reg.t=NULL){
      xi.list <- rep(list(xi),length(r))
      r.list <- as.list(r)
      sigma.list <- mapply(rho.func,r=r,
                           MoreArgs = list(parGauss=parGauss,h=h,reg=reg,reg.t=reg.t,mat=T),SIMPLIFY = F)
      logpmv <- log( pmax(1-mapply(function(xi,r,sigma){ return(pmvnorm(upper=sign(xi)*exp(log(abs(xi))-log(r)),corr=sigma)[1]) },xi=xi.list,r=r.list,sigma=sigma.list),0) )
      return( exp( logpmv + dF(r,parR,log=TRUE) ) )
    }
    val <- integrate(fun,lower=0,upper=Inf,h=h,reg=reg,reg.t=reg.t,parR=parR,parGauss=parGauss,rel.tol=10^(-8),stop.on.error=FALSE)$value
    return(val)
  }
  val <- c()
  I <- mapply(function(x) sum(is.na(x))==0,x)
  val[I] <- mapply(g,xi=x[I])
  val[!I] <- NA
  
  if(log){
    return( log(val) )	
  } else{
    return( val )
  }
}


# Multivariate exponent function V (integral of point process density over outer region)
V.dimD <- function(x,h,parR,parGauss,reg=NULL,reg.t=NULL,log=FALSE){
  D = 2;
  if(is.matrix(h)){D=ncol(h)}
  if(!is.list(x)){
    if(!is.matrix(x)){
      x <- matrix(x,nrow=length(x),ncol=D)
    }
    x <- as.list(data.frame(t(x)))
  }
  g <- function(xi){
    if (any(xi<=0)){return(Inf)}
    #function of r to be integrated (needs to be defined for different values of r (= r is a vector))
    fun <- function(r,h,parR,parGauss,reg=NULL,reg.t=NULL){
      xi.list <- rep(list(xi),length(r))
      r.list <- as.list(r)
      corr<-function(r){
        if(is.matrix(h)){
          D = ncol(h);
          pairs<-as.data.frame(combn(D,2))
          reg.list = sapply(as.list(pairs),function(x){return(reg[x,])},simplify = F)
          h.list = h[t(pairs)] 
          cor.val <- diag(1,D)
          cor.val[t(pairs)]<-cor.val[t(pairs[2:1,])] <-mapply(rho.func,h=h.list,reg=reg.list,MoreArgs = list(parGauss=parGauss,r=r,reg.t=reg.t,mat=F),SIMPLIFY = T)
          return(cor.val)
        }else{
          return(rho.func(h=h,r=r,parGauss = parGauss,reg = reg,reg.t=reg.t,mat = T))
        }
      }
      sigma.list <- mapply(corr,r=r,SIMPLIFY = F)
      temp = function(xi,r,sigma){ 
        oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
        set.seed(19873436)
        res <- pmvnorm(upper=sign(xi)*exp(log(abs(xi))-log(r)),corr=sigma)[1]
        assign(".Random.seed", oldSeed, envir=globalenv())
        return(res) 
      }
      logpmv <- log( pmax(1-mapply(temp,xi=xi.list,r=r.list,sigma=sigma.list),0) )
      return( exp( logpmv + dF(r,parR,log=TRUE) ) )
    }
    val <- integrate(fun,lower=0,upper=Inf,h=h,reg=reg,reg.t=reg.t,parR=parR,parGauss=parGauss,rel.tol=10^(-8),stop.on.error=FALSE)$value
    return(val)
  }
  val <- c()
  I <- mapply(function(x) sum(is.na(x))==0,x)
  val[I] <- mapply(g,xi=x[I])
  val[!I] <- NA
  
  if(log){
    return( log(val) )	
  } else{
    return( val )
  }
}

# Multivariate level-dependent extremal coefficient: take unit Frechet scale
Theta.dimD <- function(z,h,parR,parGauss,reg=NULL,reg.t=NULL){
  z.RW <- qG(exp(-1/z),parR=parR,log=FALSE)
  return( z*V.dimD(z.RW,h=h,parR=parR,parGauss=parGauss,reg=reg,reg.t = reg.t,log=FALSE) )
} 

#################################################################################################################################################
### marginal distribution of the max-id model, its density and quantile function: pG, dG, qG
#################################################################################################################################################
pG <- function(x,parR,log=FALSE){
  g <- function(xi){
    if (xi<=0){return(Inf)}
    #function of r to be integrated (needs to be defined for different values of r (= r is a vector))
    fun <- function(r,parR){
      logp <- log( 1-pnorm(sign(xi)*exp(log(abs(xi))-log(r))) )
      return( exp( logp + dF(r,parR,log=TRUE) ) )
    }
    val <- integrate(fun,lower=0,upper=Inf,parR=parR,rel.tol=10^(-8),stop.on.error=FALSE)$value
    return(val)
  }
  val <- c()
  I <- !is.na(x) 
  val[I] <- exp(-apply(matrix(x[I],ncol=1),1,g))
  val[!I] <- NA
  logval <- log(val)
  
  if(log){
    return( logval )	
  } else{
    return( exp(logval) )
  }	
}

dG <- function(x,parR,log=FALSE){
  g <- function(xi){
    if (xi<=0){return(0)}
    fun <- function(r,parR){
      return(exp(dnorm(sign(xi)*exp(log(abs(xi))-log(r)),log=TRUE)-log(r)+dF(r,parR,log=TRUE)))
    }
    val <- integrate(fun,lower=0,upper=Inf,parR=parR,rel.tol=10^(-8),stop.on.error=FALSE)$value
    return(val)
  }
  I <- !is.na(x) 
  val0 <- c()
  val0[I] <- apply(matrix(x[I],ncol=1),1,g)
  val0[!I] <- NA
  logval <- log(val0)+pG(x,parR,log=TRUE)  #density is -V_1 exp(-V)
  if(log){
    return( logval )
  } else{
    return( exp(logval) )
  }	
}

qG <- function(p,parR,log=FALSE){
  fun <- function(x,p,parR){
    return( pG(exp(x),parR,log=TRUE)-log(p) )
  }
  I <- !is.na(p) & (p>0) & (p<1)
  logval <- c()
  logval[!is.na(p) & p==0] <- -Inf
  logval[!is.na(p) & p==1] <- Inf
  logval[I] <- apply(matrix(p[I],ncol=1),1,function(pi) uniroot(fun,interval=c(-3,3),p=pi,parR=parR,extendInt='yes')$root)
  if(log){
    return( logval )
  } else{
    return( exp(logval) )
  }
}


##################################################################################################################################
### spatial simulations of max-id model (multivariate based on a stable correlation function, using coordinates of stations) 
##################################################################################################################################
# if type=1,2,4 then coord should be the distance matrix, otherwise its the coordinates
rmaxidspat <- function(n, coord, parR, parGauss,reg=NULL,reg.t=NULL, N=1000, ncores=1){
  D <- nrow(coord)
  Z <- matrix(nrow=n,ncol=D)
  pairs <- combn(1:D,2)
  Sigma <- diag(1,D)
  if(is.null(parGauss$nu)){
    fun <- function(pair,t=NULL){
      if(parGauss$type %in% c(1,2,4)){
        h = coord[pair[1],pair[2]]
      }else{h = abs(coord[pair[1],]-coord[pair[2],])}
      reg.new=NULL;reg.t.new=NULL
      if(!is.null(reg)){reg.new=reg[pair,]}
      if(!is.null(reg.t)){reg.t.new=reg.t[t,]}
      rho <- rho.func(h=h,r=NULL,parGauss =parGauss,reg=reg.new,reg.t=reg.t.new,mat=F)
    }
    if(is.null(reg.t)){
      Sigma[t(pairs)]<-mcmapply(fun, pair=as.list(as.data.frame(pairs)),mc.cores = ncores)
      Sigma[t(pairs)[,2:1]]<-Sigma[t(pairs)]
      A <- t(chol(Sigma)) }
    for (k in 1:n){
      if(!is.null(reg.t)){
        Sigma[t(pairs)]<-mcmapply(fun, pair=as.list(as.data.frame(pairs)),MoreArgs = list(t=k),mc.cores = ncores)
        Sigma[t(pairs)[,2:1]]<-Sigma[t(pairs)]
        A <- t(chol(Sigma))   
      }
      R <- rF(N,parR)
      Z[k,] <- apply(R*t(A%*%matrix(rnorm(D*N),ncol=N,nrow=D)),2,max)
    }
    return(Z)
  } else{	
    for (k in 1:n){
      R <- rF(N,parR)
      tmp <- function(r,t=k){
        fun <- function(pair){
          if(parGauss$type %in% c(1,2,4)){
            h = coord[pair[1],pair[2]]
          }else{h = abs(coord[pair[1],]-coord[pair[2],])}
          if(!is.null(reg)){reg.new=reg[pair,]}else{reg.new=NULL}
          if(!is.null(reg.t)){reg.t.new=reg.t[t,]}else{reg.t.new=NULL}
          rho <- rho.func(h=h,r=r,parGauss=parGauss,reg=reg.new,reg.t=reg.t.new,mat=F)
        }
        Sigma[t(pairs)]<-mcmapply(fun, pair=as.list(as.data.frame(pairs)),mc.cores = ncores)
        Sigma[t(pairs)[,2:1]]<-Sigma[t(pairs)]
        A <- t(chol(Sigma))
        return( r*t(A%*%matrix(rnorm(D),nrow=D,ncol=1)) )
      }
      Z[k,] <- apply(matrix(unlist(mclapply(R,tmp,mc.cores=ncores)),nrow=D,ncol=N),1,max)
    }
    return(Z)
  }
}


#####################################
### Pairwise likelihood inference ###
#####################################
### These functions were used for the data application

#pairwise copula likelihood computed in parallel
pw.nllik.parallel <- function(parR,parGauss,datU,new.pair,reg=NULL,reg.t=NULL,coord,doSum=TRUE,print.par.file=NULL,ncores=1){ #negative pairwise log likelihood (for the copula)
  if(parGauss$type==1){
    if (parR[1]<=0  |parR[1]>10 |parR[2]<0 | parR[2]>20 | parGauss$lambda<=0){return(Inf)} 
  }
  if(parGauss$type==2){
    if (parR[1]<=0  |parR[1]>10 |parR[2]<0 | parR[2]>20 ){return(Inf)} 
  }
  if(parGauss$type==3){
    if (parR[1]<=0  |parR[1]>10 |parR[2]<0 | parR[2]>20  | parGauss$lambda<=0 | parGauss$theta<0 | parGauss$theta > pi | parGauss$a<0){return(Inf)} 
  }
  if(parGauss$type==4){
    if (parR[1]<=0  |parR[1]>10 |parR[2]<0 | parR[2]>20  ){return(Inf)} 
  }
  if(parGauss$type==5){
    if (parR[1]<=0  |parR[1]>10 |parR[2]<0 | parR[2]>20  | parGauss$theta < 0 | parGauss$theta > pi | parGauss$a<0){return(Inf)} 
  }
  if(is.null(reg.t)){reg.t<-NULL}else{reg.t<-as.matrix(reg.t)}
  makeXDAT <- function(k){ #k is index of column
    A0 <- !is.na(datU[,k])
    res <- rep(NA,nrow(datU))
    if(sum(A0)==0){return(res)}
    res[A0] <- qG(datU[A0,k],parR)
    return(res)
  }
  
  contrib.pair <- function(k){ #k is index of pair
    dat0 <- XDAT[,new.pair[k,]]
    A0 <- (!is.na(dat0[,1])) & (!is.na(dat0[,2])) # time without missing value for the pair
    n0 <- sum(A0)
    if (n0==0) return(rep(NA,nrow(XDAT)))
    v <- c();v[!A0] <- NA
    if(!is.null(reg)){reg.new=reg[new.pair[k,],]}else{reg.new=NULL}
    if(parGauss$type %in% c(1,2,4)){
      h = coord[new.pair[k,1],new.pair[k,2]]
    }else{h <- abs(coord[new.pair[k,1],]-coord[new.pair[k,2],])}
    if(is.null(reg.t)){
      v[A0] <- V(XDAT[A0,new.pair[k,]],h=h,parR,parGauss,reg=reg.new)-log( mVk(XDAT[A0,new.pair[k,]],k=rep(1,n0),h=h,parR,parGauss,reg=reg.new)*mVk(XDAT[A0,new.pair[k,]],k=rep(2,n0),h=h,parR,parGauss,reg=reg.new) + mV12(XDAT[A0,new.pair[k,]],h=h,parR,parGauss,reg=reg.new) ) + apply(apply(XDAT[A0,new.pair[k,]],c(1,2),dG,parR=parR,log=TRUE),1,sum) #EMERIC: I think this is OK (when missing values I remove the pair from the pairwise likelihood)
    }else{
      t = which(A0==TRUE)
        fun <- function(i){
          V(XDAT[i,new.pair[k,]],h=h,parR,parGauss,reg=reg.new,reg.t=reg.t[i,])-log( mVk(XDAT[i,new.pair[k,]],k=1,h=h,parR,parGauss,reg=reg.new,reg.t=reg.t[i,])*mVk(XDAT[i,new.pair[k,]],k=2,h=h,parR,parGauss,reg=reg.new,reg.t=reg.t[i,]) + mV12(XDAT[i,new.pair[k,]],h=h,parR,parGauss,reg=reg.new,reg.t=reg.t[i,])) + sum(dG(XDAT[i,new.pair[k,]],parR=parR,log=TRUE))#EMERIC: I think this is OK (when missing values I remove the pair from the pairwise likelihood)
        }
    v[A0] <- mapply(fun,t,SIMPLIFY = T)
    return(v)
  }}
  #marginal transformation
  oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
  set.seed(19873436)
  XDAT <- matrix(unlist(mclapply(1:ncol(datU),makeXDAT,mc.cores=ncores)),ncol=ncol(datU))
  #pairwise log-likelihood contributions
  val <- matrix(unlist(mclapply(1:nrow(new.pair),contrib.pair,mc.cores=ncores)),ncol=nrow(new.pair))
  #val <- matrix(unlist(mapply(contrib.pair,1:nrow(new.pair))),ncol=nrow(new.pair))
  assign(".Random.seed", oldSeed, envir=globalenv())
  if(!is.null(print.par.file)){
    cat(c(par,sum(val),"\n"),file=print.par.file,append=TRUE)
  }
  
  if(doSum){
    return(sum(val,na.rm=T)) ### sums all contributions over replicates and pairs
  } else{
    return(apply(val,1,sum,na.rm=T)) ### sums all contributions over pairs, but NOT over replicates
  }
}

# separate the parameters into two parts.
# including one spatial covariate and one temporal covariate
get.par <- function(par,type=4){
  parR <- exp(par[1:2])
  parGauss<-list(type=type,lambda=NULL,lambda.t=NULL,a=NULL,theta=NULL,nu=NULL)
  if(type==1){parGauss$lambda=exp(par[3])}
  if(type==2){parGauss$lambda=par[3:4];parGauss$lambda.t=par[5]}
  if(type==3){parGauss$lambda=exp(par[3]);parGauss$a=exp(par[4]);parGauss$theta=par[5]}
  if(type==4){parGauss$nu=par[6];parGauss$lambda=par[3:4];parGauss$lambda.t=par[5]}
  if(type==5){parGauss$a=exp(par[3]);parGauss$theta=par[4];parGauss$nu=par[5];parGauss$lambda=par[6:7];parGauss$lambda.t=par[8]}
  return(list(parR=parR,parGauss=parGauss))
} 

#fit the copula model using pairwise likelihood in parallel
fit.pw.parallel <- function(init,datU,coord,reg=NULL,reg.t=NULL,type=4,cutoff=Inf,proppairs=1,fixed=rep(FALSE,3),optim=TRUE,gd=F,hessian=TRUE,sandwich=TRUE,eps=10^(-2),print.par.file=NULL,ncores=1,fit.save=FALSE,fit.load=FALSE,fit.file=NULL,maxit=300,...){ #fit max-id copula model by negative pairwise log likelihood
  parGauss = list(lambda=NULL,lambda.t=NULL,a=NULL,theta=NULL,nu=NULL,type=type)
  #pairs and corresponding rho: put everything in a matrix to which I will apply the likelihood of bivariate data
  D <- nrow(coord)
  pair <- expand.grid(1:D,1:D)
  pair <- pair[,2:1]
  pair <- pair[pair[,1]<pair[,2],]
  pair <- matrix(c(pair[[1]],pair[[2]]),ncol=2)
  if(parGauss$type %in% c(1,2,4)){dist.pair <- coord[pair]}else{dist.pair <- as.matrix(dist(coord))[pair]}
  ind <- sort(which(dist.pair<cutoff))
  new.pair <- pair[ind,]
  if(optim==TRUE){
    #par=c(alpha,beta,lambda,a,theta,nu)
    pw.nllik2 <- function(par,pair2,index=nrow(datU)){
      par2 <- init
      par2[which(!fixed)] <- par
      par2.list <- get.par(par2,type)
      if(is.null(reg.t)){reg.t2=NULL}else{reg.t2=reg.t[1:index,]}
      val <- pw.nllik.parallel(par2.list$parR,par2.list$parGauss,datU[1:index,],new.pair=pair2,reg=reg,reg.t=reg.t2,coord,doSum=TRUE,print.par.file=print.par.file,ncores=ncores)
      return(val)
    }
    init2 <- init[which(!fixed)]
    if(fit.load){
      load(file=fit.file)
    } else{
      init2 <- init[which(!fixed)]
      if(!gd){
        fit <- optim(init2,pw.nllik2,pair2=new.pair,control = list(trace=1,maxit=maxit),hessian=hessian,...)
      }else{
        pair.tmp = new.pair[sort(sample(1:nrow(new.pair),ceiling(nrow(new.pair)*proppairs))),]
        fit <- optim(init2,pw.nllik2,pair2=pair.tmp,method="L-BFGS-B",lower=rep(-Inf,sum(!fixed)),upper=c(log(10),log(10),rep(Inf,sum(!fixed)-2)),control=list(trace=1),hessian = hessian,...)
      }
      if(fit.save){
        save(fit,file=fit.file)
      }
    }
    
    res <- list()
    mle <- c(); mle[!fixed] <- fit$par; mle[fixed] <- init[fixed]; res$par <- mle
    res$pw.nllik <- fit$val
    res$convergence <- fit$convergence
    res$counts <- fit$counts
    if(sandwich & hessian){
      hessian.pw.nllik <- function(par,datU,coord){
        par2 <- init
        par2[which(!fixed)] <- par
        hess.mat <- matrix(nrow=length(init),ncol=length(init))
        
        for(i in which(!fixed)){
          for(j in which(!fixed)){
            if(j==i){
              par2exp <- par2
              par2p <- par2; par2p[i] <- par2p[i]+par2[i]*(2*eps)
              par2m <- par2; par2m[i] <- par2m[i]-par2[i]*(2*eps)
              hess.mat[i,j] <- (pw.nllik2(par2p,new.pair)-2*pw.nllik2(par2exp,new.pair)+pw.nllik2(par2m,new.pair))/(4*eps^2*par2[i]^2)
            } else if(j>i){
              par2pp <- par2; par2pp[i] <- par2pp[i]+par2[i]*eps
              par2pm <- par2; par2pm[i] <- par2pm[i]+par2[i]*eps
              par2mp <- par2; par2mp[i] <- par2mp[i]-par2[i]*eps
              par2mm <- par2; par2mm[i] <- par2mm[i]-par2[i]*eps
              hess.mat[i,j] <- hess.mat[j,i] <- (pw.nllik2(par2pp,new.pair)-pw.nllik2(par2pm,new.pair)-pw.nllik2(par2mp,new.pair)+pw.nllik2(par2mm,new.pair))/(4*eps^2*par2[i]*par2[j])
            }
          }
        }
        hess.mat <- hess.mat[which(!fixed),which(!fixed)]
        return(hess.mat)
      }
      #res$hessian <- hessian.pw.nllik(fit$par,datU,coord)
      res$hessian <- fit$hessian
      res$inv.hess <- solve(res$hessian)
      grad.pw.nllik <- function(par,datU,coord){
        par2 <- init
        par2[which(!fixed)] <- par
        grad.mat <- matrix(nrow=length(init),ncol=nrow(datU)-9)
        for(i in which(!fixed)){
          par2p <- par2; par2p[i] <- par2p[i]+par2[i]*eps
          par2m <- par2; par2m[i] <- par2m[i]-par2[i]*eps
          for(j in 10:nrow(datU)){
          grad.mat[i,j-9] <- (pw.nllik2(par2p,new.pair,index = j)-pw.nllik2(par2m,new.pair,index = j))/(2*eps*par2[i])}
          message(grad.mat[i,j-9])
          }
        grad.mat <- grad.mat[which(!fixed),]
        return(grad.mat)
      }
      res$grad <- grad.pw.nllik(fit$par,datU,coord)
      res$var.grad <- cov(t(res$grad))
      res$sandwich <- res$inv.hess%*%res$var.grad%*%res$inv.hess
      res$sd.mle <- sqrt(diag(res$sandwich))
      res$clic <- 2*res$pw.nllik+2*sum(diag(res$inv.hess%*%res$var.grad))
      D <- nrow(coord)
      pair <- expand.grid(1:D,1:D)
      pair <- pair[,2:1]
      pair <- pair[pair[,1]<pair[,2],]
      pair <- matrix(c(pair[[1]],pair[[2]]),ncol=2)
      dist.pair <- as.matrix(dist(coord))[pair]
      npairs <- new.pair
      npairs <- sum(dist.pair<cutoff)
      res$clic.star <- res$clic/(2*npairs/D)
    }
    return(res)
  } else{
    par2 <- init
    pw.nllik2(par2,new.pair)	
  }
}

########################################################
##### Paramertric Copula Model           ###############
########################################################
pw.copula.parallel <- function(parGauss,datU,new.pair,reg=NULL,reg.t=NULL,coord,doSum=TRUE,copula.type=1,log.df=5,print.par.file=NULL,ncores=1){ #negative pairwise log likelihood (for the copula)
  if(parGauss$type==1){
    if (parGauss$lambda<=0){return(Inf)} 
  }
  if(parGauss$type==3){
    if ( parGauss$lambda<=0 | parGauss$theta<0 | parGauss$theta > pi | parGauss$a<0){return(Inf)} 
  }
  if(parGauss$type==5){
    if (parGauss$theta < 0 | parGauss$theta > pi | parGauss$a<0){return(Inf)} 
  }
  if(is.null(reg.t)){reg.t<-NULL}else{reg.t<-as.matrix(reg.t)}
  
  dcop <- function(rho,y){
    if(copula.type==1){
      cop <- normalCopula(rho,dim=2)
    }else{
      cop <- tCopula(rho,dim=2,df=ceiling(exp(log.df)))
    }
    val<- -dCopula(y,cop,log=T)
    return(val)
  }
   
  contrib.pair <- function(k){ #k is index of pair
    dat0 <- datU[,new.pair[k,]]
    A0 <- (!is.na(dat0[,1])) & (!is.na(dat0[,2])) # time without missing value for the pair
    n0 <- sum(A0)
    if (n0==0) return(rep(NA,nrow(datU)))
    v <- c();v[!A0] <- NA
    if(!is.null(reg)){reg.new=reg[new.pair[k,],]}else{reg.new=NULL}
    if(parGauss$type %in% c(1,2,4)){
      h = coord[new.pair[k,1],new.pair[k,2]]
    }else{h <- abs(coord[new.pair[k,1],]-coord[new.pair[k,2],])}
    if(is.null(reg.t)){
      rho <- rho.func(h=h,parGauss = parGauss,reg = reg.new)
      v[A0] <- dcop(rho,datU[,new.pair[k,]])
    }else{
      rho <- mapply(rho.func,reg.t=reg.t[A0,],MoreArgs = list(parGauss=parGauss,reg=reg.new,mat=F,h=h),SIMPLIFY = T)
      
      v[A0] <- mapply(dcop,rho=rho,y=as.list(as.data.frame(t(datU[A0,new.pair[k,]]))),SIMPLIFY = T)
      return(v)
    }}
    #pairwise log-likelihood contributions
    val <- matrix(unlist(mclapply(1:nrow(new.pair),contrib.pair,mc.cores=ncores)),ncol=nrow(new.pair))
    #val <- matrix(unlist(mapply(contrib.pair,1:nrow(new.pair))),ncol=nrow(new.pair))
    if(!is.null(print.par.file)){
      cat(c(par,sum(val),"\n"),file=print.par.file,append=TRUE)
    }
    
    if(doSum){
      return(sum(val,na.rm=T)) ### sums all contributions over replicates and pairs
    } else{
      return(apply(val,1,sum,na.rm=T)) ### sums all contributions over pairs, but NOT over replicates
    }
  }

## lambda_0,lambda_s,lambda_t,log.df
fit.copula.parallel <- function(init,datU,coord,reg=NULL,reg.t=NULL,type=1,cutoff=Inf,proppairs=1,copula.type=1,fixed=rep(FALSE,4),gd=F,eps=10^(-2),ncores=1,...){ 
  #pairs and corresponding rho: put everything in a matrix to which I will apply the likelihood of bivariate data
  D <- nrow(coord)
  pair <- expand.grid(1:D,1:D)
  pair <- pair[,2:1]
  pair <- pair[pair[,1]<pair[,2],]
  pair <- matrix(c(pair[[1]],pair[[2]]),ncol=2)
  if(type %in% c(1,2,4)){dist.pair<-coord[pair]}else{dist.pair <- as.matrix(dist(coord))[pair]}
  ind <- sort(which(dist.pair<cutoff))
  new.pair <- pair[ind,]
  pw.cop2 <- function(par,ind){
    par2 <- init
    par2[which(!fixed)] <- par
    parGauss <- get.par(c(0,0,par2[1:3]),2)$parGauss
    log.df = par2[4]
    val <- pw.copula.parallel(parGauss=parGauss,datU=datU,new.pair=new.pair,reg=reg,reg.t=reg.t,coord=coord,ncores=ncores,copula.type=copula.type,log.df = log.df,doSum=TRUE)
    return(val)
  }
  init2 <- init[which(!fixed)]
  if(!gd){
    fit <- optim(init2,pw.cop2,ind=1:nrow(new.pair),control = list(trace=1),hessian=FALSE,...)
  }else{
    ind.tmp =sort(sample(1:nrow(new.pair),ceiling(nrow(new.pair)*proppairs))) 
    fit <- optim(init2,pw.cop2,ind=ind.tmp,method="BFGS",control=list(trace=1),hessian = F,...)
  }
  res <- list()
  mle <- c(); mle[!fixed] <- fit$par; mle[fixed] <- init[fixed]; res$par <- mle
  res$pw.nllik <- fit$val
  res$convergence <- fit$convergence
  res$counts <- fit$counts
  return(res)
}

########################################################
##### Extremal Coefficients Least Square ###############
########################################################

#fit the copula model using pairwise extremal coefficients in parallel
## In this model, we can not include temporal covariates
fit.extremal.parallel <- function(init,datU,coord,reg=NULL,type=1,cutoff=Inf,proppairs=1,q=seq(0.05,0.95,0.1),fixed=rep(FALSE,5),gd=F,eps=10^(-2),ncores=1,...){ 
  #pairs and corresponding rho: put everything in a matrix to which I will apply the likelihood of bivariate data
  D <- nrow(coord)
  pair <- expand.grid(1:D,1:D)
  pair <- pair[,2:1]
  pair <- pair[pair[,1]<pair[,2],]
  pair <- matrix(c(pair[[1]],pair[[2]]),ncol=2)
  if(type %in% c(1,2,4)){dist.pair<-coord[pair]}else{dist.pair <- as.matrix(dist(coord))[pair]}
  ind <- sort(which(dist.pair<cutoff))
  new.pair <- pair[ind,]
  #pairwise extremal coefficients computed in parallel
  z = qfrechet(q,loc=0,scale=1,shape=1)
  makeXDAT <- function(k){ #k is index of column
    A0 <- !is.na(datU[,k])
    res <- rep(NA,nrow(datU))
    res[A0] <- qfrechet(datU[A0,k],loc=0,scale=1,shape=1)
    return(res)
  }
  #marginal transformation to unit frechet
  XDAT <- matrix(unlist(mclapply(1:ncol(datU),makeXDAT,mc.cores=ncores)),ncol=ncol(datU))
  emp.pair <- function(k){# index of the pairs
    dat0 <- XDAT[,new.pair[k,]]
    A0 <- (!is.na(dat0[,1])) & (!is.na(dat0[,2])) # time without missing value for the pair
    n0 <- sum(A0)
    if (n0==0) return(rep(NA,nrow(XDAT)))
    pf <- ecdf(pmax(dat0[A0,1],dat0[A0,2]))
    v <- -z*log(pf(z))
    return(v)
  }
  theta.emp <- matrix(unlist(mclapply(1:nrow(new.pair),emp.pair,mc.cores=ncores)),ncol=nrow(new.pair))
  pw.extremal.parallel <- function(parR,parGauss,ind,doSum=TRUE){ 
    new.pair2 = new.pair[ind,]
    if(parGauss$type==1){
      if (parR[1]<=0 | parR[1]>20 | parR[2]<0 | parGauss$lambda<=0){return(Inf)} 
    }
    if(parGauss$type==2){
      if (parR[1]<=0 |  parR[2]<0 ){return(Inf)} 
    }
    if(parGauss$type==3){
      if (parR[1]<=0 | parR[2]<0  | parGauss$lambda<=0 | parGauss$theta<0 | parGauss$theta > pi | parGauss$a<0){return(Inf)} 
    }
    if(parGauss$type==4){
      if (parR[1]<=0 | parR[2]<0 ){return(Inf)} 
    }
    if(parGauss$type==5){
      if (parR[1]<=0 | parR[2]<0 | parGauss$theta < 0 | parGauss$theta > pi  | parGauss$a<0){return(Inf)} 
    }
    
    contrib.pair <- function(k){ #k is index of pair
      if(!is.null(reg)){reg.new=reg[new.pair2[k,],]}else{reg.new=NULL}
      if(type %in% c(1,2,4)){
        h = coord[new.pair2[k,1],new.pair2[k,2]]
      }else{h <- abs(coord[new.pair2[k,1],]-coord[new.pair2[k,2],])}
      v <- Theta.dimD(z,h,parR,parGauss,reg.new,reg.t=NULL)
      return(v)
    }
    #pairwise log-likelihood contributions
    val <- matrix(unlist(mclapply(1:nrow(new.pair2),contrib.pair,mc.cores=ncores)),ncol=nrow(new.pair2))
    ls <- (val-theta.emp[,ind])^2
    
    if(doSum){
      return(sum(ls)) ### sums all contributions over replicates and pairs
    } else{
      return(apply(ls,1,sum)) ### sums all contributions over pairs, but NOT over replicates
    }
  }
  
  #parGauss = list(lambda=NULL,lambda.t=NULL,a=NULL,theta=NULL,nu=NULL,type=type)
    #par=c(alpha,beta,lambda,a,theta,nu)
    pw.extremal2 <- function(par,ind){
      par2 <- init
      par2[which(!fixed)] <- par
      par2.list <- get.par(par2,type)
      val <- pw.extremal.parallel(par2.list$parR,par2.list$parGauss,ind=ind,doSum=TRUE)
      return(val)
    }
    init2 <- init[which(!fixed)]
    if(!gd){
      fit <- optim(init2,pw.extremal2,ind=1:nrow(new.pair),control = list(trace=1),hessian=FALSE,...)
    }else{
      ind.tmp =sort(sample(1:nrow(new.pair),ceiling(nrow(new.pair)*proppairs))) 
      fit <- optim(init2,pw.extremal2,ind=ind.tmp,method="BFGS",control=list(trace=1),hessian = F,...)
    }
    res <- list()
    mle <- c(); mle[!fixed] <- fit$par; mle[fixed] <- init[fixed]; res$par <- mle
    res$pw.nllik <- fit$val
    res$convergence <- fit$convergence
    res$counts <- fit$counts
    return(res)
}




