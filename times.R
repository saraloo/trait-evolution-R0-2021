
fun.R0 <- function(r=as.numeric(init.tr["r"]),del=as.numeric(init.tr["delta"])){
  
  # initial values
  delta <- del
  M0 <- as.numeric(init["M"])
  rP <- r
  P0 <- as.numeric(init["P"])
  c <- as.numeric(init.tr["c"])
  u <- as.numeric(init.tr["u"])
  
  tk <- (1/rP) * log(nu /P0)  # immune activation time
  tm <- tk + (1/rM)*log(KM/M0)
  
  # Gives solution to rhat (dividing weak from stronger pathogens)
  Pasterisk <- function(x){ 
    b <- (KM/M0 -1)*(nu/P0)^(rM/x)
    PP<- P0 *(((1+b)*(delta*KM-x)/(delta*KM*b))^(delta*KM/rM))*(b*x/(delta*KM-x))^(x/rM)
    return(PP-K)
  }
  rhat <-uniroot(Pasterisk,lower = 0.05, upper = 0.2, tol = 1e-9)$root 
  
  t1 <- (1/rP)*log(eta*(K-P0)/(P0*(K-eta)))
  
  if(rP < rhat){ ## weak pathogens
    b <- (KM/M0 -1)*(nu/P0)^(rM/rP)
    tP <- (1/rM)*log(b*rP/(delta*KM-rP))
    Past <- P0*(((1+b)*(delta*KM-rP)/(delta*KM*b))^(delta*KM/rM))*(b*rP/(delta*KM-rP))^(rP/rM) 
    t2 <- (1/(delta*KM-rP))*log((P0*(1+b)^(delta*KM/rM))/eta)
  }else{ ## strong pathogens
    t2 <- tm + (1/(rP-delta*KM))*log(eta/K)
    Past <- K
    tP <- tm
  }
  
  if(t2 < t1){
    t2<- t1
  }
  
  Delta <- t2 - t1
  alpha <- u * Past * Delta^1
  beta <- c * Past 
  
  chi <- 1 - exp(-alpha*Delta)
  
  R0 <- (beta/alpha)*(1-exp(-alpha * Delta))
  if(is.na(R0)){R0 <- 0}
  
  list(beta=beta,alpha=alpha,chi=chi,
       Delta=Delta,R0=R0,
       Past = Past,t1=t1,t2=t2, tP=tP,
       tk=tk,tm=tm,rhat=rhat)
  
}

Papprox <- function(t,r){
  if(r<rhat){  # weak regime
    b <- (KM/M0 -1)*(nu/P0)^(rM/r)
    P <- P0*(exp((r-delta*KM)*t))*((b+1)/(b*exp(-rM*t)+1))^(delta*KM/rM)
  }else{       # stronger regime
    if(t<tm){
      P <- K/(1+(K/P0 -1)*exp(-r*t))
    }else{
      P <- K*exp((r-delta*KM)*(t-tm))
    }
  }
  return(P)
}





