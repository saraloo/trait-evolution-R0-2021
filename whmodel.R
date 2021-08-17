set.p <- function(
  ## default parameters for within host model Loo & Tanaka 21
  K    = 10^9,     #  carrying capacity of pathogen P
  rM   = 0.1,      #  growth rate of immune system M
  KM   = 10^7,     #  carrying capacity of immune system M
  nu = 1000,       # immune threshold for activation
  eta = 10000,     # threshold for infectious period
  init.tr =  c(r=0.13,delta=3e-8,c=4e-12,u=2.8e-13) ,  # r=growth, delta=clearance by immune system, c=transmission constant, u=virulence constant
  tdiv = 1,        # divide time for discrete time simulation 
  T.end= 2500,     # time for end of simulation (WH)
  inoc.size=10,    # initial infection size
  init = c(P=inoc.size,M=10)  # initial conditions 
) {
  list(init.tr=init.tr/tdiv,
    K=K, rM=rM, KM=KM,
    nu = nu,eta = eta,
    tdiv=tdiv,
    T.end=T.end, 
    init=init
  )
}


wh.euler <- function(p = parameters){
  ## Function simulates within-host model by implementing Euler's method for ODE in Loo & Tanaka 21
  
  imm.fun <- function(t){ KM/(1+(KM/M0-1)*exp(-rM*(t-tk)*ifelse(t>tk,1,0)))}
  
  list2env(p,.GlobalEnv)
  r <- as.numeric(init.tr["r"])*tdiv
  delta <- as.numeric(init.tr["delta"])*tdiv
  P0 <- as.numeric(init["P"])
  M0 <- as.numeric(init["M"])
  
  times <- seq(0,T.end,by=1/tdiv) 

  tk <- (1/r) * log(nu /P0) # immune activation
  
  sim.out <- data.frame(t=times,M=c(M0,rep(NA,length(times)-1)),P=c(P0,rep(NA,length(times)-1)),row.names=NULL)   
  
  for(j in 2:length(times)){
    if(times[j]<tk){   ## immunity OFF
      imm <- 0
    }else{    ## immunity ON
      imm <- imm.fun(times[j])
    }
    sim.out$P[j] <- sim.out$P[j-1] + ( r*sim.out$P[j-1]*(1-sim.out$P[j-1]/K) - delta * imm * sim.out$P[j-1] )*(1/tdiv)
    sim.out$M[j] <- imm.fun(times[j])
  }
  
  
  return(sim.out)
  
}
