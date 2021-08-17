## This code runs within-host model and plots alongside approximate times, as in Figure 1 of Loo & Tanaka 21
library(ggplot2)
library(reshape2)
library(scales)

param <- set.p(init.tr=c(r=0.13,delta=3e-8,c=4.5e-12,u=5e-17),rM=0.1,T.end=20*24,tdiv=1)

list2env(param,.GlobalEnv)

r <- as.numeric(init.tr["r"])*tdiv
delta <- as.numeric(init.tr["delta"])*tdiv
P0 <- as.numeric(init["P"])
M0 <- as.numeric(init["M"])

soln.euler <- wh.euler(p = param)

sim.out <- soln.euler
imdf <- data.frame(t=sim.out$t,Im=sim.out$M)
# m.d <- sim.out[,-2]
m.d <- melt(sim.out,id="t")

R0times <- fun.R0(r=r,del=delta)
times <- c(t1=R0times$t1,tp=R0times$tP,t2=R0times$t2)
list2env(R0times,.GlobalEnv)
approxP <- data.frame(t=seq(0,20*24),P=sapply(seq(0,20*24),Papprox,r=r))

hlabels <- data.frame(x = c(22,22),y=c(Past,eta),label=c("K","eta"))
vlabels <- data.frame(y = c(K+5000,K+5000),x=c(times["t1"]/24,times["t2"]/24),label=c("t[1]","t[2]"))

plot1 <- ggplot(m.d)+
  geom_line(aes(t/24,value,colour=variable))+
  geom_line(data=approxP,aes(t/24,P),linetype="dashed",colour="black")+
  theme_classic() +
  theme(text = element_text(size=11),panel.background = element_rect(colour = "black"),
        aspect.ratio = 1)+
  scale_color_manual(values=c("black","#F8766D"))+
  labs(x = "Time (days)", y="Population size") +
  guides(colour=FALSE)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(10^0,10^10),expand=c(0,0))+
  coord_cartesian(expand=0,xlim = c(0,20),clip = 'off')+
  geom_text(data=hlabels,aes(x,y,label=label), parse = TRUE,hjust = 0.3)+ #,size = 6
  geom_text(data=vlabels,aes(x,y,label=label), parse = TRUE,vjust = -2)+
  geom_vline(xintercept = times["t1"]/24,colour="#00b159",linetype="dashed")+
  geom_vline(xintercept = times["t2"]/24,colour="#00b159",linetype="dashed")+
  geom_hline(yintercept = Past,colour="blue",linetype="dashed")+
  geom_hline(yintercept = eta,colour="blue",linetype="dashed")

plot1
