library("Sim.DiffProc") 


# Simulation parameters
bigT = 20.; 
Dt = 0.01;

# model parameters
C = 250.; # membrane capacitance
EL = 0.;# -70; # resting membrane potential
gL = 1./60.; # membrane leak conductance
EE = 70.;# 0.; # excitation synaptic reversal potential
EI = -5.;# -75.; # inhibitory synaptic reversal potential
tauE = 0.2; # excitatory conducatnce time constant
tauI = 2.; # inhibitory conducatnce time constant
cE = 7.1; # excitatory peak conductance
cI = 3.7; # inhibitory peak conductance

# RE = 300; RI = 160; # balance excitation inhibition
# RE = 150; RI = 80; # balance excitation inhibition
# rates = data.frame(gEq = c(300. , 5.) , gIq = c(160. , 0.65))
# rates.mod = lm(gIq ~ gEq , data=rates)
REs = 10^seq(log10(1),log10(200),length.out=20)
# REs = seq(1,5,length.out=20)
vars <- c()
means <- c()
ress <- list()

for (ii in 1:length(REs)){
  RE = REs[ii];
  RI = 0.533*RE;
  gE0 = cE*tauE*RE; gI0 = cI*tauI*RI; # Richardson & Gerstner 2005
  sigE = cE*sqrt(tauE*RE/2.); sigI = cI*sqrt(tauI*RI/2.);
  
  print(RE)
  print(sigE/gE0) # << 1 for diffusion approximation
  print(sigI/gI0) #  to hold
  
  fx <- expression((1./C)*(-gL*(x - EL) - y*(x - EE) - z*(x - EI)) , 
                   -(1./tauE)*(y-gE0) , -(1./tauI)*(z-gI0) )
  gx <- expression(0 , sqrt(2)*sigE , sqrt(2)*sigI)
  
  res <- snssde3d(x0=c(15. , gE0 , gI0),drift=fx,diffusion=gx,
                  M=250,t0=0,T=bigT,Dt=Dt,method="euler",N=bigT/Dt)
  
  summary(res,at=bigT)
  vars <- append(vars , var(t(tail(na.omit(res$X),n=1))))
  means <- append(means , mean(t(tail(na.omit(res$X),n=1))))
  save(res, file=sprintf("dat/resRE%f.rda",RE))
}
print(means)
print(vars)

plot(log10(REs),means)
plot(log10(REs),vars)
