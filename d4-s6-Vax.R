##
## Tutorial: A Simple Vaccine Intervention
## Day 4 | Network Modeling for Epidemics
##

library(EpiModel)

# Estimate simple temporal ERGM
nw <- network_initialize(n = 100, directed = FALSE)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Epidemic model: SI model with highly effective intervention
param <- param.net(inf.prob = 0.5, inter.eff = 0.96, inter.start = 25)
init <- init.net(i.num = 5)
control <- control.net(type = "SI", nsteps = 100, nsims = 10, ncores = 5)
sim <- netsim(est, param, init, control)
plot(sim)

# Epidemic model: SIS model with slightly less effective intervention
param <- param.net(inf.prob = 0.5, inter.eff = 0.8, inter.start = 100, 
                   rec.rate = 0.07)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 250, nsims = 10, ncores = 5)
sim <- netsim(est, param, init, control)
plot(sim)

