##
## Tutorial: Epidemic Models with Demography
## Day 4 | Network Modeling for Epidemics
##

library(EpiModel)
library(ndtv)


# Network Model -----------------------------------------------------------

# Initialize network
nw <- network_initialize(n = 1000)
#nw <- set_vertex_attribute(nw, attrname = "risk", value = rep(0:1, each = 250))

# Model 1: Random Mixing

# Formation model
formation <- ~edges + concurrent 

# Input the appropriate target statistics for each term
target.stats <- c(250, 100)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 40, d.rate = 0.0015)
coef.diss

# Dissolution model with error
## coef.diss.high <- dissolution_coefs(dissolution = ~offset(edges),
##                                duration = 40, d.rate = 0.05)

# Fit model 1
est1 <- netest(nw, formation, target.stats, coef.diss)

# Model 1 diagnostics
dx1 <- netdx(est1, nsims = 10, nsteps = 1000, ncores = 4,
             nwstats.formula = ~edges + concurrent)
dx1
plot(dx1)

par(mfrow= c(1,2))

# Epidemic model: SI model with highly effective intervention
param <- param.net(inf.prob = 0.5, inter.eff = 0.50, inter.start = 15)
init <- init.net(i.num = 5)
control <- control.net(type = "SI", nsteps = 100, nsims = 10, ncores = 5)
sim <- netsim(est1, param, init, control)
plot(sim)


### counter factual - adjust implementation time

param <- param.net(inf.prob = 0.5, inter.eff = 0.50, inter.start = 55)
init <- init.net(i.num = 5)
control <- control.net(type = "SI", nsteps = 100, nsims = 10, ncores = 5)
sim2 <- netsim(est1, param, init, control)
plot(sim2)



### create animation

nw <- get_network(sim)
nw <- color_tea(nw, verbose = FALSE)

par(mfrow = c(1,1))
plot(nw, col = nw$val$status)

slice.par <- list(start = 1, end = 10, interval = 1, 
                  aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))

compute.animation(nw, slice.par = slice.par, verbose = TRUE)

render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))

# question 1:





# Model 2: Assortative Mixing and Differential Degree

# Formation model
formation <- ~edges + nodefactor("risk") 

# Input the appropriate target statistics for each term
target.stats <- c(125, 187.5)

# Fit model 2
est2 <- netest(nw, formation, target.stats, coef.diss)

# Review coefficients for nodefactor and nodematch terms
summary(est2)

# Model 2 diagnostics
dx2 <- netdx(est2, nsims = 10, nsteps = 1000, ncores = 4,
             nwstats.formula = ~edges +
                                nodefactor("risk", levels = NULL)
                                )
dx2
plot(dx2)


# Epidemic Simulation -----------------------------------------------------

# Parameterize model
param <- param.net(inf.prob = 0.1, act.rate = 5,
                   a.rate = 0.001, ds.rate = 0.001, di.rate = 0.002)

# Initial conditions
init <- init.net(i.num = 50)

# Control settings with additional parameters
# Without tergmLite
# control.full <- control.net(type = "SI", nsteps = 300, nsims = 1, ncores = 1,
#                             resimulate.network = TRUE, epi.by = "risk",
#                             tergmLite = FALSE)
# With tergmLite
control.tl <- control.net(type = "SI", nsteps = 300, nsims = 1, ncores = 1,
                          resimulate.network = TRUE, epi.by = "risk",
                          tergmLite = TRUE)

# Simulation with and without tergmLite
sim1.full <- netsim(est1, param, init, control.full)
sim1.tl <- netsim(est1, param, init, control.tl)

sim1.full
sim1.tl

# Main Epidemic Simulation

# Control settings
control <- control.net(type = "SI", nsteps = 300, nsims = 10, ncores = 8,
                       resimulate.network = TRUE, epi.by = "risk", tergmLite = T)

# Simulations for Model 1 and Model 2
sim1 <- netsim(est1, param, init, control)
sim2 <- netsim(est2, param, init, control)

# Print Model 1 simulation
sim1

# Plot model diagnostics
par(mfrow = c(1, 2))
plot(sim1, type = "formation", stats = "edges", ylim = c(0, 250), qnts = 1)
plot(sim2, type = "formation", stats = "edges", ylim = c(0, 250), qnts = 1)

# Plot epidemiologic outcomes

# Overall prevalence
par(mfrow = c(1, 1))
plot(sim1, y = "i.num", qnts = 1, main = "Total Prevalence", mean.col = 3, qnts.col = 3)
plot(sim2, y = "i.num", qnts = 1, mean.col = 6, qnts.col = 6, add = TRUE)
legend("topleft", c("Model 1", "Model 2"), lwd = 3, col = c(3, 6), bty = "n")

# Prevalence by risk group
par(mfrow = c(1, 2))
plot(sim1, y = c("i.num.risk0", "i.num.risk1"),  legend = TRUE, qnts = 1,
     ylim = c(0, 500), main = "M1: Prevalence by Group")
plot(sim2, y = c("i.num.risk0", "i.num.risk1"), legend = TRUE,  qnts = 1,
     ylim = c(0, 500), main = "M2: Prevalence by Group")

