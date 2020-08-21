##
## Tutorial: Feedback Through Nodal Attribute Changes
## Day 4 | Network Modeling for Epidemics
##

library(EpiModel)


# Network Model -----------------------------------------------------------

# Initialize network
n <- 500
nw <- network_initialize(n)

# Set infection status attribute
prev <- 0.2
infIds <- sample(1:n, n*prev)
nw <- set_vertex_attribute(nw, "status", "s")
nw <- set_vertex_attribute(nw, "status", "i", infIds)
get_vertex_attribute(nw, "status")

# Mean degree by infection status
mean.deg.inf <- 0.3
inedges.inf <- mean.deg.inf * n * prev
inedges.inf

mean.deg.sus <- 0.8
inedges.sus <- mean.deg.sus * n * (1 - prev)
inedges.sus

# Calculate number of edges as sum of the nodefactor statistics divided by 2
edges <- (inedges.inf + inedges.sus)/2
edges

# Exact solution for nodematch term under random mixing
p <- inedges.sus/(edges*2)
q <- 1 - p
nn <- p^2
np <- 2*p*q
pp <- q^2
round(nn + pp, 3)

# Simulation based solution for nodematch term under random mixing
fit <- netest(nw,
              formation = ~edges + nodefactor("status"),
              target.stats = c(edges, inedges.sus),
              coef.diss = dissolution_coefs(~offset(edges), duration = 1))
sim <- netdx(fit, dynamic = FALSE, nsims = 1e4,
             nwstats.formula = ~edges + nodematch("status"))
stats <- get_nwstats(sim)
head(stats)
round(mean(stats$nodematch.status/stats$edges), 3)

# Nodematch term with slightly more within-group ties
nmatch <- edges * 0.91

# Model 1 estimation
formation <- ~edges + nodefactor("status") + nodematch("status")
target.stats <- c(edges, inedges.sus, nmatch)

coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 50)

est <- netest(nw, formation, target.stats, coef.diss)

# Network diagnostics
dx <- netdx(est, nsims = 10, nsteps = 500, ncores = 5,
            nwstats.formula = ~edges + 
                               meandeg + 
                               nodefactor("status", levels = NULL) + 
                               nodematch("status"), verbose = FALSE)
dx
plot(dx)
plot(dx, type = "dissolution")

# Model 2 estimation
est2 <- netest(nw, formation = ~edges, target.stats = edges, coef.diss)

# Network diagnostics
dx2 <- netdx(est2, nsims = 10, nsteps = 1000, ncores = 5,
             nwstats.formula = ~edges + 
                                meandeg + 
                                nodefactor("status", levels = NULL) + 
                                nodematch("status"), verbose = FALSE)
dx2
plot(dx2)


# Epidemic Model ----------------------------------------------------------

# Parameterize model
param <- param.net(inf.prob = 0.03)

# Initial conditions
init <- init.net()

# Control settings
control <- control.net(type = "SI", nsteps = 500, nsims = 5, ncores = 5,
                       resimulate.network = TRUE, tergmLite = FALSE,
                       nwstats.formula = ~edges + 
                                          meandeg + 
                                          nodefactor("status", levels = NULL) + 
                                          nodematch("status"))

# Model 1 simulation
sim <- netsim(est, param, init, control)

# Model 2 simulation
sim2 <- netsim(est2, param, init, control)

# Epidemic Results
par(mfrow = c(1,2))
plot(sim, main = "Seroadaptive Behavior")
plot(sim2, main = "No Seroadaptive Behavior")


par(mfrow = c(1,1))
plot(sim, y = "i.num", popfrac = TRUE, sim.lines = FALSE, qnts = 1)
plot(sim2, y = "i.num", popfrac = TRUE, sim.lines = FALSE, qnts = 1, 
     mean.col = 2, qnts.col = 2, add = TRUE)
legend("topleft", c("Serosort", "Non-Serosort"), lty = 1, lwd = 3,
       col = c(4, 2), cex = 0.9, bty = "n")

# Network results
plot(sim, type = "formation")

plot(sim, type = "formation", 
     stats = c("nodefactor.status.s", 
               "nodefactor.status.i"))

plot(sim, type = "formation", 
     stats = c("edges", "nodematch.status"), 
     sims = 1, sim.lwd = 2)

