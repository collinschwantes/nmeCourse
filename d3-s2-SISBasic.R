##
## Tutorial: SIS Epidemic Across a Dynamic Network
## Day 3 | Network Modeling for Epidemics
##

library(EpiModel)


# Network model estimation ------------------------------------------------

# Initialize network
nw <- network_initialize(n = 500)

# Define the formation model
formation <- ~edges + concurrent + degrange(from = 4)

# Input the appropriate target statistics for each term
target.stats <- c(175, 110, 0)

# Parameterize the dissolulation model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

# Review the arguments for the netest function
args(netest)

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 10, nsteps = 1000,
            nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent,
            keep.tedgelist = TRUE)
dx

# Use multiple cores
parallel::detectCores()

dx <- netdx(est, nsims = 10, nsteps = 1000, ncores = 4,
            nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent,
            keep.tedgelist = TRUE)
dx

# Plot the formation diagnostics
plot(dx)

nwstats1 <- get_nwstats(dx, sim = 1)
head(nwstats1, 20)

# Plot the dissolution diagnostics
par(mfrow = c(1, 2))
plot(dx, type = "duration")
plot(dx, type = "dissolution")

tel <- as.data.frame(dx, sim = 1)
head(tel, 20)

# Static (cross-sectional ERGM) model diagnostics (if the dynamic dx are not 
# looking good)
dx.static <- netdx(est, nsims = 10000, dynamic = FALSE)
dx.static

par(mfrow = c(1,1))
plot(dx.static, sim.lines = TRUE, sim.lwd = 0.1)

nwstats2 <- get_nwstats(dx.static)
head(nwstats2, 20)

#### no concurrency 


# Network model estimation modified network ----------------------------------------

# Initialize network
nw <- network_initialize(n = 500)

# Define the formation model
formation <- ~edges + concurrent + degrange(from = 4)

# Input the appropriate target statistics for each term
target.stats <- c(175,0,0)

# Parameterize the dissolulation model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

# Review the arguments for the netest function
args(netest)

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
# dx <- netdx(est, nsims = 10, nsteps = 1000,
#             nwstats.formula = ~edges + meandeg + degree(0:4),
#             keep.tedgelist = TRUE)
# dx

# Use multiple cores
parallel::detectCores()

dx <- netdx(est, nsims = 10, nsteps = 1000, ncores = 4,
            nwstats.formula = ~edges + meandeg + degree(0:4),
            keep.tedgelist = TRUE)
dx

# Plot the formation diagnostics
plot(dx)

nwstats1 <- get_nwstats(dx, sim = 1)
head(nwstats1, 20)

# Plot the dissolution diagnostics
par(mfrow = c(1, 2))
plot(dx, type = "duration")
plot(dx, type = "dissolution")

tel <- as.data.frame(dx, sim = 1)
head(tel, 20)

# Static (cross-sectional ERGM) model diagnostics (if the dynamic dx are not 
# looking good)
dx.static <- netdx(est, nsims = 10000, dynamic = FALSE)
dx.static

par(mfrow = c(1,1))
plot(dx.static, sim.lines = TRUE, sim.lwd = 0.1)

nwstats2 <- get_nwstats(dx.static)
head(nwstats2, 20)



# Epidemic model simulation -----------------------------------------------

# Parameterizing an SIS epidemic 
param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.1)

# Initial conditions
init <- init.net(i.num = 10)

# Control settings
#control <- control.net(type = "SIS", nsims = 5, nsteps = 500) 
control <- control.net(type = "SIS", nsims = 5, nsteps = 500, ncores = 5)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)

# Print the output to show the model contents
sim

# Plot the output to epi stats
par(mfrow = c(1, 1))
plot(sim)

summary(sim,at = 113)

# Many different arguments for plot.netsim
par(mfrow = c(1, 2))
plot(sim, sim.lines = TRUE, mean.line = FALSE, qnts = FALSE, popfrac = TRUE)
plot(sim, mean.smooth = FALSE, qnts = 1, qnts.smooth = FALSE, popfrac = TRUE)

# Us the y argument to pull out non-default stats, such as incidence
par(mfrow = c(1,1))
plot(sim, y = c("si.flow", "is.flow"), qnts = FALSE, 
     ylim = c(0, 25), legend = TRUE, main = "Flow Sizes")

# Static network plot from one sim at two time points
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
plot(sim, type = "network", col.status = TRUE, at = 50, sims = 1)
plot(sim, type = "network", col.status = TRUE, at = 100, sims = 1)

# Summary stats
summary(sim, at = 500)

# Convert model to a data frame for further analysis
df <- as.data.frame(sim)
head(df, 10)
tail(df, 10)

# Extracting mean values also possible
df <- as.data.frame(sim, out = "mean")
head(df, 10)
tail(df, 10)

# Extract the fully dynamic network for further analysis
nw1 <- get_network(sim, sim = 1)
nw1

# Temporal edgelist
nwdf <- as.data.frame(nw1)
head(nwdf, 25)

# A transmission matrix contains the time-ordered chain of tranmission
tm1 <- get_transmat(sim, sim = 1)
head(tm1, 10)

# Plotting with ggplot
df <- as.data.frame(sim)
df.mean <- as.data.frame(sim, out = "mean")

library(ggplot2)
ggplot() +
  geom_line(data = df, mapping = aes(time, i.num, group = sim), alpha = 0.25,
            lwd = 0.25, color = "firebrick") +
  geom_bands(data = df, mapping = aes(time, i.num),
             lower = 0.1, upper = 0.9, fill = "firebrick") +
  geom_line(data = df.mean, mapping = aes(time, i.num)) +
  theme_minimal()

