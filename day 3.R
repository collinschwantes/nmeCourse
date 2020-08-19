## day 3 

library(EpiModel)

nw <- network_initialize(n = 500)

#define formation model
# edges almost always there, like intercept in a linear model
# edges is a function of mean degree of persons in pop time pop size /2
# concurrent is a count of nodes that have 2 or more edges at any one time
# degrange parameterizes ranges of degree distribution. From = starting point,
#  to = end, default is infinity

formation <- ~edges + concurrent + degrange(from = 4)

# edges, concurrent, degrange
# in this case we are constraining degrange to less than 4 degrees
target.stats <- c(175, 110, 0)

## edge persistance formula
## offset - analytically calculated, not estimated in model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)

## edapprox - approximation method for est. when mean degree is low and duration is long
# 
est <- netest(nw,formation, target.stats, coef.diss)
## dyadic dependent model - nodal attributes ij but also relationship between jk


## can have a difference between ergm formula and dx formula
# degree(n:m) monitor nodes with these degrees, same as degree(n) + ... + degree(m)
# 
dx <- netdx(est, nsims = 10, nsteps = 1000, nwstats.formula = ~edges + meandeg + degree(0:4)+ concurrent, keep.tedgelist = T)

# use multicores

parallel::detectCores()

dx <- netdx(est, nsims = 10, nsteps = 1000, nwstats.formula = ~edges + meandeg + degree(0:4)+ concurrent, keep.tedgelist = T, ncores = 4)

plot(dx)

tel <- as.data.frame(dx, sim = 1)

head(tel,20)

## could look at cross sectional diagnostics to check that MCMC fitting process is working as expected
# essentially turning ergm to tergm 

## see code in tutorial 
## model mispecification can result in MCMC not converging 





