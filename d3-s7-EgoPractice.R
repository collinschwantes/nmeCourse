
##
## PRACTICE WITH EGOCENTRIC DATA AND TARGET STATISTICS
## Day 3 | Network Modeling for Epidemics
##

### Create network

library(EpiModel)
mynet <- network_initialize(2000)

sex <- c(rep(1, 1000), rep(2, 1000))
mynet <- set_vertex_attribute(mynet, 'group', sex)

cmty <- c(rep(1,400), rep(2,600), rep(1,400), rep(2,600))
mynet <- set_vertex_attribute(mynet, 'cmty', cmty)

### Fit and simulate model

formation <- ~edges+nodematch('cmty')+degrange(from=3)+
  degree(1, by='group')

target.stats <- c(850, 600, 425, 550, 350)

myfit <- netest(mynet,
                formation=formation,
                target.stats = target.stats,
                coef.diss = dissolution_coefs(~offset(edges), 60))

mydx <- netdx(myfit, nsims=10, nsteps=100)
mydx
get_nwstats(mydx)
plot(mydx)

## Conduct disease sim

mycontrol <- control.net("SIS", ncores = 4, nsteps = 1000, nsims = 1,	nwstats.formula = ~edges+nodematch('cmty') + 
                           degree(0:5, by='group'), verbose = TRUE)
myinit <- init.net(i.num = 10, i.num.g2 = 10)
myparam <-param.net(inf.prob = 0.8, inf.prob.g2 = 0.2,
                    act.rate = 0.6,
                    rec.rate = 0.05, rec.rate.g2 = 0.05)

mySIS <- netsim(myfit, param = myparam, control = mycontrol,
                init = myinit)

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

par(resetPar())

plot(mySIS)

sisDF <- as.data.frame(mySIS, out = "mean")

head(sisDF)

## crude cumulative incidence 
sum(sisDF$si.flow,na.rm = NA)
sum(sisDF$si.flow.g2,na.rm = NA)

par(mfrow = c(1,2))

## prev 
plot(mySIS, y = c("i.num", "i.num.g2"), popfrac = TRUE,
     qnts = FALSE, ylim = c(0, 0.4), legend = F)

## inc
plot(mySIS, y = c("si.flow", "si.flow.g2"), 
     qnts = FALSE, ylim = c(0, 10), legend = F)

## Conduct diagnostics

get_nwstats(mySIS)
plot(mySIS, type = "formation", sim.lines = TRUE)

### netest = approximation false in sim
## 

