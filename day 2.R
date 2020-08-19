library(EpiModel)

set.seed(0)


### modeling msm in steady partnerships

## avg msm in 0.4 parternships
## average partnership lasts 90 time steps

net1 <- network_initialize(100)

coef.diss.1 <- dissolution_coefs(~offset(edges), 90)

coef.diss.1

## dissolution coefficient
log(89)

### target statistic 
## can form 20 dyads (100*0.4)/2 

fit1 <- netest(net1,formation =  ~edges, target.stats = 20, coef.diss = coef.diss.1)

summary(fit1)

## diagnostic simulation

sim1 <- netdx(x = fit1,nsteps =  1000, nsims = 10, keep.tedgelist = T)


## diagnostic plots 

plot(sim1, type = "formation")

plot(sim1, type = "duration")

plot(sim1, type = "dissolution")

tel <- as.data.frame(sim1)

head(tel)

hist(tel$duration)

mean(tel$duration[tel$onset < 100])

sum(tel$terminus.censored == TRUE)


plot(tel$onset, tel$terminus)

table(c(tel$head,tel$tail))

hist(table(c(tel$head,tel$tail)))


############### scenario 2

## avg msm in 0.4 parternships
## average partnership lasts 90 time steps

net1 <- network_initialize(1000)

coef.diss.1 <- dissolution_coefs(~offset(edges), 90)

coef.diss.1

## dissolution coefficient
log(89)

### target statistic 
## can form 20 dyads (100*0.4)/2 

fit1 <- netest(net1,formation =  ~edges, target.stats = 200, coef.diss = coef.diss.1)

summary(fit1)

## diagnostic simulation

sim1 <- netdx(x = fit1,nsteps =  1000, nsims = 10, keep.tedgelist = T)


## diagnostic plots 

plot(sim1, type = "formation")

plot(sim1, type = "duration")

plot(sim1, type = "dissolution")

tel <- as.data.frame(sim1)

head(tel)

hist(tel$duration)

mean(tel$duration[tel$onset < 100])

sum(tel$terminus.censored == TRUE)


plot(tel$onset, tel$terminus)

table(c(tel$head,tel$tail))

hist(table(c(tel$head,tel$tail)))

############### scenario 3

## msm community 50% black 50% white

# control different aspects of momentary degree distribution

# ego centric data 
# mean degree = .9
# 83% of relationships are racially homogenious 

# homogenious = 100, heterogenous = 200

n <- 500
net3 <- network_initialize(n)
net3 <- set_vertex_attribute(net3, "race", rep(c("B","W"), each = n/2))
net3

## formation model 

# concurrent -- number of nodes of degree 2 or more

form.formula.3 <-  ~edges + nodematch('race') + degree(0) + concurrent #why no degree 1? to avoid issues with colinearity

target.stats.3 <- c(0.9*n/2, #
                    (.9*n/2)*(5/6), #
                    0.36*n, # degree 0
                    0.18*n) # degree 2+

diss.formula.3 <- ~offset(edges) + offset(nodematach("race"))

## fit model for 
fit3 <- netest(net3, formation = form.formula.3, target.stats = target.stats.3, coef.diss = (dissolution_coefs(~offset(edges) + offset(nodematch("race")), c(200,100))))

sim3 <- netdx(fit3, nsteps = 1000, nsims = 10, keep.tedgelist = T)


plot(sim3, "formation")







