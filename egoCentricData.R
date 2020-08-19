library(EpiModel)

## create our empty network of nodes

mynet <- network_initialize(2000)

## add group attributes

# male -female
# location

sex <- c(rep(1, 1000), rep(2, 1000))
mynet <- set_vertex_attribute(mynet, 'group', sex)

cmty <- c(rep(1,400), rep(2,600), rep(1,400), rep(2,600))
mynet <- set_vertex_attribute(mynet, 'cmty', cmty)

## count edges in sample

## 8 edges with f and 9 edges with m -- taking average 8.5

numEdges <- 850

## specify mean degree

meanF <- mean(c(1,0,2,1,1,1,1,0,0,1))
meanM <- mean(c(2,0,0,1,1,0,2,1,2,0))

## limit degree range

## concurrency

### formation functions
## dissolution -- duration
# nodematch - no homophily witin sex, but homophily within location

formation <- ~edges+nodematch('cmty')+degrange(from=3)+
  degree(1, by='group')

# add target stats
target.stats <- c(850, 600, 425, 550, 350)

# model fit
# model diagnostics

# simulate model