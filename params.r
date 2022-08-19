# size of approximated matrix
n.bins <- 400
# number of steps in simulation
n.steps <- 200
# number of age classes
n.ages <- 61
# range of mid-point values
z <- seq(0.1,35,length.out=n.bins)
# survival parameters -- intercept, size slope, density slope
surv.a <- -0.875
surv.b <- 0.15
surv.c <- 0
# reproductive parameters -- intercept, size slope, density slope
repr.a <- 1
repr.b <- 0
repr.c <- -0.001
# inheritance parameters -- intercept, size slope, variance
offs.a <- 4
offs.b <- 0
offs.c <- 1
# growth parameters (most defined elsewhere -- growth variance)
grow.c <- 1
# range of life histories -- defined by growth parameters
grow.b.vals <- seq(0.3,0.8,length.out=20)
  
# for table 1
surv.a2 <- 0.25
surv.b2 <- 0.125
surv.c2 <- -0.001
repr.a2 <- -1
repr.b2 <- 0
repr.c2 <- 0