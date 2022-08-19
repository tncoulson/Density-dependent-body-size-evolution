setwd("/Users/timcoulson/Applications/OneDrive - Nexus365/Guppy resource accrual evolution paper/figure 3")
# load in parameter values
source('params.r')
# load in functions
source('functions.r')
# define array for storing results
outp <- array(NA,c(20,11)) # scalar stats
struct <- array(NA,c(n.bins,n.bins,20)) # population structure
rv <- array(NA,c(n.bins,20)) # reproductive values
r.mats <- s.mats <- array(NA,c(n.bins,n.bins,n.steps))
surv.rate.by.age <- array(NA,c((n.ages-1),20))

# main loop
for (j in 1:20){ # loop through the 20 life histories
  # specify the growth function
  grow.b <- grow.b.vals[j]
  grow.a <- 8-grow.b*4
  # specify an array of population structure by time for the simulation
  N <- array(NA,c(n.bins,n.bins,n.steps))
  # random starting population structure
  N[,,1] <- diag(1,n.bins,n.bins)
  # loop through the simulation for each life history
  for (i in 2:n.steps){
    # run the IPM
    n <- sum(N[,,i-1])
    s <- diag(surv(z,surv.a,surv.b,surv.c,n))
    z.a <- (grow.a/(1-grow.b)) * 0.8
    r <- repr(z,repr.a,repr.b,repr.c,n)
    r <- ifelse(z<z.a,0,r)
    r <- diag(r)
    g <- outer(z,z,grow,grow.a,grow.b,grow.c)
    g <- g/matrix(as.vector(apply(g,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    d <- outer(z,z,offs,offs.a,offs.b,offs.c)
    d <- d/matrix(as.vector(apply(d,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    m <- d%*%r + g%*%s
    N[,,i] <- m%*%N[,,i-1]
  }
  # calculate the mean survival rate
  sad <- apply(N[,,n.steps],1,sum)
  rvv <- apply(N[,,n.steps],2,sum)
  mean.s <- sum(diag(s)*sad)/sum(sad)
  # calculate the mean reproductive rate
  mean.r <- sum(diag(r)*sad)/sum(sad)
  # population structure
  struct[,,j] <- N[,,i]
  # reproductive value structure
  temp <- get.eigen.stuff(m)
  # reproductive value
  rv <- temp[[3]]
  # extract survival and reproductive matrices
  r.mats[,,j] <- d%*%r
  s.mats[,,j] <- g%*%s
  # proportion of the population that are sexually mature
  prop.sad <- sad/sum(sad)
  prop.alive <- sum(sad[z>=z.a])
  n <- r.mats[,,j]%*%sad
  n <- n/sum(n)
  # survival rates by class
  adults <- ifelse(z>=z.a,prop.sad,0)
  juveni <- ifelse(z<z.a,prop.sad,0)
  adults.surv <- sum(s%*%adults)
  juveni.surv <- sum(s%*%juveni)
  adults.repr <- sum(r%*%adults)
  prop.mature <- 0
  for (k in 1:50){
    n <- s.mats[,,j]%*%n
    prop.mature <- prop.mature + sum(n[z>=z.a])
    n <- ifelse(z>=z.a,0,n)
  }
  print(j)
  outp[j,] <- c(grow.a,grow.b,z.a,sum(sad),mean.s,mean.r,prop.mature,prop.alive,adults.surv,juveni.surv,adults.repr)
}

sad.array <- array(NA,c(n.bins,20))
prop.above <- rep(NA,20)
for (i in 1:20) sad.array[,i] <- apply(struct[,,i],1,sum)
for (i in 1:20) prop.above[i] <- sum(ifelse(z<outp[i,3],0,sad.array[,i]))
prop.above <- prop.above/outp[,4]

plotter<- TRUE
if (plotter==TRUE) {
  source('size to age-stage.r')
  # do the age-structured calculations
  # extract age-specific survival rates and lx schedules
  # first calculate survival rates by age from size-specific functions
  for (i in 1:20){
    temp <- age.mat.quants[[i]][[2]]
    sad <- tapply(temp,age,sum)
    surv.rate.by.age[,i] <- sad[-1]/sad[-n.ages]
  }
  # set up plotting grid
  # calculate survivorship functions
  survivorships <- array(NA,c(n.ages,20))
  for (i in 1:20) survivorships[,i] <- c(1,cumprod(surv.rate.by.age[,i]))
  source('final figure 1.r')
}

save.image("figure 3 code run.RData")
