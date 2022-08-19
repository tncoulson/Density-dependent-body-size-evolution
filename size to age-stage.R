# store age.structured matrices at equilibrium
#age.mat.res <- array(NA,c(n.bins*n.ages,n.bins*n.ages,20))
age.mat.quants <- NULL
age.spec.wt <- array(NA,c(20,n.ages))
props <- rep(NA,20)
# converts size structured model into age-stage model
for (j in 1:20){ # loop through the 20 life histories
  age.mat <- matrix(0,nrow=n.bins*n.ages,ncol=n.bins*n.ages) # 15 age classses
  for (i in 1:n.ages) age.mat[1:n.bins,((i-1)*n.bins+1):(i*n.bins)] <- r.mats[,,j]
  for (i in 2:n.ages) age.mat[((i-1)*n.bins+1):(i*n.bins),((i-2)*n.bins+1):((i-1)*n.bins)] <- s.mats[,,j]
  age.mat[((n.ages-1)*n.bins+1):(n.ages*n.bins),((n.ages-1)*n.bins+1):(n.ages*n.bins)] <- s.mats[,,j]
#  age.mat.res[,,j] <- age.mat
  age.mat.quants[[j]] <- get.eigen.stuff(age.mat)
}
# calculate key quantities (rv, stable stage distribution, pop growth)
# for (i in 1:20) age.mat.quants[[i]] <- get.eigen.stuff(age.mat.res[,,i])

# extract 
age <- rep(1:n.ages,each=n.bins)
z.l <- rep(z,times=n.ages)
age.struct.array <- matrix(NA,nrow=n.ages,ncol=20)
for (i in 1:20){
  sad <- age.mat.quants[[i]][[2]]
  age.struct <- tapply(sad,age,sum)
  age.struct.array[,i] <- age.struct*outp[i,4]
  num <- z.l*sad
  age.spec.wt[i,] <- tapply(num,age,sum)/tapply(sad,age,sum)
  props[i] <- sum(sad[z.l>outp[i,3]])
}

