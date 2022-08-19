# where to plot figure labels
place.letter <- function(wh,a,c){
  temp <- par("usr")
  x.range <- temp[2]-temp[1]
  y.range <- temp[4]-temp[3]
  if (wh=="tl") text(temp[1]+x.range*0.05,temp[4]-y.range*0.1,a,cex=c,adj=c(0,0))
  if (wh=='bl') text(temp[1]+x.range*0.05,temp[3]+y.range*0.1,a,cex=c,adj=c(0,1))
  if (wh=='tr') text(temp[2]-x.range*0.05,temp[4]-y.range*0.1,a,cex=c,adj=c(0,0))
  if (wh=='br') text(temp[2]-x.range*0.005,temp[3]+y.range*0.1,a,cex=c,adj=c(0,1))
}

# survival function
surv <- function(z,a,b,c,n) 1/(1+exp(-(a+b*z+c*n)))

# reproduction function
repr <- function(z,a,b,c,n){
  repr <- exp(a+b*z+c*n)
}

# inheritance function
offs <- function(z1,z,a,b,c) dnorm(z1,mean=a+b*z,sd=c)

# growth function
grow <- function(z1,z,a,b,c) dnorm(z1,mean=a+b*z,sd=c)

# calculate dominant eigenvalue and vectors
get.eigen.stuff <- function(mat){ 
  sz <- dim(mat)[1]
  t.now <- runif(sz)
  t.now <- t.now/sum(t.now)
  t.next <- mat%*%t.now
  t.next <- t.next/sum(t.next)
  i <- 0
  while (sum(abs(t.next-t.now))>0.0000001){
    i <- i+1
    t.now <- t.next
    t.next <- mat%*%t.now
    lambda <- sum(t.next)/sum(t.now)
    t.next <- t.next/sum(t.next)
  }
  r.now <- runif(sz)
  r.now <- r.now/sum(r.now)
  r.next <- r.now%*%mat
  r.next <- r.next/sum(r.next)
  while (sum(abs(r.next-r.now))>0.0000001){
    r.now <- r.next
    r.next <- r.now%*%mat
    r.next <- r.next/sum(r.next)
  }
  r.next <- as.vector(r.next)
  t.next <- as.vector(t.next)
  return(list(lambda,t.next,r.next))
}

covar <- function(x,y,n){
  sum(x*y*n)/sum(n)-sum(x*n)/sum(n)*sum(y*n)/sum(n)
}