####### MODEL 1
age.x <- 0:13
surv.y <- rep(1/(1+exp(-1)),13)
lx.m1 <- c(1,cumprod(surv.y))
par(mfrow=c(2,3))
plot(age.x,lx.m1[1:14],type='l',xlab='Age',ylab='Survivorship')
text(11,0.95,'(A)',cex=1.5)
plot(1:13,surv.y,type='l',xlab='Age',ylab='Survival rate')
text(11,0.995,'(B)',cex=1.5)
lep <- 1/(1-surv.y[1])
plot(z,rep(lep,400),type='l',xlab='Body size',ylab='Life expectancy')
points(outp[,3]*1.25,rep(lep,20),col=1:20,pch=15,cex=1.5)
text(29,5.05,'(C)',cex=1.5)
####### MODEL 2
# extract age-specific survival rates and lx schedules
# first calculate survival rates by age from size-specific functions
surv.rate.by.age <- array(NA,c((n.ages-1),20))
for (i in 1:20){
  temp <- age.mat.quants[[i]][[2]]
  sad <- tapply(temp,age,sum)
  surv.rate.by.age[,i] <- sad[-1]/sad[-n.ages]
}
# set up plotting grid
# calculate surviroship functions
survivorships <- array(NA,c(15,20))
for (i in 1:20) survivorships[,i] <- c(1,cumprod(surv.rate.by.age[,i]))
# plot the lx schedules from the survival rates
plot(0:13,survivorships[1:14,1],type='l',ylim=c(0,1),xlab='Age',ylab='Surivorship')
for (i in 2:20) lines(0:13,survivorships[1:14,i],col=i)
# plot survival by age
text(11,0.95,'(D)',cex=1.5)
plot(1:(n.ages-2),surv.rate.by.age[1:(n.ages-2),1],type='l',ylim=c(0.5,1),xlab='Age',ylab='Survival rate')
for (i in 2:20){
  lines(1:(n.ages-2),surv.rate.by.age[1:(n.ages-2),i],col=i)
}
text(2.5,0.975,'(E)',cex=1.5)
# plot life expectancy against size for each life history
le <- rep(NA,20)
I <- diag(1,n.bins)
plot(z,apply(solve(I-s.mats[,,1]),2,sum),type='l',col=i,ylim=c(0,30),xlab='Body size',ylab='Life expectancy')
for (i in 2:20) lines(z,apply(solve(I-s.mats[,,i]),2,sum),col=i)
for (i in 1:20){
  xx <- match(1,ifelse(z<outp[i,3],0,ifelse(z>(outp[i,3]+0.0875),0,1)))
  yy <- apply(solve(I-s.mats[,,i]),2,sum)[xx]
  le[i] <- yy
  points(z[xx],yy,col=i,pch=16)
  xx <- match(1,ifelse(z<(outp[i,3]*1.25),0,ifelse(z>((outp[i,3]+0.0875)*1.25),0,1)))
  yy <- apply(solve(I-s.mats[,,i]),2,sum)[xx]
  points(z[xx],yy,col=i,pch=15,cex=1.5)
}
text(4,29,'(F)',cex=1.5)
