quartz()

# code to draw the life histories
par(mfrow=c(3,3),mar=c(4,4,1,1),oma=c(1,1,1,1))

# define colours for each life history
rbPal <- colorRampPalette(c('blue','red'))
Col <- rbPal(10)[as.numeric(cut(outp[,3],breaks=10))]

# plot the size-survival function 
plot(z,surv(z,surv.a,surv.b,surv.c,outp[1,4]),type='l',xlab='Body size',ylab='Survival rate',ylim=c(0,1),col=Col[1])
for (i in 2:20) lines(z,surv(z,surv.a,surv.b,surv.c,outp[i,4]),col=Col[i])
place.letter('tl','(A)',1.5)
for (i in 1:20) points(outp[i,3],surv(outp[i,3],surv.a,surv.b,surv.c,outp[i,4]),pch=16,col=Col[i],cex=2)
text(20,0.06,"Survival function")

# plot the inheritance function
contour(z,z,t(d),xlab='Adult body size t',ylab='Offspring body size t+1')
points(outp[,3],rep(4,20),pch=16,col=Col,cex=2)
place.letter('tl','(B)',1.5)
text(20,10,'Inheritance function')

# plot the reproduction function
for (i in 1:20){
  x <- ifelse(z<outp[i,3],0,exp(repr.a + repr.b*z +repr.c*outp[i,4]))
  if (i==1) plot(z,x,type='l',ylim=c(0,1.75),xlab='Body size',ylab='Reproductive rate',xlim=c(5,25)) else lines(z,x,col=Col[i],cex=2)
}
place.letter('tl','(C)',1.5)
text(18,1.65,'Reproduction function')

# plot the growth function -- size t versus size t+1
outpz <- array(NA,c(20,3))
for (i in 1:20){
  grow.b <- grow.b.vals[i]
  grow.a <- 8-grow.b*4
  z.a <- (grow.a/(1-grow.b)) * 0.8
  outpz[i,] <- c(grow.a,grow.b,z.a)
}

plot(z,z,type='l',lty=2,xlab='Size t',ylab='Size t+1')
for(i in 1:20) abline(outpz[i,1],outpz[i,2],col=Col[i])
points(outp[,3],outp[,1]+outpz[,2]*outp[,3],pch=16,col=Col,cex=2)
place.letter('tl','(D)',1.5)
text(20,2,"Growth function")

# plot the von Bertalanffy versions of the growth functions
size.at.age <- array(NA,c(20,20))
size.at.age[1,] <- 4
for (i in 1:20){
  for (j in 2:20){
    size.at.age[j,i] <- outp[i,1]+outp[i,2]*size.at.age[j-1,i]
  }
}

plot(1:20,size.at.age[,20],col=Col[20],type='l',xlab='Age',ylab='Size')
for (i in 19:1){ 
  lines(1:20,size.at.age[,i],col=Col[i])
}
place.letter('tl','(E)',1.5)

asymp.sz <- outp[,3]/0.8
K <- -log(outp[,2])
t0 <- 3+((1/K)*log((asymp.sz-size.at.age[3,])/asymp.sz))
asm <- t0 - log(0.2)/K
points(asm,outp[,3],pch=16,col=Col)

# surv rate by age
surv.rate.by.age <- array(NA,c((n.ages-1),20))
for (i in 1:20){
  temp <- age.mat.quants[[i]][[2]]
  sad <- tapply(temp,age,sum)
  surv.rate.by.age[,i] <- sad[-1]/sad[-n.ages]
}

plot(1:(n.ages-2),surv.rate.by.age[1:(n.ages-2),1],type='l',ylim=c(0,1),xlab='Age',ylab='Survival rate',col=Col[1],xlim=c(0,30))
for (i in 2:20){
  lines(1:(n.ages-2),surv.rate.by.age[1:(n.ages-2),i],col=Col[i])
}
place.letter('tl','(F)',1.5)  

lx <- apply(surv.rate.by.age,2,cumprod)
lx <- rbind(rep(1,20),lx)
lx <- lx[-21,]
plot(0:19,lx[,1],xlab='Age',ylab='`Survivorship',type='l',ylim=c(0,1.1))
lines(0:19,lx[,20],col='red')

x1 <- c(0,asm[1],asm[1],19)
temp <- exp(repr.a + repr.b*z[1] +repr.c*outp[1,4])
y1 <- c(0,0,temp,temp)
lines(x1,y1)
mx <- ifelse((1:20)<(asm[1]+1),0,temp)
lxmx <- lx[,1]*mx

x <- c(0,asm[,1],asm[,1])

polygon(x=c(1:20,20:1),y=c(lxmx,rep(0,20)),col='gray50')
lines(0:19,lx[,1])

x2 <- c(0,asm[20],asm[20],19)
temp <- exp(repr.a + repr.b*z[1] +repr.c*outp[20,4])
y2 <- c(0,0,temp,temp)
lines(x2,y2,col='red')




# examine the survival rates by age
I <- diag(1,n.bins)
leza <- le <- rep(NA,20)
for (i in 1:20){
  yy <- apply(solve(I-s.mats[,,i]),2,sum)
  diff <- z[2]-z[1]
  t1 <- z<(outp[i,3]-diff)
  t2 <- z>outp[i,3]
  print(match(0,(t1+t2))-1)
  leza[i] <- yy[match(0,(t1+t2))]
  if (i==1) plot(z[z<(outp[i,3]/0.8)],yy[z<(outp[i,3]/0.8)],xlab='Body size',ylab='Life expectancy',ylim=c(0,20),xlim=c(0,25),col=Col[i],type='l') else lines(z[z<(outp[i,3]/0.8)],yy[z<(outp[i,3]/0.8)],col=Col[i])
}
place.letter('tl','(G)',1.5)
Z.m <- outp[,3]
points(Z.m,leza,col=Col,pch=16,cex=2)


# draws the age-by-size structure distributions for fastest and slowest life histories
lh <- 1
plot(z,age.mat.quants[[lh]][[2]][age==1],type='l',ylim=c(0,0.2),col='white',xlab='Body size',ylab='',axes=FALSE,xlim=c(0,30))
polygon(x=c(z,rev(z)),y=c(age.mat.quants[[lh]][[2]][age==1]*2.5,rep(0,n.bins)),col=rgb(0,0,1,1),border=FALSE)
for (i in 2:n.ages){
  polygon(x=c(z,rev(z)),y=c(age.mat.quants[[lh]][[2]][age==i]*2.5+(i-1)/100,rep((i-1)/100,n.bins)),col=rgb(0,0,1,1),border=FALSE)
}
abline(v=outp[lh,3],lty=2)

lh <- 20
polygon(x=c(z,rev(z)),y=c(age.mat.quants[[lh]][[2]][age==1]*2.5,rep(0,n.bins)),col=rgb(1,0,0,1),border=FALSE)
for (i in 2:n.ages){
  polygon(x=c(z,rev(z)),y=c(age.mat.quants[[lh]][[2]][age==i]*2.5+(i-1)/100,rep((i-1)/100,n.bins)),col=rgb(1,0,0,1),border=FALSE)
}
abline(v=outp[lh,3],lty=2)
box()
axis(1,labels=TRUE)
mtext(side=2,line=1,'Density by age',cex=0.7)
place.letter('tl','(H)',1.5)

n.array <- t(array(outp[,4],c(20,20)))
surv.rate.by.age <- surv(size.at.age,surv.a,surv.b,surv.c,n.array)
repr.rate.by.age <- repr(size.at.age,repr.a,repr.b,repr.c,n.array)
age.array <- array(1:20,c(20,20))
asm.array <- t(array(asm,c(20,20)))
repr.rate.by.age <- ifelse(age.array<asm.array,0,repr.rate.by.age)
fitness.by.age <- surv.rate.by.age+repr.rate.by.age
age.struct.array[20,] <- age.struct.array[20,]+age.struct.array[21,]
k.fitness.by.age <- age.struct.array[1:20,] * fitness.by.age
adder.array <- ifelse((age.array-1)<asm.array,1,2)
ffs <- array(NA,c(20,2))
for (i in 1:20) ffs[i,] <- tapply(k.fitness.by.age[,i],adder.array[,i],sum)
plot(1:20,ffs[,1],ylim=c(250,640),col='blue')
points(1:20,ffs[,2],col='red')

quartz()
newcol <- c(rep('blue',13),'black',rep('red',6))
lwd.array <- rep(1,20)
lwd.array[14] <- 2
plot(1:20,k.fitness.by.age[,1],type='l',col=newcol[1],lwd=lwd.array[1])
for(i in 2:20) lines(1:20,k.fitness.by.age[,i],col=newcol[i],lwd=lwd.array[i])
lines(1:20,k.fitness.by.age[,14],col=newcol[14],lwd=lwd.array[14])



# proportion of population that are adults
P <- outp[,8]/outp[,4]
# number of juveniles
juv.n <- (outp[,4]-outp[,8])
# number of adults
adu.n <- outp[,8]
# per capita juvenile survival rate
juv.s <- outp[,10]/(1-P)
# per capita adult survival rate
adu.s <- outp[,9]/P
# per capita adult reproductive rate
adu.r <- outp[,11]/P
# per capita adut fitness
adu.f <- adu.s+adu.r

# number of surviving juveniles
juv.s.term <- juv.n*juv.s
# number of surviving adults
adu.s.term <- adu.n*adu.s
# number of reproducing adults
adu.r.term <- adu.n*adu.r
# total adult fitness
adu.fitness <- adu.s.term+adu.r.term

di.term <- juv.s.term + adu.s.term
dd.term <- adu.r.term

rbPal <- colorRampPalette(c('blue','red'))
Col <- rbPal(10)[as.numeric(cut(outp[,4],breaks=10))]

quartz()
par(mfrow=c(3,3),mar=c(4,4,1,1),oma=c(1,1,1,1))
plot(outp[,3],outp[,4],pch=16,col=Col,xlab='Size at sexual maturity',ylab='Carrying capacity (strategy fitness)',cex=2)
place.letter('tl','(A)',1.5)
abline(v=outp[8,3],col='green',lty=2,lwd=2)

plot(outp[,3],juv.n,ylim=c(100,800),col=Col,pch=16,cex=2,xlab='Size at sexual maturity',ylab='Number of iindividuals')
points(outp[,3],adu.n,col=Col,cex=2,pch=17)
leg.txt <- c('Juveniles','Adults')
leg.pch <- c(16,17)
leg.col <- c('blue','blue')
legend('topleft',legend=leg.txt,pch=leg.pch,bty='n',cex=1.5,col=leg.col)
place.letter('br','(B)',1.5)
abline(v=outp[8,3],col='green',lty=2,lwd=2)

plot(juv.n,outp[,4],xlim=c(150,750),col=Col,cex=2,pch=17,xlab='Number of individuals',ylab='Carrying capacity')
points(adu.n,outp[,4],col=Col,cex=2,pch=16)
place.letter('tl','(C)',1.5)

plot(outp[,3],juv.s,ylim=c(0.2,1.5),col=Col,pch=16,cex=2,xlab='Size at sexual maturity',ylab='Per capita per time step rate')
points(outp[,3],adu.f,col=Col,pch=17,cex=2)
points(outp[,3],adu.s,col=Col,pch=2,cex=1.5)
points(outp[,3],adu.r,col=Col,pch=2,cex=1.5)
abline(v=outp[8,3],col='green',lty=2,lwd=2)
place.letter('tl','(D)',1.5)

plot(outp[,3],outp[,7],xlab='Size at sexual maturity',ylab='Proportion surviving to maturity',pch=16,cex=2,col=Col)
abline(v=outp[8,3],col='green',lty=2,lwd=2)
place.letter('tr','(E)',1.5)

plot(outp[,3],leza,xlab='Size at sexual maturity',ylab='Life expectancy at maturity',pch=16,cex=2,col=Col)
abline(v=outp[8,3],col='green',lty=2,lwd=2)
place.letter('tl','(F)',1.5)











r0.res <- rep(NA,20)
for (j in 1:20){
sad <- apply(struct[,,j],1,sum)
sad[z<outp[j,3]] <- 0
sad <- sad/sum(sad)
r0 <- 0
for (i in 1:300){
  r0 <- r0 + sum(r.mats[,,j]%*%sad)
  sad <- s.mats[,,j]%*%sad
}
r0.res[j] <- r0
}

plot(outp[,3],r0.res,xlab='Size at sexual maturity',ylab='Adult LRS',pch=16,cex=2,col=Col)
abline(v=outp[8,3],col='green',lty=2,lwd=2)
place.letter('tl','(G)',1.5)

plot(outp[,3],r0*outp[,7],pch=16,col=Col,cex=2,xlab='Size at sexual maturity',ylab='Juvenile survvial x adult fitness')
abline(v=outp[8,3],col='green',lty=2,lwd=2)
place.letter('tr','(H)',1.5)

plot(r0*outp[,7],outp[,4],pch=16,col=Col,cex=2,xlab='Proportion achieving maturity x adult fitness',ylab='Carrying capacity')
place.letter('tr','(I)',1.5)


rbPal <- colorRampPalette(c('blue','red'))
Col <- rbPal(10)[as.numeric(cut(outp[,4],breaks=10))]
cexy <- seq(1,1.75,length.out=20)
quartz()
sorter <- 1:20
sorter[14] <- 21
terms <- cbind(outp[,4],juv.n,adu.n,juv.s,adu.s,adu.r,sorter)
terms <- terms[order(sorter),]
Col[14:19] <- Col[15:20]
cexy[20] <- 3
Col[20] <- 'green'

labs <- c('K',expression('N'['J']),expression('N'['A']),expression(bar('S'['J'])),expression(bar('S'['A'])),expression(paste(bar('R'['A']),'(K)')))
pairs(terms[,1:6],col=Col,pch=16,labels=labs,cex=cexy,lower.panel=NULL)


