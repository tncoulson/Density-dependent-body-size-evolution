
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
juv.s.term <- (1-P)*juv.s
# number of surviving adults
adu.s.term <- P*adu.s
# number of reproducing adults
adu.r.term <- P*adu.r
# total adult fitness
adu.fitness <- adu.s.term+adu.r.term


quartz()

# code to draw the life histories
par(mfrow=c(4,3),mar=c(4,4,1,1),oma=c(1,1,1,1))

# define colours for each life history
rbPal <- colorRampPalette(c('blue','red'))
Col <- rbPal(10)[as.numeric(cut(outp[,4],breaks=10))]

# plot the inheritance function
contour(z,z,t(d),xlab='Parent body size t',ylab='Offspring body size t+1')
points(outp[,3],rep(4,20),pch=16,col=Col,cex=2)
place.letter('tl','(A)',1.5)
text(20,10,'Inheritance function')

# plot the growth function -- size t versus size t+1
outpz <- array(NA,c(20,3))
for (i in 1:20){
  grow.b <- grow.b.vals[i]
  grow.a <- 8-grow.b*4
  z.a <- (grow.a/(1-grow.b)) * 0.8
  outpz[i,] <- c(grow.a,grow.b,z.a)
}

# plot the growth function
plot(z,z,type='l',lty=2,xlab='Body size t',ylab='Body size t+1')
for(i in 1:20) abline(outpz[i,1],outpz[i,2],col=Col[i])
points(outp[,3],outp[,1]+outpz[,2]*outp[,3],pch=16,col=Col,cex=2)
place.letter('tl','(B)',1.5)
text(20,2,"Growth function")

# plot the von Bertalanffy versions of the growth functions
size.at.age <- array(NA,c(20,20))
size.at.age[1,] <- 4
for (i in 1:20){
  for (j in 2:20){
    size.at.age[j,i] <- outp[i,1]+outp[i,2]*size.at.age[j-1,i]
  }
}

plot(1:20,size.at.age[,20],col=Col[20],type='l',xlab='Age',ylab='Body size')
for (i in 19:1){ 
  lines(1:20,size.at.age[,i],col=Col[i])
}
place.letter('tl','(C)',1.5)

asymp.sz <- outp[,3]/0.8
K <- -log(outp[,2])
t0 <- 3+((1/K)*log((asymp.sz-size.at.age[3,])/asymp.sz))
asm <- t0 - log(0.2)/K
points(asm,outp[,3],pch=16,col=Col,cex=2)

# plot the size-survival function 
plot(z,surv(z,surv.a,surv.b,surv.c,outp[1,4]),type='l',xlab='Body size time t',ylab='Survival rate time t',ylim=c(0,1),col=Col[1])
for (i in 2:20) lines(z,surv(z,surv.a,surv.b,surv.c,outp[i,4]),col=Col[i])
place.letter('tl','(D)',1.5)
for (i in 1:20) points(outp[i,3],surv(outp[i,3],surv.a,surv.b,surv.c,outp[i,4]),pch=16,col=Col[i],cex=2)
text(20,0.06,"Survival function")

# Survivorship function
lx <- apply(surv.rate.by.age,2,cumprod)
lx <- rbind(rep(1,20),lx)
lx <- lx[-61,]
plot(0:59,lx[,1],xlab='Age',ylab='`Survivorship',type='l',ylim=c(0,0.2),xlim=c(0,40),col=Col[1])
for (i in 2:20) lines(0:59,lx[,i],col=Col[i])
place.letter('tr','(E)',1.5)

# survivorship to sexual maturity
plot(outp[,3],outp[,7],xlab='Size at sexual maturity',ylab='Probability of surviving to maturity',pch=16,cex=2,col=Col)
place.letter('tr','(F)',1.5)
abline(v=outp[14,3],lty=2,col='green',lwd=2)

# life expectancy at sexual maturity
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
abline(v=outp[14,3],lty=2,col='green',lwd=2)


# size at sexual maturity vs proportion population sexually mature
# proportion of population that are adults
P <- outp[,8]/outp[,4]
plot(outp[,3],P,col=Col,xlab='Size at sexual maturity',ylab='Proportion sexually mature',pch=16,cex=2)
place.letter('tr','(H)',1.5)
abline(v=outp[14,3],lty=2,col='green',lwd=2)
# size at sexual maturity vs reproductive output
plot(outp[,3],adu.r,col=Col,xlab='Size at sexual maturity',ylab='Per capita reproductive rate',pch=16,cex=2)
place.letter('tl','(I)',1.5)
abline(v=outp[14,3],lty=2,col='green',lwd=2)

plot(outp[,3],outp[,4],pch=16,col=Col,cex=2,xlab='Size at sexual maturity',ylab='Carrying capacity')
place.letter('tl','(J)',1.5)

plot(log(leza*outp[,7]),log(adu.r),col=Col,cex=2,pch=16,xlab=expression(paste('Log(L'['z'[a]],'* life expectancy at z'[a],')')),ylab='Log(reproductive rate)')
y <- log(adu.r)
x <- log(leza*outp[,7])
m1 <- lm(y~x)
summary(m1)
abline(m1)
place.letter('tr','(K)',1.5)
arrows(x0=-0.2,x1=0,y0=0.3,y1=0.1,length=0.15)

plot(log(leza*outp[,7]),(outp[,4]),col=Col,cex=2,pch=16,xlab=expression(paste('Log(L'['z'[a]],'* life expectancy at z'[a],')')),ylab='Carrying capacity')
y <- outp[,4]
m2 <- lm(y~x)
abline(m2)
place.letter('tl','(L)',1.5)









quartz()
par(mfrow=c(3,3),mar=c(4,4,1,1),oma=c(1,1,1,1))
plot(outp[,3],adu.n,pch=16,col=Col,cex=2)
plot(outp[,3],juv.n,pch=16,col=Col,cex=2)
plot(outp[,3],adu.s,pch=16,col=Col,cex=2)
plot(outp[,3],adu.r,pch=16,col=Col,cex=2)
plot(outp[,3],juv.s,pch=16,col=Col,cex=2)


quartz()
par(mfrow=c(3,3),mar=c(4,4,1,1),oma=c(1,1,1,1))
plot(1-P,outp[,4],pch=16,col=Col,cex=2,xlab='Proportion juvenile',ylab='Carrying capacity')
points(1-P[14],outp[14,4],col='green',pch=16,cex=2.5)
place.letter('tl','(A)',1.5)

plot(P,outp[,4],pch=16,col=Col,cex=2,xlab='Proportion adult',ylab='Carrying capacity')
points(P[14],outp[14,4],col='green',pch=16,cex=2.5)
place.letter('tl','(B)',1.5)

plot(juv.s,outp[,4],pch=16,col=Col,cex=2,xlab='Mean juvenile survival',ylab='Carrying capacity')
points(juv.s[14],outp[14,4],col='green',pch=16,cex=2.5)
place.letter('tl','(C)',1.5)

plot(adu.s,outp[,4],pch=16,col=Col,cex=2,xlab='Mean adult survival',ylab='Carrying capacity')
points(adu.s[14],outp[14,4],col='green',pch=16,cex=2.5)
place.letter('tl','(D)',1.5)

plot(adu.r,outp[,4],pch=16,col=Col,cex=2,xlab='Mean adult reproduction',ylab='Carrying capacity')
points(adu.r[14],outp[14,4],col='green',pch=16,cex=2.5)
place.letter('tr','(E)',1.5)

plot(juv.s,adu.s,pch=16,col=Col,cex=2,xlab='Mean juvenile survival',ylab='Mean adult survival')
points(juv.s[14],adu.s[14],col='green',pch=16,cex=2.5)
place.letter('tl','(F)',1.5)

plot(juv.s,adu.r,pch=16,col=Col,cex=2,xlab='Mean juvenile survival',ylab='Mean adult reproduction')
points(juv.s[14],adu.r[14],col='green',pch=16,cex=2.5)
place.letter('tl','(G)',1.5)

plot(adu.s,adu.r,pch=16,col=Col,cex=2,xlab='Mean adult survival',ylab='Mean adult reproduction')
points(adu.s[14],adu.r[14],col='green',pch=16,cex=2.5)
place.letter('tr','(H)',1.5)


pairs(cbind(juv.s.term,adu.s.term,adu.r.term,adu.fitness))
plot(juv.s.term,adu.fitness)
points(juv.s.term[14],adu.fitness[14],pch=16,col='green',cex=2.5)

plot(juv.s.term,outp[,4],pch=16,col=Col,cex=2,xlab='Juvenile survival term',ylab='Carrying capacity')
place.letter('tl','(F)',1.5)
plot(adu.s.term,outp[,4],pch=16,col=Col,cex=2,xlab='Adult survival term',ylab='Carrying capacity')
place.letter('tl','(G)',1.5)
plot(adu.r.term,outp[,4],pch=16,col=Col,cex=2,xlab='Adult recruitment term',ylab='Carrying capacity')
place.letter('tl','(H)',1.5)
plot(adu.fitness,outp[,4],pch=16,col=Col,cex=2,xlab='Adult demographic term',ylab='Carrying capacity')
place.letter('tl','(I)',1.5)



mx <- matrix(adu.r,nrow=60,ncol=20,byrow=TRUE)
asm.mat <- matrix(asm,nrow=60,ncol=20,byrow=TRUE)
age.mat <- matrix(1:60,nrow=60,ncol=20)
multiplier.up <- ifelse(age.mat<(asm.mat-1),0,1)
multiplier.down <- ifelse(age.mat>(asm.mat),0,1)
border.line <- ifelse(multiplier.up+multiplier.down==2,1,0)
repo <- 1-(asm-floor(asm))
repo <- matrix(repo,nrow=60,ncol=20,byrow=TRUE)
props <- border.line*repo
multiplier <- ifelse(age.mat<asm.mat,0,1)+props

lx.mx <- lx*mx*multiplier

age.c <- seq(1,25,length.out=1000)
rectangle <- array(NA,c(1000,20))
for (i in 1:20){
  rectangle[,i] <- ifelse(age.c>asm[i] & age.c<(asm[i]+leza[i]),(outp[i,7]*adu.r[i]),0)
}

quartz()
par(mfrow=c(2,2),mar=c(1,1,1,1),oma=c(1,1,1,1))
ribbon3D(x=1:20,y=outp[,3],z=lx.mx[1:20,],xlab='Age',ylab='Size at sexual matuirty',zlab='L(a)M(a)',xlim=c(1,20),theta=200,phi=10,colkey=FALSE)
#quartz()
ribbon3D(z=rectangle,theta=200,xlab='Age',ylab='Size at sexual maturity',zlab='L(a)M(a)',colkey=FALSE,phi=10)


pairs(cbind(outp[,3],leza,outp[,7],adu.r))

plot(outp[,7],leza)
points(outp[14,7],leza[14],col='green',pch=16,cex=2)
abline(0,1,col='red')
points(outp[14,7],1/leza[14],col='green',pch=16,cex=2)

quartz()
plot(outp[,3],outp[,7],pch=16,ylim=c(0,0.5))
points(outp[,3],1/leza,col='red',pch=16)
abline(v=outp[14,3],col='red',lty=2)

y <- (outp[,7]-mean(outp[,7]))/sd(outp[,7])
x <- outp[,3]
x2 <- outp[,3]^2
x3 <- outp[,3]^3
x4 <- outp[,3]^4
m1 <- lm(y~x+x2+x3+x4)
summary(m1)


y2 <- (leza-mean(leza))/sd(leza)
m2 <- lm(y2~x+x2+x3+x4)
summary(m2)

qq <- outp[15,3]
m1$coef[2]+m1$coef[3]*2*qq+m1$coef[4]*3*qq*qq+m1$coef[5]*4*qq*qq*qq
m2$coef[2]+m2$coef[3]*2*qq+m2$coef[4]*3*qq*qq+m2$coef[5]*4*qq*qq*qq






-3.1564476+0.1912330*2*outp[13,3]-0.0039279*3*outp[13,3]*outp[13,3]
-1.3552620+0.1256527*2*outp[13,3]-0.0051092*3*outp[13,3]*outp[13,3]
