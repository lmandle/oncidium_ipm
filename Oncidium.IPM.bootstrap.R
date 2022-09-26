#########IPM for Oncidium brachyandrum in Oaxaca######

#load the libraries
library(lattice)
library(MASS)
library(popbio)
library(gridBase)
library(ggplot2)
library(glmmTMB)
library(nlme)
library(MuMIn)
library(tidyverse)

#import the file
onc<-read.table("oncidium.txt", header=T, sep="\t")

#check the file
head(onc)
str(onc)
length(onc)

#change host number to factor
onc$host_number<-as.factor(onc$host_number)
onc$host_species<-as.factor(onc$host_species)
onc$year<-as.factor(onc$year)

###########SURVIVAL############

#Remove the NAs
oncs<-subset(onc, initial.parea>0)

surv10<-glmmTMB(survival~log(initial.parea)+year*host_species+(1|host_number/id),family="binomial",data=oncs)
summary(surv10)

s1<-fixef(surv10)


############GROWTH##########

oncf<-subset(oncs,final.parea>0)

#added id
grow<-glmmTMB(log(final.parea)~log(initial.parea)+host_species*year+(1|host_number/id),data=oncf)
summary(grow)

g1<-fixef(grow)

###GROWTH VARIANCE
g.sig2<-summary(grow)$sigma^2     ##overall variance


#### REPRODUCTION####
###PROBABILITY OF REPRODUCTION

#Limit to plants or reproductive size
oncsr<-subset(oncs, initial.parea>0.7854)

rep3<-glmmTMB(reproduce~log(initial.parea)+host_species+(1|host_number/id),family="binomial",data=oncsr)
summary(rep3)

r1<-fixef(rep3)


###NUMBER OF CAPSULES
oncsc<-subset(oncs, capsules>0)

cap11b<-glmmTMB(capsules~log(initial.parea), family="nbinom1", data=oncsc)
summary(cap11b)

c1<-fixef(cap11b)


## ESTIMATING THE CAPSULE/SEEDLING RATIO "p.est".

#number of plantulas in 2018+2019+2020/no. total capsules in 2017+2018+2019

sd.cap.m<-(23+15+8)/(23+14+21) #0.79 mart
sd.cap.r<- (8+2+3)/(1+5+5) #1.18 rug

###SIZE DISTRIBUTION OF NEW SEEDLINGS

sdlg.mean.m<-log(0.12)
sdlg.mean.r<-log(0.12)

#sdlg.std - std of logged values

sdlg.std.m<-0.40
sdlg.std.r<-0.38


##IPM FUNCTION##################3
nhost_species<-length(unique(onc$host_species))
nyears<-length(unique(onc$year))
ncoef<-6
p.vec<-array(0,c(4,ncoef,nhost_species,nyears))

p.vec[1,1,1,1]<- s1$cond[1] #intercept for survival mart yr1
p.vec[1,1,1,2]<- s1$cond[1]+s1$cond[3] #mart yr2
p.vec[1,1,1,3]<- s1$cond[1]+s1$cond[4] #mar yr3
p.vec[1,1,2,1]<- s1$cond[1]+s1$cond[5] #rug yr1
p.vec[1,1,2,2]<- s1$cond[1]+s1$cond[3] +s1$cond[5]+s1$cond[6]#rugyr2
p.vec[1,1,2,3]<- s1$cond[1]+s1$cond[4] +s1$cond[5]+s1$cond[7]#rugyr3


p.vec[1,2,,]<- rep(s1$cond[2],nyears*nhost_species) #slope for survival all years

p.vec[2,1,1,1]<- g1$cond[1]#intercept for growth mart yr1
p.vec[2,1,1,2]<- g1$cond[1]+g1$cond[4]#intercept for growth mart yr2
p.vec[2,1,1,3]<- g1$cond[1]+g1$cond[5]#intercept for growth mart yr3
p.vec[2,1,2,1]<- g1$cond[1]+g1$cond[3]#intercept for growth rug yr1
p.vec[2,1,2,2]<- g1$cond[1]+g1$cond[3]+g1$cond[4]+g1$cond[6]#intercept for growth rug yr2
p.vec[2,1,2,3]<- g1$cond[1]+g1$cond[3]+g1$cond[5]+g1$cond[7]#intercept for growth rug yr3

p.vec[2,2,,]<- rep(g1$cond[2],nyears*nhost_species)#slope for growth
p.vec[2,3,,]<-rep(g.sig2,nyears*nhost_species)#g.sigma2 overall variation

#I changed this bc model was wrong
p.vec[3,1,1,1]<- r1$cond[1] #intercept for prob fruiting mart yr1
p.vec[3,1,1,2]<- r1$cond[1] # mart yr2
p.vec[3,1,1,3]<- r1$cond[1] # mart yr3
p.vec[3,1,2,1]<- r1$cond[1]+r1$cond[3] #  rug yr1
p.vec[3,1,2,2]<- r1$cond[1]+r1$cond[3] #  rug yr2
p.vec[3,1,2,3]<- r1$cond[1]+r1$cond[3] #  rug yr3

p.vec[3,2,,]<- rep(r1$cond[2],nyears*nhost_species)# slope for prob fruiting


p.vec[4,1,,]<- rep(c1$cond[1],nyears)#intercept for #caps/fruiting indiv, mart and rug allyrs
p.vec[4,2,,]<- rep(c1$cond[2],nyears*nhost_species)#slope for no.caps/fruit individuals

p.vec[4,4,1,]<- rep(sd.cap.m,nyears)#mean.sdlgs/capsule, mart years
p.vec[4,4,2,]<- rep(sd.cap.r,nyears)#mean.sdlgs/capsule, rug years


p.vec[4,5,1,]<- rep(sdlg.mean.m,nyears) #mean seedling size
p.vec[4,5,2,]<- rep(sdlg.mean.r,nyears) #mean seedling size

p.vec[4,6,1,]<- rep(sdlg.std.m,nyears) #mean seedling SD
p.vec[4,6,2,]<- rep(sdlg.std.r,nyears) #mean seedling S

p.vec


## The IPM FUNCTIONS: s(x), g(y,x), f(y,x), p(y,x), K(y,x)####### 
####ADULTS
####A.SURVIVAL function s(x)

sx<-function(x, pvec,host_species,year) {
  xbeta<-pvec[1,1,host_species, year] +pvec[1,2,host_species,year]*x;
  s<-exp(xbeta)/(1+exp(xbeta))
  return(s);
}


#graph to make sure function is correct
plot(sx(seq(log(0.01),log(10), length.out=50),p.vec,1,1)~(seq(log(0.01),log(10),length.out=50)), type="l",lty=1)
points(sx(seq(log(0.01),log(10), length.out=50),p.vec,1,2)~(seq(log(0.01),log(10),length.out=50)), type="l",lty=2)
points(sx(seq(log(0.01),log(10), length.out=50),p.vec,1,3)~(seq(log(0.01),log(10),length.out=50)), type="l",lty=3)
points(sx(seq(log(0.01),log(10), length.out=50),p.vec,2,1)~(seq(log(0.01),log(10),length.out=50)), type="l",col="blue",lty=1)
points(sx(seq(log(0.01),log(10), length.out=50),p.vec,2,2)~(seq(log(0.01),log(10),length.out=50)), type="l",col="blue",lty=2)
points(sx(seq(log(0.01),log(10), length.out=50),p.vec,2,3)~(seq(log(0.01),log(10),length.out=50)), type="l",col="blue",lty=3)


### B. GROWTH function g(y,x) 
gyx<-function(y, x, pvec, host_species, year) {
  mux<-pvec[2,1,host_species,year] + pvec[2,2,host_species,year]*x 
  sigmax2<-pvec[2,3,host_species,year]
  sigmax<-sqrt(sigmax2)
  g<-dnorm(y, mux, sigmax)
  return(g)
}

##Check function with graph
gx<-function(x, pvec, host_species, year) {
  mux<-pvec[2,1,host_species,year] + pvec[2,2,host_species,year]*x 
  return(mux)
}

#graph
plot(gx(seq(log(0.01),log(10),length.out=50),p.vec,1,1)~(seq(log(0.01),log(10),length.out=50)), type="l",col="black",lty=1,lwd=2,xlab="Adult parea a t (log)", ylab="Adult parea at time (t+1)")
points(gx(seq(log(0.01),log(10),length.out=50),p.vec,1,2)~(seq(log(0.01),log(10),length.out=50)), type="l",col="black",lty=2,lwd=2)
points(gx(seq(log(0.01),log(10),length.out=50),p.vec,1,3)~(seq(log(0.01),log(10),length.out=50)), type="l",col="black",lty=3,lwd=2)
points(gx(seq(log(0.01),log(10),length.out=50),p.vec,2,1)~(seq(log(0.01),log(10),length.out=50)), type="l",col="blue",lty=1,lwd=2)
points(gx(seq(log(0.01),log(10),length.out=50),p.vec,2,2)~(seq(log(0.01),log(10),length.out=50)), type="l",col="blue",lty=2,lwd=2)
points(gx(seq(log(0.01),log(10),length.out=50),p.vec,2,3)~(seq(log(0.01),log(10),length.out=50)), type="l",col="blue",lty=3,lwd=2)

points(log(oncf$initial.parea),log(oncf$final.parea), pch=as.numeric(oncf$host_species), col=as.numeric(oncf$year),xlab="Size parea at t (log)", ylab="Size parea at t+1 (log)")
legend(0.3, 4.3, c("Year1", "Year2", "Year3"), bg="white", bty="n", lwd=c(3,3,3),lty=c(1,2,3))


####C.The SURVIVAL-GROWTH function P(y, x)
pyx<-function(y,x, pvec, host_species, year) { 
  p<-sx(x, pvec,host_species,year)*gyx(y, x, pvec,host_species, year)
  return(p) 
}


###D.FERTILITY FUNCTION
fyx<-function(y, x, pvec,host_species,year) {
  xbeta<-pvec[3,1,host_species, year] +pvec[3,2,host_species,year]*x;
  p.cap<-ifelse(x<log(0.7854),0,(exp(xbeta)/(1+exp(xbeta))))
  no.cap<-exp(pvec[4,1,host_species,year] + pvec[4,2,host_species,year]*x); ##log(mu)=a+bx
  no.sdls.cap<-pvec[4,4,host_species,year]
  scd<-dnorm(y, pvec[4,5,host_species,year], pvec[4,6,host_species,year])
  f<-p.cap*no.cap*no.sdls.cap*scd    
  return(f)
}


#check graph wrt fruit production
fx<-function(x, pvec,host_species,year) {
  xbeta<-pvec[3,1,host_species, year] +pvec[3,2,host_species,year]*x;
  p.cap<-ifelse(x<log(0.7854),0,(exp(xbeta)/(1+exp(xbeta))))
  no.cap<-exp(pvec[4,1,host_species,year] + pvec[4,2,host_species,year]*x); ##log(mu)=a+bx
  f<-p.cap*no.cap
  return(f)
}  


plot(fx(seq(log(0.7),log(10),length.out=50),p.vec,1,1)~(seq(log(0.7),log(10),length.out=50)), type="l",col="black",lty=1,lwd=2,xlab="Size parea at (log)", ylab="No.fruits produced")
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,1,2)~(seq(log(0.7),log(10),length.out=50)), type="l",col="black",lty=2,lwd=2)
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,1,3)~(seq(log(0.7),log(10),length.out=50)), type="l",col="black",lty=3,lwd=2)
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,2,1)~(seq(log(0.7),log(10),length.out=50)), type="l",col="blue",lty=1,lwd=2)
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,2,2)~(seq(log(0.7),log(10),length.out=50)), type="l",col="blue",lty=2,lwd=2)
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,2,3)~(seq(log(0.7),log(10),length.out=50)), type="l",col="blue",lty=3,lwd=2)

points(log(oncsc$initial.parea), oncsc$capsules, type="p",lwd=1,pch=as.numeric(oncsc$host_species))

#Check graph - this is mean number of seedlings produced

fx<-function(x, pvec,host_species,year) {
  xbeta<-pvec[3,1,host_species, year] +pvec[3,2,host_species,year]*x;
  p.cap<-ifelse(x<log(1.44),0,(exp(xbeta)/(1+exp(xbeta))))
  no.cap<-exp(pvec[4,1,host_species,year] + pvec[4,2,host_species,year]*x); ##log(mu)=a+bx
  no.sdls.cap<-pvec[4,4,host_species,year]
  f<-p.cap*no.cap*no.sdls.cap   
  return(f)
}


plot(fx(seq(log(0.7),log(10),length.out=50),p.vec,1,1)~(seq(log(0.7),log(10),length.out=50)), type="l",col="black",lty=1,lwd=2,xlab="Size parea at (log)", ylab="No.fruits produced")
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,1,2)~(seq(log(0.7),log(10),length.out=50)), type="l",col="black",lty=2,lwd=2)
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,1,3)~(seq(log(0.7),log(10),length.out=50)), type="l",col="black",lty=3,lwd=2)
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,2,1)~(seq(log(0.7),log(10),length.out=50)), type="l",col="blue",lty=1,lwd=2)
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,2,2)~(seq(log(0.7),log(10),length.out=50)), type="l",col="blue",lty=2,lwd=2)
points(fx(seq(log(0.7),log(10),length.out=50),p.vec,2,3)~(seq(log(0.7),log(10),length.out=50)), type="l",col="blue",lty=3,lwd=2)

##Looking at probability density of new seedlings
plot(fyx(seq(-2,2,length.out=50),seq(2.5,3.5,length.out=50),p.vec,1,1)~(seq(-2,0.5,length.out=50)), type="l",col="black")
points(fyx(seq(-2,2,length.out=50),seq(2.5,3.5,length.out=50),p.vec,1,2)~(seq(-2,0.5,length.out=50)), type="l",col="black")
points(fyx(seq(-2,2,length.out=50),seq(2.5,3.5,length.out=50),p.vec,1,3)~(seq(-2,0.5,length.out=50)), type="l",col="black")
points(fyx(seq(-2,2,length.out=50),seq(2.5,3.5,length.out=50),p.vec,2,1)~(seq(-2,0.5,length.out=50)), type="l",col="green")
points(fyx(seq(-2,2,length.out=50),seq(2.5,3.5,length.out=50),p.vec,2,2)~(seq(-2,0.5,length.out=50)), type="l",col="green")
points(fyx(seq(-2,2,length.out=50),seq(2.5,3.5,length.out=50),p.vec,2,3)~(seq(-2,0.5,length.out=50)), type="l",col="green")


## The (master) KERNEL: K(y,x)= p(y,x) + f(y,x) + c(y,x)
Kyx<-function(y, x, pvec, host_species,year) {
  k<-pyx(y, x, pvec,host_species,year)+fyx(y, x, pvec, host_species,year)
  return(k) 
}


bigmat<-function(bigM, pvec, host_species, year){
  ## Set matrix size and convergence tolerance 
  min.sz<-1*min((range(log(onc$initial.parea),na.rm=T)))#check to remove outlier #0.9
  max.sz<-1*max((range(log(onc$initial.parea),na.rm=T)))#check to remove outlier #1.1
  
  # Compute meshpoints iteration matrix KD 
  h=(max.sz-min.sz)/(bigM+1); 
  y=(seq(min.sz,max.sz,length=bigM)+seq(min.sz+h,max.sz+h,length=bigM))/2;  
  
  ## Apply Kyx funct for y and y, pvec (=p.vec), species (choose!!)
  K=outer(y,y,Kyx, pvec, host_species, year);
  KD=h*K;
  return(KD);  		  ## This is your have the matrix
}

############################# PART (7) ##########################
## Eigenvalues and eigenvectors analysis, sensitivity,...

##Loop to extract lambda values
lam<-matrix(0,nhost_species,nyears)
A<-array(0, c(800, 800, nhost_species,nyears)) 

for(i in 1:nhost_species){
  for(j in 1:nyears){
    library(popbio)
    #A[,,i,j]<-bigmat(500, pvec=p.vec, year=i,site=j)
    lam[i,j]<-eigen.analysis(bigmat(800, pvec=p.vec, host_species=i,year=j))$lambda1
  }
}


#Extract lambda values
lam[1,1]#mart year 1 0.84
lam[1,2]#mart, year2 0.83 
lam[1,3]#mart, year3 0.88

lam[2,1]#rug yr1 0.82
lam[2,2]#rug yr2 0.97
lam[2,3]#rug yr3 1.01



#########################################################
### Bootstrap lambda (code adapted from Kuss et al. 2008)
#########################################################

n.boot=100 #Kuss used 5000
nhost_species<-2
nyears<-2
  lambda.boot=array(NA,c(nhost_species,nyears,n.boot))
  for(b.samp in 1:n.boot){
  
    #Adult survival
  sample.boot=c(sample(1:nrow(oncs),replace=T)) #generate bootstrapped sample
  oncs.boot<-data.frame( survival=oncs$survival[sample.boot],#create bootstrapped dataset
                           initial.parea=oncs$initial.parea[sample.boot],
                           year=oncs$year[sample.boot],
                           host_species=oncs$host_species[sample.boot],
                           host_number=oncs$host_number[sample.boot],
                           id=oncs$id[sample.boot])
  
  surv10.boot<-update(surv10, data=oncs.boot) #refit model
  s1.boot<-fixef(surv10.boot)
  
  
  #Adult growth
  sample.boot=c(sample(1:nrow(oncf),replace=T)) #generate bootstrapped sample
  oncf.boot<-data.frame(final.parea=oncf$final.parea[sample.boot],#create bootstrapped dataset
                        initial.parea=oncf$initial.parea[sample.boot],
                        year=oncf$year[sample.boot],
                        host_species=oncf$host_species[sample.boot],
                        host_number=oncf$host_number[sample.boot],
                        id=oncf$id[sample.boot])
  grow.boot<-update(grow, data=oncf.boot) #refit model
  g1.boot.boot<-fixef(grow.boot)
  g.sig2.boot<-summary(grow.boot)$sigma^2 
  
  
  #Fecundity
  sample.boot=c(sample(1:nrow(oncsr),replace=T))
  oncsr.boot<-data.frame(reproduce=oncsr$reproduce[sample.boot], #create bootstrapped dataset
                        initial.parea=oncsr$initial.parea[sample.boot],
                        host_species=oncsr$host_species[sample.boot],
                        host_number=oncsr$host_number[sample.boot],
                        id=oncsr$id[sample.boot])
  
  rep3.boot<-glmmTMB(reproduce~log(initial.parea)+host_species, data=oncsr.boot) #remove individual-level random effect
  r1.boot.boot<-coef(rep3.boot)
  
  #Capsules produced
  
  sample.boot=c(sample(1:nrow(oncsc),replace=T))
  oncsr.boot<-data.frame(capsules=oncsc$capsules[sample.boot], #create bootstrapped dataset
                         initial.parea=oncsc$initial.parea[sample.boot])
                         
  cap11b.boot<-update(cap11b, data=oncsr.boot) #refit model
  c1.boot<-fixef(cap11b.boot)
  
  
  oncsc<-subset(oncs, capsules>0)
  
  cap11b<-glmmTMB(capsules~log(initial.parea), family="nbinom1", data=oncsc)
  summary(cap11b)
  
  c1<-fixef(cap11b)
  
  
  #rebuild p.vec from bootstrapped params
  p.vec.boot<-array(0,c(4,ncoef,nhost_species,nyears))
  
    #params constant across scenarios
  p.vec.boot.[1,1,1,1]<- s1.boot$cond[1] #intercept for survival mart yr1.boot
  p.vec.boot.boot[1,1,1,2]<- s1.boot$cond[1]+s1.boot$cond[3] #mart yr2
  p.vec.boot[1,1,1,3]<- s1.boot$cond[1]+s1.boot$cond[4] #mar yr3
  p.vec.boot[1,1,2,1]<- s1.boot$cond[1]+s1.boot$cond[5] #rug yr1.boot
  p.vec.boot[1,1,2,2]<- s1.boot$cond[1]+s1.boot$cond[3] +s1.boot$cond[5]+s1.boot$cond[6]#rugyr2
  p.vec.boot[1,1,2,3]<- s1.boot$cond[1]+s1.boot$cond[4] +s1.boot$cond[5]+s1.boot$cond[7]#rugyr3
  
  
  p.vec.boot[1,2,,]<- rep(s1.boot$cond[2],nyears*nhost_species) #slope for survival all years
  
  p.vec.boot[2,1,1,1]<- g1.boot.boot$cond[1]#intercept for growth mart yr1.boot
  p.vec.boot[2,1,1,2]<- g1.boot.boot$cond[1]+g1.boot$cond[4]#intercept for growth mart yr2
  p.vec.boot[2,1,1,3]<- g1.boot$cond[1]+g1.boot$cond[5]#intercept for growth mart yr3
  p.vec.boot[2,1,2,1]<- g1.boot$cond[1]+g1.boot$cond[3]#intercept for growth rug yr1.boot
  p.vec.boot[2,1,2,2]<- g1.boot$cond[1]+g1.boot$cond[3]+g1.boot$cond[4]+g1.boot$cond[6]#intercept for growth rug yr2
  p.vec.boot[2,1,2,3]<- g1.boot$cond[1]+g1.boot$cond[3]+g1.boot$cond[5]+g1.boot$cond[7]#intercept for growth rug yr3
  
  p.vec.boot[2,2,,]<- rep(g1.boot$cond[2],nyears*nhost_species)#slope for growth
  p.vec.boot[2,3,,]<-rep(g.sig2.boot,nyears*nhost_species)#g.sigma2 overall variation
  
  p.vec.boot[3,1,1,1]<- r1.boot$cond[1] #intercept for prob fruiting mart yr1.boot
  p.vec.boot[3,1,1,2]<- r1.boot$cond[1] # mart yr2
  p.vec.boot[3,1,1,3]<- r1.boot$cond[1] # mart yr3
  p.vec.boot[3,1,2,1]<- r1.boot$cond[1]+r1.boot$cond[3] #  rug yr1.boot
  p.vec.boot[3,1,2,2]<- r1.boot$cond[1]+r1.boot$cond[3] #  rug yr2
  p.vec.boot[3,1,2,3]<- r1.boot$cond[1]+r1.boot$cond[3] #  rug yr3
  
  p.vec.boot[3,2,,]<- rep(r1.boot$cond[2],nyears*nhost_species)# slope for prob fruiting
  
  
  p.vec.boot[4,1,,]<- rep(c1.boot$cond[1],nyears)#intercept for #caps/fruiting indiv, mart and rug allyrs
  p.vec.boot[4,2,,]<- rep(c1.boot$cond[2],nyears*nhost_species)#slope for no.caps/fruit individuals
  
  p.vec.boot[4,4,1,]<- rep(sd.cap.m,nyears)#mean.sdlgs/capsule, mart years
  p.vec.boot[4,4,2,]<- rep(sd.cap.r,nyears)#mean.sdlgs/capsule, rug years
  
  
  p.vec.boot[4,5,1,]<- rep(sdlg.mean.m,nyears) #mean seedling size
  p.vec.boot[4,5,2,]<- rep(sdlg.mean.r,nyears) #mean seedling size
  
  p.vec.boot[4,6,1,]<- rep(sdlg.std.m,nyears) #mean seedling SD
  p.vec.boot[4,6,2,]<- rep(sdlg.std.r,nyears) #mean seedling S
  
 
  #calculate bootstrapped lambda across all scenarios
  ##LISA NEXT 2 PARTS IS WHERE I NEED HELP W THE CODE! :)
  
  for(i in 1:nhost_species){
    for(j in 1:nyears){
      lambda.boot[i,j,b.samp]<-tryCatch(
      expr<-eigen.analysis(bigmat(400, pvec=p.vec.boot[host_species=i,year=j]))$lambda1,
      error=function(cond){
        message(cond)
        return(NA)},
      warning=function(cond){
        message(cond)
        return(NULL)})
    }
  }
  }
 
  #convert lambda.boot array to dataframe
lambda.boot.df<-data.frame(scenario=rep(c(1:n_host*nyears), each=n.boot),
                           rep=rep(c(1:n.boot), teams=nhost_species*nyears),
                           lambda=as.vector(t(lambda.boot))) #turns lambda.boot lambda values into a vector, ordered by scenario

date = gsub(":","-",Sys.time()) #get date and time to append to filename  
date = gsub(" ","_",date)
write.csv(lambda.boot.df, file=paste0("lama.lambda.boot", "_", date, ".csv"))

lambda.summary <- lambda.boot.df %>%
  group_by(host_species,year) %>%
  summarize(mean_l = mean(lambda, na.rm=T), median_l=median(lambda,na.rm=T),
            lower_95ci=quantile(lambda, p=0.025, na.rm=T), 
            upper_95ci=quantile(lambda, p=0.975, na.rm=T))

ggplot(data=lambda.boot.df, aes(lambda))+  
  geom_histogram(bins=50)+
  facet_grid(cols=vars(scenario))+
  geom_vline(data  = lambda.summary, aes(xintercept = mean_l), color = "green")+
  geom_vline(data  = lambda.summary, aes(xintercept = median_l), color = "blue")+
  geom_vline(data  = lambda.summary, aes(xintercept = lower_95ci), color = "red")+
  geom_vline(data  = lambda.summary, aes(xintercept = upper_95ci), color = "red")
