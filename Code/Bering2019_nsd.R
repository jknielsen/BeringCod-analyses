
# Fit NSD data from migrating fish to dispersal NSD curve

#First, calculate NSD for all 
library(lubridate)

vit.2019<- read.csv("Z:/_Research/Projects/pcod/Analysis/R/pcod/Bering/BS_2019_Viterbi.csv")
vit.2019$date<- as.Date(vit.2019$date)

#For each tag
tags.2019<- unlist(dimnames(table(vit.2019$Tag)))
i=1
for (i in 1:length(tags.2019)){
  tag.temp<- vit.2019[vit.2019$Tag== tags.2019[i],]
  tag.X.0<- tag.temp$X[1]
  tag.Y.0<- tag.temp$Y[1]
  vit.2019$disp[vit.2019$Tag== tags.2019[i]]<- sqrt((tag.temp$X-tag.X.0)^2 + (tag.temp$Y-tag.Y.0)^2)
  vit.2019$rec.mo[vit.2019$Tag== tags.2019[i]]<- tag.temp$month[length(tag.temp$month)]
  vit.2019$rec.yr[vit.2019$Tag== tags.2019[i]]<- year(tag.temp$date[length(tag.temp$date)])
  vit.2019$DAL[vit.2019$Tag== tags.2019[i]]<- 1:nrow(tag.temp)
}

vit.2019$r2n.km<- (vit.2019$disp/1000)^2
vit.2019$days.datum<- as.numeric(vit.2019$date - min(vit.2019$date))

conv.tags<- read.csv("Bering2019_conv_recap.csv",header=T)
conv.tags$REL_TIME_UTC<- as.POSIXct(conv.tags$REL_TIME_UTC, tz="GMT", format="%m/%d/%Y %H:%M")
conv.tags$REC_TIME_UTC<- as.POSIXct(conv.tags$REC_TIME_UTC, tz="GMT", format="%m/%d/%Y %H:%M")
conv.tags$Rel_date<- as.Date(conv.tags$REL_TIME_UTC)
conv.tags$Rec_date<- as.Date(conv.tags$REC_TIME_UTC)

library(terra)
rel.locs.pr<- terra::project(vect(cbind(conv.tags$REL_LONGDD,conv.tags$REL_LATDD),type="points", crs="+proj=longlat +datum=WGS84"),
                             "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs ")
rec.locs.pr<- terra::project(vect(cbind(conv.tags$REC_LNG,conv.tags$REC_LAT),type="points", crs="+proj=longlat +datum=WGS84"),
                             "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs ")

conv.tags$X.0<- geom(rel.locs.pr)[,3]
conv.tags$Y.0<- geom(rel.locs.pr)[,4]
conv.tags$X.end<- geom(rec.locs.pr)[,3]
conv.tags$Y.end<- geom(rec.locs.pr)[,4]
conv.tags$disp<- sqrt((conv.tags$X.end-conv.tags$X.0)^2 + (conv.tags$Y.end-conv.tags$Y.0)^2)
conv.tags$r2n.km<- (conv.tags$disp/1000)^2

conv.tags$days.datum<- as.numeric(conv.tags$Rec_date - min(vit.2019$date))


bering2019.nsd<- vit.2019[,c("Tag","days.datum","r2n.km")]
bering2019.nsd$Tag<- factor(bering2019.nsd$Tag)

conv.tag.nsd<- conv.tags[,c("TagID","days.datum","r2n.km")]
conv.tag.nsd$TagID<- factor(conv.tag.nsd$TagID)
dimnames(conv.tag.nsd)[[2]][1]<- "Tag"
bering2019.nsd<- rbind(bering2019.nsd,conv.tag.nsd)

summer.2019<- vit.2019$rec.mo < 12 & vit.2019$rec.yr == "2019"
winter.2020<- vit.2019$rec.mo < 5 & vit.2019$rec.yr == "2020"
summer.2020<- vit.2019$rec.mo > 4 & vit.2019$rec.yr == "2020"
winter.tags<- vit.2019[vit.2019$rec.yr=="2020",]

plot(vit.2019$days.datum,vit.2019$r2n.km,ylab="Net Squared Displacement (km^2)",type="n",xlab="# days since 8/1/2019")
points(vit.2019$days.datum[winter.2020],vit.2019$r2n.km[winter.2020],col="blue",pch=19,cex=.4)
points(vit.2019$days.datum[summer.2020],vit.2019$r2n.km[summer.2020],col="cyan",pch=19,cex=.4)
points(vit.2019$days.datum[summer.2019],vit.2019$r2n.km[summer.2019],col="red",pch=19,cex=.4)
points(conv.tags$days.datum,conv.tags$r2n.km,pch=19)
legend("topleft",c("Summer 2019 pop-up","Winter 2020 pop-up","Summer 2020 pop-up",
                   "End Location only"),
       col=c("red","blue","cyan","black"),pch=19,pt.cex=c(.5,.5,.5,1),cex=.8)
abline(v = 225) #Cut off for dispersal model

summer.2019.nsd<- vit.2019[summer.2019,]
winter.2020.nsd<- vit.2019[winter.2020,]
summer.2020.nsd<- vit.2019[summer.2020,]

#Look at individual fish
tag.no<- "183962"
plot(summer.2019.nsd$days.datum,summer.2019.nsd$r2n.km)
points(summer.2019.nsd$days.datum[summer.2019.nsd$Tag==tag.no],
       summer.2019.nsd$r2n.km[summer.2019.nsd$Tag==tag.no],col="blue")

#########################################
#*** Fitting NSD curves to individual fish***
#########################################
#Remove data after days.datum 225 (after March 21, when the number of tagged fish drops off)
winter.tags<- winter.tags[winter.tags$days.datum < 226,]
winter.tag.ids<- unlist(dimnames(table(winter.tags$Tag)))

i=3
animal<- winter.tags[winter.tags$Tag == winter.tag.ids[i],]
plot(animal$days.datum,animal$r2n.km)

#-----------------------------------
# First model: home range (NSD ~ c) - 1 parameter
#----------------------------------
library(nlme)
mod1<- lm(r2n.km ~ 1, data=animal)
summary(mod1)
plot(animal$days.datum,animal$r2n.km,main=animal$Tag[1],xlab="Days since Aug 1",ylab="NSD (m^2)")
points(animal$days.datum,mod1$fitted, type="l",col="yellow",lwd=2)

plot(mod1)
sqrt(coef(mod1))

AIC(mod1)

#-------------------------------------
# Second model: random movement (linear relationship with time) - 2 parameters
#-----------------------------------
mod2<- lm(r2n.km ~ days.datum, data=animal)

summary(mod2)

plot(animal$days.datum,animal$r2n.km,main=animal$Tag[1],xlab="Days since Aug 1",ylab="NSD (m^2)")
points(animal$days.datum,mod2$fitted, type="l",col="red")
plot(mod2)

points(animal$days.datum,2*predict(mod2, se.fit=T)$se.fit+predict(mod2),type="l",lty=2,col="red")
points(animal$days.datum,-2*predict(mod2, se.fit=T)$se.fit+predict(mod2),type="l",lty=2,col="red")

AIC(mod1,mod2)

#-------------------------------------
# Third model: dispersal (non-linear) - 3 parameters
#-------------------------------------
NSD.disp <- function(day,asym,xmid,scale) {
  asym /(1+exp((xmid-day)/scale))
}

plot(animal$days.datum,animal$r2n.km,main=animal$Tag[1],xlab="Days since Aug 1",ylab="NSD (m^2)")

a1<- 600000
mid1<- 150
s1<- 30
curve(NSD.disp(x, asym=a1,xmid=mid1,scale=s1), 0, 225, col="red", lwd=2, add=T)

mod3<- nls(r2n.km ~   asym /(1+exp((xmid-days.datum)/scale)), 
           start = c(asym=a1,xmid=mid1,scale=s1),
           data=animal,control = list(maxiter = 500))

AIC(mod1,mod2,mod3)
summary(mod3)

plot(animal$days.datum,animal$r2n.km,main=animal$Tag[1],xlab="Days since Aug 1",ylab="NSD (m^2)")

mod3.out<- summary(mod3)$coef
curve(NSD.disp(x, asym=mod3.out[1],xmid=mod3.out[2],scale=mod3.out[3]), 0, 225, col="blue", lwd=2, add=T)
AIC(mod1,mod2,mod3)


################################################################################
#***Then fit mixed effects models for those fish that fit dispersal NSD model***
################################################################################

#Only fish that were able to fit dispersal model (model 3)
disp.fish<- winter.tag.ids[c(3,4,5,9,10,12,14,15,16,17)]
winter.disp<- winter.tags[winter.tags$Tag %in% disp.fish,]
winter.disp$Tag<- factor(winter.disp$Tag)

#Find starting values
#Parameters estimated individually with nls (not accounting for autocorrelation!)
indiv.asym<- c(641100, 509800, 277200, 377100, 313400, 396500, 268400, 189300)
hist(indiv.asym)
mean(indiv.asym)#371600

indiv.xmid<- c(154, 159, 166.5, 155.1, 141.8, 164.3, 129.1, 166.5)
hist(indiv.xmid)
mean(indiv.xmid) #154.5

indiv.scale<- c(28.59, 25.07, 18.38, 16.65, 23.3, 16.13, 6.96, 16.95)
hist(indiv.scale)
mean(indiv.scale)# 19.0

#############   ** Mixed effects model, dispersal **  ################
NSD.disp <- function(day,asym,xmid,scale) {
  asym /(1+exp((xmid-day)/scale))
}

a1<- 371600
mid1<- 150
s1<- 15
plot(winter.disp$days.datum,winter.disp$r2n.km)
curve(NSD.disp(x, asym=a1,xmid=mid1,scale=s1), 0, 225, col="red", lwd=2, add=T)

library(nlme)
Model3<-nlme(r2n.km ~  asym /(1+exp((mid1-days.datum)/s1)), 
             fixed = list(asym+mid1+s1~1),  
             random= asym ~ 1|Tag, 
             start = c(asym=a1,mid1=mid1,s1=s1),                  
             data=winter.disp,
             correlation=corAR1(form=~days.datum|Tag))

mod<- Model3
summary(mod)

plot(mod)
plot(mod, Tag ~ resid(.),abline=0)
qqnorm(mod)
qqnorm(mod, ~ranef(.))


plot(winter.disp$days.datum,winter.disp$r2n.km)
mod.out<- summary(mod)$coef$fixed
curve(NSD.disp(x, asym=mod.out[1],xmid=mod.out[2],scale=mod.out[3]), 0, 225, col="blue", lwd=2, add=T)

aX<- summary(mod)$tTable[1,1]
aX.se<- summary(mod)$tTable[1,2]

midX<- summary(mod)$tTable[2,1] #181 days
midX.se<- summary(mod)$tTable[2,2]

scaleX<- summary(mod)$tTable[3,1] #25.8 days
scaleX.se<- summary(mod)$tTable[3,2]

range(sqrt(coef(mod)[1])) #464-1166 km
sqrt(aX) #775 km 
sqrt(aX + 1.96* aX.se)#919 km
sqrt(aX - 1.96* aX.se)#597 km

## Migration parameters
#Scale estimates time to travel between midpoint and 3/4 of distance
#95% of the distance between midpoint and arrival occurs at 3 * scale (Borger and Fryxell 2012)

#Date of departure:
day.dep<- midX-(scaleX*3)
earliest<- (midX - 1.96*midX.se) - (scaleX + 1.96*scaleX.se)*3
latest<- (midX + 1.96*midX.se) - (scaleX - 1.96*scaleX.se)*3
as.Date("2019-08-01")+day.dep
as.Date("2019-08-01")+earliest
as.Date("2019-08-01")+latest

##Plot output
png("FigX_BeringNSD_rev.png",width = 6,height=5,res= 300,units="in",pointsize = 10)
par(mfrow=c(2,1))
par(mar=c(3,2,1,1))
plot(winter.disp$days.datum,winter.disp$r2n.km,type="n",xlab="",ylab="")
     #xlab="# days since 8/1/2019")#,
     #ylab="Net Squared Displacement (km^2)")
poly.winter.dispx<- c(seq(0,225,by =1),seq(225,0,by=-1))
poly.winter.dispy<- c((aX+(1.96*aX.se))/(1+exp(((midX-(1.96*midX.se))-poly.winter.dispx[1:226])/(scaleX+1.96*scaleX.se))),
             (aX-(1.96*aX.se))/(1+exp(((midX+(1.96*midX.se))-poly.winter.dispx[227:452])/(scaleX-1.96*scaleX.se))))
polygon(c(earliest,latest,latest,earliest),c(1300000,1300000,0,0),col=gray(.9),border=NA)
polygon(poly.winter.dispx,poly.winter.dispy,col="gray",border=NA)
points(winter.disp$days.datum,winter.disp$r2n.km,pch=19, cex=.3,col=gray(.3))
curve(NSD.disp(x, asym=aX,xmid=midX,scale=scaleX), 0, 225, col="black", lwd=2,lty=2,add=T)
abline(v = day.dep,lty=2)
legend("topleft",c("Fish NSD","NSD migr. curve (95% c.i.)","Migr. init. date (95% c.i.)"),
       pch=c(19,NA,NA),pt.cex=c(.3,NA,NA),lty=c(NA,2,2),lwd=c(NA,2,1),
       fill=c(NA,"gray",gray(.9)),border=NA,cex=.6)
box()

#All fish
plot(vit.2019$days.datum,vit.2019$r2n.km,ylab="",type="n",xlab="",xaxt="n")
polygon(c(earliest,latest,latest,earliest),c(1300000,1300000,0,0),col=gray(.9),border=NA)

polygon(poly.winter.dispx,poly.winter.dispy,col="gray",border=NA)

points(vit.2019$days.datum[winter.2020],vit.2019$r2n.km[winter.2020],col="blue",pch=19,cex=.2)
points(vit.2019$days.datum[summer.2020],vit.2019$r2n.km[summer.2020],col="cyan",pch=19,cex=.2)
points(vit.2019$days.datum[summer.2019],vit.2019$r2n.km[summer.2019],col="red",pch=19,cex=.2)
points(conv.tags$days.datum,conv.tags$r2n.km,pch=19)
abline(v = day.dep,lty=2)
curve(NSD.disp(x, asym=aX,xmid=midX,scale=scaleX), 0, 225, col="black", lwd=2,lty=2,add=T)
axis(side=1,at=c(0,31,61,92,122,153,184,213,244,274,305,335,366,397),labels = month.abb[c(8:12,1:9)])
box()
legend("topleft",c("Summer 2019 NSD","Winter 2020 NSD","Summer 2020 NSD",
                   "End Loc only","NSD migr curve (95% C.I.)","Migr. init. date (95% c.i.)"),
       col=c("red","blue","cyan","black","black","black"),pch=c(19,19,19,19,NA,NA),
       pt.cex=c(.5,.5,.5,1,NA,NA),cex=.6,lty=c(NA,NA,NA,NA,2,2),lwd=c(NA,NA,NA,NA,2,1),
       fill=c(NA,NA,NA,NA,"gray",gray(.9)),border=NA)
box()
dev.off()

