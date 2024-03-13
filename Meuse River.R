rm(list=ls())

#install.packages("spData")

library(sp)

data(meuse.grid)
d<-meuse.grid
d1<-SpatialPoints(cbind(d$x,d$y))
summary(d1)
plot(d1) 

d2<-SpatialPixels(d1)
summary(d2)
plot(d2) 
image(d2,col="grey")

coordinates(d)<-c("x","y")
d<-as(d,"SpatialPixels")
image(d,col="grey")
title("grid")

data(meuse)
e<-meuse
coordinates(e) <- c("x", "y")
plot(e)
title("points")

data(meuse.riv)
b<-meuse.riv
b1<-list(Polygons(list(Polygon(b)),"b"))
b2<-SpatialPolygons(b1)
plot(b2,col="grey")
title("polygons")

image(d2,col="lightgrey")
plot(b2,col="grey",add=TRUE)
plot(e,add=TRUE)

layout(matrix(c(1, 2), 1, 2))
plot(b2,axes=TRUE)
plot(b2,axes=FALSE)
axis(1,at=c(178000 + 0:2 * 2000), cex.axis = 0.7)
axis(2,at=c(326000 + 0:3 * 4000), cex.axis = 0.7)
box()

oldpar = par(no.readonly = TRUE)
layout(matrix(c(1, 2), 1, 2))
plot(e,axes=TRUE,cex=0.6)
plot(b2, add = TRUE)
title("Sample locations")
par(mar = c(0, 0, 0, 0) + 0.1)
plot(e, axes = FALSE, cex = 0.6)


image(d2,col="green")
plot(b2,col="blue",add=TRUE,xlim="")
plot(e,add=TRUE,pch=1)
axis(1,at=seq(178100,181800,length=5),cex.axis = 0.7)
axis(2,at=seq(329000,333800,length=7),cex.axis = 0.7)
title(main="", sub="",
      xlab="X Coord", ylab="Y Coord")

box()


image(d2,col="lightgrey")
plot(b2, add = TRUE)
plot(e, pch = 1, cex = dist, add = TRUE)
plot(e, pch = 2, cex = lime/2, add = TRUE)
(max(dist)-min(dist))/4
legVals <- c(0, 0.22, 0.44, 0.66, 0.88)
legend("left", legend = c("0-0.22", "0.22-0.44", "0.44-0.66",
                          "0.66-0.88"), pch = 1, pt.cex = legVals,bty = "n", title = "measured")
axis(1,at=seq(178100,181800,length=5),cex.axis = 0.7)
axis(2,at=seq(329000,333800,length=7),cex.axis = 0.7)
box()


data(meuse)
# Continuous variables in meuse data:
# "cadmium" "copper"  "lead"    "zinc"    "elev"    "dist"    "om"   "dist.m" 

attach(meuse)

##
reg1<-glm(lime~cadmium ,family="binomial")
reg1
summary(reg1)
AIC(reg1)
BIC(reg1) 
reg2<-glm(lime~lead+zinc ,family="binomial")
reg2
summary(reg2)
AIC(reg2)
BIC(reg2) 
##


#install.packages("aod")
#install.packages("ggplot2")

library(aod)
library(ggplot2)
wald.test(coef(reg1),Sigma=vcov(reg1),Terms=1:2)
L1<-cbind(-1,-1)
wald.test(coef(reg1),Sigma=vcov(reg1),L=L1)

##
wald.test(coef(reg2),Sigma=vcov(reg2),Terms=1:3)
L2<-cbind(0,1,3)
wald.test(coef(reg2),Sigma=vcov(reg2),L=L2)
##


newdata2<-data.frame(lead=seq(min(lead),max(lead),length=1000),zinc=seq(min(zinc),max(zinc),length=1000)) 

newlimep2<-predict(reg2,newdata=newdata2,type="response")
names(newlimep2)<-NULL
newlime2<-rep(1,1000)
newlime2[newlimep2<0.5]<-0


ggplot(newdata2, aes(x = lead, y = newlimep2)) + geom_ribbon(aes(ymin = zinc,
                                                                 ymax = lead, fill = rank), alpha = 0.2) + geom_line(aes(colour = rank),
                                                                                                                     size = 1)

library(maps)
library(maptools)
library(lattice)
data(meuse)
d<-meuse


