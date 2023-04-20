## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
 fig.width = 8 ,
 fig.height = 12,
 fig.align ='center'
)

## -----------------------------------------------------------------------------
library("cellWise")

## -----------------------------------------------------------------------------
library(robustHD)
data(TopGear)
dim(TopGear)
rownames(TopGear) = paste(TopGear[,1],TopGear[,2])
colnames(TopGear)
TopGear = TopGear[c(5,7,9:17)] # objective numerical variables
colnames(TopGear)
colnames(TopGear) = c("price", "displacement", "horsepower",
                      "torque", "acceleration", "top speed",
                      "miles/gallon", "weight", "length",
                      "width", "height")

## -----------------------------------------------------------------------------
out = checkDataSet(TopGear, silent = TRUE)
Xorig = out$remX
dim(Xorig) 
Xorig[,1] = Xorig[,1]/1000
head(Xorig)

## -----------------------------------------------------------------------------
X = Xorig
X[,c(1:4,6)] = log(Xorig[,c(1:4,6)])

## -----------------------------------------------------------------------------
X = as.matrix(X)
colSums(is.na(X))

## -----------------------------------------------------------------------------
out = estLocScale(X)
Z = scale(X,center=out$loc, scale=out$scale)
cutf = sqrt(qchisq(p=0.99,df=1))
boxplot(Z)

## -----------------------------------------------------------------------------
out <- robustbase::covMcd(X) 
plot(sqrt(out$mah))
abline(h = sqrt(qchisq(0.99, df = 11)))

## -----------------------------------------------------------------------------
cellout <- cellMCD(X)

## -----------------------------------------------------------------------------
Zres = cellout$Zres
boxplot(Zres) 
qqnorm(as.vector(Zres)); abline(0,1) 

## -----------------------------------------------------------------------------
W = cellout$W
sum(as.vector(1-W)) 
rowS = rowSums(1-W)
sum(rowS == 0) # 151 cars do not have a single flagged cell.

## -----------------------------------------------------------------------------
j = 1 # price
# Index plot of Zres:
# (ids = plot.cellMCD(cellout, j, type="index", identify=T)$ids)
ids = c(6,54,50,195,212,222,221)
Zres[ids,j,drop=F]
plot_cellMCD(cellout, j, type="index", ids=ids)

## -----------------------------------------------------------------------------
# Plot Zres versus X:
plot_cellMCD(cellout, j, type="Zres/X", ids=ids)
# Plot Zres versus predicted cell:
plot_cellMCD(cellout, j, type="Zres/pred", ids=ids)
# Plot observed cell versus predicted cell:
# (ids = plot.cellMCD(cellout, j, type="X/pred", identify=T)$ids)
ids = c(179,3,212,218,54,6,222,221,50,195)
Xorig[ids,j,drop=F]
plot_cellMCD(cellout, j, type="X/pred", vband=T, hband=T, ids=ids)

## ----fig.height=16,fig.width=8------------------------------------------------
###########
## Figure 2
###########

# pdf(file="TopGear_fig2_test.pdf",width=5.2,height=10.4)
par(mfrow=c(2,1)) # (nrows,ncols)
rn = rownames(X) # for object labels
#
######### 1. HP
par(mar = c(3.6, 3.5, 2.8, 1.0))
j = 3
ids = c(52,59,70,218)
labs = rn[ids]; labs[1] = "Caterham"
xvals = c(52, 59, 70, 218) 
yvals = c(-5.70, -7.70, -4.20, -4.35)
plot_cellMCD(cellout, j, type="index") 
text(x = xvals, y = yvals, labels = labs, cex = 0.9)
#
######### 2. Length
par(mar = c(3.6, 3.5, 2.8, 1.0))
#
j = 9
ids = c(3,50,51,119,144,195,218,232,249)
labs = rn[ids]
labs[1] = "Aston Cygnet"; labs[3] = "Caterham"
xvals = c(2560, 5050, 3100, 5252, 4245,
          4605, 2360, 2700, 3330) 
yvals = c(-4.11, -3.88, -3.60, 3.87, -4.30,
          -5.00, -5.00, -5.90, -4.85)
plot_cellMCD(cellout, j, type = "Zres/X",
             main = "standardized residual versus X for length", xlab = "length (mm)", xlim = c(1970, 6000),  ylim = c(-6, 4.5)) #,ids=ids)
text(x = xvals, y = yvals, labels = labs, cex = 0.9)
# dev.off()

## ----fig.height=16,fig.width=8------------------------------------------------

# pdf(file="TopGear_fig3_test.pdf",width=5.2,height=10.4)
par(mfrow=c(2,1))
rn = rownames(X) # for object labels
#
######### 1. weight
par(mar = c(3.6, 3.5, 2.8, 1.0))
j = 8
ids = c(29,51,144,163,183,197,195)
pred = cellout$preds
labs = rn[ids]
labs[2] = "Caterham"
labs[4] = "Mercedes-Benz G"
xvals = c(1490, 1000, 1890, 1575, 1757,  810, 2250) 
yvals = c(5.22,-4.70,-4.73, 4.19,-5.65,-6.05,-7.60)
plot_cellMCD(cellout, j, type="Zres/pred", # vband=F,
             xlab = "predicted weight (kg)") # , ids=ids)
text(x = xvals, y = yvals, labels = labs, cex = 0.9)
#
######### 2. top speed
par(mar = c(3.6, 3.5, 2.8, 1.0))
j = 6
ids = c(50,42,52,79,195,218,219,232,258)
labs = rn[ids]
labs[3] = "Caterham"
xvals = c(5.145, 4.937, 5.083, 5.065, 5.076, 4.90, 4.75, 4.63, 5.038)
yvals = c(5.53, 4.53, 4.72, 5.345, 5.44, 3.97, 4.43, 4.30, 4.61)
plot_cellMCD(cellout, j, type="X/pred", vband=F, vlines=F,
             xlab="predicted top speed (transformed)")#, ids=ids)
text(x = xvals, y = yvals, labels = labs, cex = 0.9)
# dev.off()

## ----fig.height=16,fig.width=8------------------------------------------------

############
### Figure 4
############

# pdf(file="TopGear_fig4_test.pdf",width=5.2,height=10.4)
par(mfrow=c(2,1)) # (nrows,ncols)
rn = rownames(X) # for object labels
#
######### 1. MPG versus torque
par(mar = c(3.6, 3.5, 2.8, 1.0))
ids=c(42,258,50) 
labs = rn[ids]
xvals = c(5.56, 6.28, 7.27) 
yvals = c(470, 235, -16)
plot_cellMCD(cellout, type = "bivariate", horizvar = 4, 
             vertivar = 7, ids = ids, 
             xlim=c(3.5,7.7), ylim = c(-27,480), 
             main = "miles/gallon versus torque",
             xlab = "torque (transformed)",
             ylab = "miles/gallon",
             opacity=0.5, labelpoints=F)
text(x = xvals, y = yvals, labels = labs, cex = 0.9)
#
######### 2. Width versus acceleration
par(mar = c(3.6, 3.5, 2.8, 1.0))
ids = c(218,144,233,51,136,179) 
labs = rn[ids]
labs[3] = "Ssangyong" # name was too long for plot 
labs[4] = "Caterham" 
labs[5] = "Land Rover"
xvals = c(0.0,  -3.0, -3.23,  0.5, 14.2, 15.9) 
yvals = c(1200, 1850, 1915, 1688, 2039, 1445)
plot_cellMCD(cellout,  type = "bivariate", horizvar=5, 
             vertivar=10, ids=ids,
             xlab = "acceleration (seconds to 62 mph)",
             ylab = "width (mm)", xlim = c(-6.1,20),
             ylim = c(1200, 2150),
             opacity=0.5, labelpoints=F)
text(x = xvals, y = yvals, labels = labs, cex = 0.9)
# dev.off()

par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
library(robustbase)

## -----------------------------------------------------------------------------
data(aircraft)
X      <- as.matrix(aircraft)[, 1:4]
pairs(X) 
X      <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd     <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975,ncol(X)))
sum(rd > cutoff)
cellout <- cellMCD(X)

cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
plot(rd, cd)
round(cor(rd, cd),3)

## -----------------------------------------------------------------------------
data(aircraft)
X <- as.matrix(aircraft)
nrow(X) / ncol(X)

## ----error=TRUE---------------------------------------------------------------
cellout <- cellMCD(X)

## -----------------------------------------------------------------------------
y <- X[, 5]
y <- transfo(y)$Y

## -----------------------------------------------------------------------------
data(alcohol, package = "robustbase")
X <- as.matrix(alcohol)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
eigen(rowout$cov)$values
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff)
cellout <- cellMCD(X)

cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd),3)  

## -----------------------------------------------------------------------------
data(Animals2)
X <- as.matrix(log(Animals2))
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
plot(rd, cd)
round(cor(rd, cd), 3) 

## -----------------------------------------------------------------------------
data(bushfire)
X <- as.matrix(bushfire)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X,center=rowout$center,cov=rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff)
cellout <- cellMCD(X) 
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3)

## -----------------------------------------------------------------------------
data(cloud)
X <- as.matrix(cloud)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cellout <- cellMCD(X) 

cd <- sqrt(mahalanobis(X,center=cellout$mu,cov=cellout$S))
round(cor(rd,cd),3) 

## -----------------------------------------------------------------------------
data(delivery)
X <- as.matrix(delivery)[, 1:2]
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff)
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, -cd), 3) 

## -----------------------------------------------------------------------------
data(delivery)
X <- as.matrix(delivery)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff) 
cellout <- cellMCD(X) 
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3) 

## -----------------------------------------------------------------------------
data(exAM)
X <- as.matrix(exAM)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff) 
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3)

## -----------------------------------------------------------------------------
data(hbk)
X <- as.matrix(hbk)[, 1:3]
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff) 
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3) 

## -----------------------------------------------------------------------------
data(hbk)
X <- as.matrix(hbk)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff)
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3)

## -----------------------------------------------------------------------------
data(kootenay)
X <- as.matrix(kootenay)
X <- transfo(X)$Y 
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff)
cellout <- cellMCD(X)
cd = sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3) 

## ---- error=TRUE--------------------------------------------------------------
data(lactic)
X <- as.matrix(lactic)
X <- transfo(X)$Y

## -----------------------------------------------------------------------------
data(milk)
X <- as.matrix(milk)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3) 

## ---- error=TRUE--------------------------------------------------------------
data(pension)
X <- as.matrix(pension)
X <- transfo(X)$Y

## -----------------------------------------------------------------------------
data(phosphor)
X <- as.matrix(phosphor)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3) 

## -----------------------------------------------------------------------------
data(pilot)
X <- as.matrix(pilot)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff) 
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3) 

## -----------------------------------------------------------------------------
cellout$mu - colMeans(X) 
n <- nrow(X)
cellout$S - (n - 1) * cov(X) / n

## -----------------------------------------------------------------------------
data(radarImage)
X <- as.matrix(radarImage)
X <- transfo(X)$Y
pairs(X) 
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975,ncol(X)))
sum(rd > cutoff) 
cellout <- cellMCD(X) 
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd,cd),3)

## -----------------------------------------------------------------------------
data(salinity)
X <- as.matrix(salinity)[, 1:3]
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff) 
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd,cd),3) 

## -----------------------------------------------------------------------------
data(salinity)
X <- as.matrix(salinity)
X <- transfo(X)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff) 
cellout <- cellMCD(X) 
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd, cd), 3) 

## -----------------------------------------------------------------------------
data(starsCYG)
X <- as.matrix(starsCYG)
X <- transfo(X, nbsteps = 1)$Y
rowout <- robustbase::covMcd(X)
rd <- sqrt(mahalanobis(X, center = rowout$center, cov = rowout$cov))
cutoff <- sqrt(qchisq(0.975, ncol(X)))
sum(rd > cutoff)
cellout <- cellMCD(X)
cd <- sqrt(mahalanobis(X, center = cellout$mu, cov = cellout$S))
round(cor(rd,cd),3) 

## ---- error = TRUE------------------------------------------------------------
data(telef)
X <- as.matrix(telef)
# plot(X)
X <- transfo(X)$Y

## ---- error = TRUE------------------------------------------------------------
data(toxicity)
X <- as.matrix(toxicity)
X <- transfo(X, nbsteps = 1)$Y

## -----------------------------------------------------------------------------
data(wood)
X <- as.matrix(wood)
X <- transfo(X, nbsteps = 1)$Y
cellout <- cellMCD(X)

