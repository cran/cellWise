## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
  fig.width  = 5 ,
  fig.height = 3.5,
  fig.align  = 'center'
)

## -----------------------------------------------------------------------------
library("cellWise")
library("gridExtra") # has grid.arrange()

## ----fig.height=3,fig.width=7-------------------------------------------------
set.seed(12345) # for reproducibility
n <- 50; d <- 10
A <- matrix(0.9, d, d); diag(A) = 1 
round(A,1) # true covariance matrix
library(MASS) # only needed for the following line:
x <- mvrnorm(n, rep(0,d), A)
x[sample(1:(n * d), 50, FALSE)] <- NA
x[sample(1:(n * d), 25, FALSE)] <- 10
x[sample(1:(n * d), 25, FALSE)] <- -10
x <- cbind(1:n, x)
# When not specifying MacroPCApars all defaults are used:
MacroPCA.out <- MacroPCA(x, k = 0)

MacroPCA.out <- MacroPCA(x, k = 1)

## ----fig.height=9,fig.width=5-------------------------------------------------
cellMap(MacroPCA.out$stdResid)
# Red cells have higher value than predicted, blue cells lower,
# white cells are missing values, and all other cells are yellow.

## ----results='hide',message=FALSE,warning=FALSE-------------------------------
library(robustHD)
data(TopGear)
dim(TopGear)

## ----fig.height=3,fig.width=7-------------------------------------------------
myTopGear = TopGear
rownames(myTopGear) = paste(myTopGear[,1],myTopGear[,2]) 
# These rownames are make and model of the cars.
rownames(myTopGear)[165] = "Mercedes-Benz G" # name was too long
myTopGear = myTopGear[,-31] # removes subjective variable `Verdict'

# Transform some skewed variables:
transTG = myTopGear
transTG$Price = log(myTopGear$Price)
transTG$Displacement = log(myTopGear$Displacement)
transTG$BHP = log(myTopGear$BHP)
transTG$Torque = log(myTopGear$Torque)
transTG$TopSpeed = log(myTopGear$TopSpeed)

# Check the data:
checkData = checkDataSet(transTG, silent = TRUE)
# With option silent = FALSE we obtain more information:
# 
# The input data has 297 rows and 31 columns.
# 
# The input data contained 19 non-numeric columns (variables).
# Their column names are:
#   
#  [1] Maker         Model              Type               Fuel              
#  [5] DriveWheel    AdaptiveHeadlights AdjustableSteering AlarmSystem       
#  [9] Automatic     Bluetooth          ClimateControl     CruiseControl     
# [13] ElectricSeats Leather            ParkingSensors     PowerSteering     
# [17] SatNav        ESP                Origin            
# 
# These columns will be ignored in the analysis.
# We continue with the remaining 12 numeric columns:
#   
# [1] Price Cylinders Displacement BHP   Torque Acceleration TopSpeed    
# [8] MPG   Weight    Length       Width Height      
# 
# The data contained 1 rows with over 50% of NAs.
# Their row names are:
#   
#   [1] Citroen C5 Tourer
# 
# These rows will be ignored in the analysis.
# We continue with the remaining 296 rows:
#   
#   [1] Alfa Romeo Giulietta             Alfa Romeo MiTo                 
# .......
# [295] Volvo XC70                       Volvo XC90                      
# 
# The data contained 1 columns with zero or tiny median absolute deviation.
# Their column names are:
# 
# [1] Cylinders
# 
# These columns will be ignored in the analysis.
# We continue with the remaining 11 columns:
# 
# [1] Price  Displacement BHP   Torque Acceleration TopSpeed MPG         
# [8] Weight Length       Width Height      
# 
# The final data set we will analyze has 296 rows and 11 columns.

# The remainder of the dataset:
remTG = checkData$remX
dim(remTG) 

colSums(is.na(remTG)) # we still have quite a few NA's, especially in Weight:

# Check robust scale of the variables:
locscale = estLocScale(remTG)
round(locscale$scale,2)

# The scales are clearly different, so we will standardize before PCA.
# This is the argument scale=TRUE in MacroPCApars below.

# Small option lists (MacroPCA will automatically extend them with the
# default choices for the other options):
DDCpars <- list(fastDDC = FALSE, silent = TRUE)
MacroPCApars <- list(DDCpars = DDCpars, scale = TRUE, silent = TRUE)
# Note that MacroPCA needs DDCpars because it first runs DDC.

# To choose the number k of principal components we can run MacroPCA with k=0:
MacroPCAtransTG0 <- MacroPCA(transTG,k=0,MacroPCApars=MacroPCApars)

MacroPCAtransTG <- MacroPCA(transTG,k=2,MacroPCApars=MacroPCApars) # takes about 1 second.
names(MacroPCAtransTG)
MacroPCAtransTG$MacroPCApars # these are all the options used, starting with DDCpars:
MacroPCAtransTG$cutoffOD # cutoff value for orthogonal distances OD:
MacroPCAtransTG$cutoffSD # cutoff value for score distances OD:
length(MacroPCAtransTG$indcells) # list of flagged cells

## -----------------------------------------------------------------------------
ICPCAtransTG <- ICPCA(remTG,k=2,scale=TRUE,tolProb=0.99)
names(ICPCAtransTG)

# Compare imputed values for missings:
MacroPCAtransTG$X.NAimp[c(94,176),8]
ICPCAtransTG$X.NAimp[c(94,176),8] 
# CAR MAGAZINE: Ford Kuga Kuga 2.0 TDCi weighs 1605kg 
# CAR MAGAZINE: MINI Coupe 1.6T Cooper  weighs 1165kg

## ----fig.height=10,fig.width=6------------------------------------------------
# Make untransformed submatrix of X for labeling the cells 
# in the residual plot:
tempTG = myTopGear[checkData$rowInAnalysis,checkData$colInAnalysis]
dim(tempTG)
tempTG$Price = tempTG$Price/1000 # to avoid printing long numbers

showrows = c(12,42,50,51,52,59,72,94,98,135,150,164,176,
             180,195,196,198,209,210,215,219,234,259,277) # these 24 cars will be shown

# Make the ggplot2 objects for the residual maps by the function cellMap:
ggpICPCA = cellMap(ICPCAtransTG$stdResid, showcellvalues="D", 
                   D=tempTG, mTitle="ICPCA residual map", 
                   showrows=showrows, sizecellvalues = 0.7)

plot(ggpICPCA)

ggpMacroPCA = cellMap(MacroPCAtransTG$stdResid, showcellvalues="D",
                      D=tempTG, mTitle="MacroPCA residual map",
                       showrows=showrows, sizecellvalues = 0.7)
plot(ggpMacroPCA)

# Creating the combined pdf:
# pdf(file="TopGear_IPCA_MacroPCA_residualMap.pdf", width=12, height=10)
# gridExtra::grid.arrange(ggpICPCA, ggpMacroPCA, ncol=2) # arranges two plots on a page
# dev.copy(pdf, file="TopGear_IPCA_MacroPCA_residualMap.pdf", width=20, height=16)
# dev.off()

## ----fig.height=5,fig.width=5-------------------------------------------------
### Creating the outlier maps

outlierMap(MacroPCAtransTG,title="MacroPCA outlier map",
           labelOut=FALSE)
rowlabels = rownames(tempTG)
plotLabs = rep("",nrow(tempTG))
plotLabs[42]  = rowlabels[42] # BMW i3
plotLabs[50]  = "Bugatti" 
plotLabs[52]  = rowlabels[52] # Caterham Super 7
plotLabs[59]  = rowlabels[59] # Chevrolet Volt
plotLabs[180] = rowlabels[180] # Mitsubishi i-MiEV
plotLabs[195] = rowlabels[195] # Noble M600
plotLabs[196] = rowlabels[196] # Pagani Huayra
plotLabs[219] = rowlabels[219] # Renault Twizy
plotLabs[234] = rowlabels[234] # Ssangyong Rodius
plotLabs[259] = rowlabels[259] # Vauxhall Ampera
textPos = cbind(MacroPCAtransTG$SD,MacroPCAtransTG$OD)
textPos[42,1]  = textPos[42,1] -0.5 # BMW i3
textPos[42,2]  = textPos[42,2] +0.8 # BMW i3
textPos[50,1]  = textPos[50,1] -0.03 # Bugatti Veyron
textPos[50,2]  = textPos[50,2] +0.05 # Bugatti Veyron
textPos[52,1]  = textPos[52,1] +0.3  # Caterham Super 7
textPos[52,2]  = textPos[52,2] -0.5  # Caterham Super 7
textPos[59,1]  = textPos[59,1] -0.05 # Chevrolet Volt
textPos[59,2]  = textPos[59,2] -0.1  # Chevrolet Volt
textPos[180,1] = textPos[180,1] -1.2 # Mitsubishi i-MiEV
textPos[180,2] = textPos[180,2] +0.6 # Mitsubishi i-MiEV
textPos[195,1] = textPos[195,1] -0.05 # Noble M600
textPos[195,2] = textPos[195,2] +0.35 # Noble M600
textPos[196,1] = textPos[196,1] -0.6 # Pagani Huayra
textPos[196,2] = textPos[196,2] +0.7 # Pagani Huayra
textPos[219,1] = textPos[219,1] -0.03 # Renault Twizy
textPos[234,1] = textPos[234,1] +0.1 # Ssangyong Rodius
textPos[234,2] = textPos[234,2] +0.6 # Ssangyong Rodius
textPos[259,1] = textPos[259,1] -1.35 # Vauxhall Ampera
textPos[259,2] = textPos[259,2] +0.65 # Vauxhall Ampera
text(textPos,plotLabs,cex=0.8,pos=4)
# dev.copy(pdf,"TopGear_MacroPCA_outlierMap.pdf",width=6,height=6)
# dev.off()

outlierMap(ICPCAtransTG,title="ICPCA outlier map",labelOut=FALSE)
plotLabs = rep("",nrow(tempTG))
plotLabs[42]  = rowlabels[42]  # BMW i3
plotLabs[50]  = rowlabels[50]  # Bugatti Veyron
plotLabs[52]  = rowlabels[52]  # Caterham Super 7
plotLabs[59]  = rowlabels[59]  # Chevrolet Volt
plotLabs[180] = rowlabels[180] # Mitsubishi i-MiEV
plotLabs[195] = rowlabels[195] # Noble M600
plotLabs[196] = rowlabels[196] # Pagani Huayra
plotLabs[219] = rowlabels[219] # Renault Twizy
plotLabs[234] = rowlabels[234] # Ssangyong Rodius
plotLabs[259] = rowlabels[259] # Vauxhall Ampera
textPos = cbind(ICPCAtransTG$SD,ICPCAtransTG$OD)
textPos[50,1]  = textPos[50,1] -0.05 # Bugatti Veyron
textPos[50,2]  = textPos[50,2] +0.05 # Bugatti Veyron
textPos[52,1]  = textPos[52,1] -0.07 # Caterham Super 7
textPos[52,2]  = textPos[52,2] -0.37 # Caterham Super 7
textPos[59,1]  = textPos[59,1] -0.05 # Chevrolet Volt
textPos[59,2]  = textPos[59,2] -0.3  # Chevrolet Volt
textPos[180,1] = textPos[180,1] -1.2 # Mitsubishi i-MiEV
textPos[180,2] = textPos[180,2] +0.23 # Mitsubishi i-MiEV
textPos[195,1] = textPos[195,1] +0.25 # Noble M600
textPos[195,2] = textPos[195,2] +0.2 # Noble M600
textPos[196,1] = textPos[196,1] -0.05 # Pagani Huayra
textPos[196,2] = textPos[196,2] +0.32 # Pagani Huayra
textPos[219,1] = textPos[219,1] -0.8 # Renault Twizy
textPos[219,2] = textPos[219,2] +0.4 # Renault Twizy
textPos[234,1] = textPos[234,1] -0.05 # Ssangyong Rodius
textPos[234,2] = textPos[234,2] +0.2 # Ssangyong Rodius
textPos[259,1] = textPos[259,1] -0.9 # Vauxhall Ampera
textPos[259,2] = textPos[259,2] +0.4 # Vauxhall Ampera
text(textPos,plotLabs, cex=0.8, pos=4)
# dev.copy(pdf,"TopGear_ICPCA_outlierMap.pdf",width=6,height=6)
# dev.off()

## ----fig.height=10,fig.width=6------------------------------------------------

# For comparison, remake the residual map of the entire dataset, but now 
# showing the values of the residuals instead of the data values:

ggpMacroPCAres = cellMap(MacroPCAtransTG$stdResid, 
                         showcellvalues="R", sizecellvalues = 0.7,
                         mTitle="MacroPCA residual map",
                         showrows=showrows) 
plot(ggpMacroPCAres)

# Define the "initial" dataset as the rows not in these 24:
initX = remTG[-showrows,]
dim(initX) 

# Fit initX:
MacroPCAinitX = MacroPCA(initX,k=2,MacroPCApars=MacroPCApars) 

# Define the "new" data to predict:
newX = remTG[showrows,]
dim(newX) 

# Make predictions by MacroPCApredict. 
# Its inputs are:
#
# Xnew            : the new data (test data), which must be a 
#                   matrix or a data frame. 
#                   It must always be provided.
# InitialMacroPCA : the output of the MacroPCA function on the 
#                   initial (training) dataset. Must be provided.
# MacroPCApars    : the input options to be used for the prediction.
#                   By default the options of InitialMacroPCA
#                   are used. For the complete list of options
#                   see the function MacroPCA.

predictMacroPCA = MacroPCApredict(newX,MacroPCAinitX)
# We did not need to specify the third argument because it is taken
# from the initial fit MacroPCAinitX .

names(predictMacroPCA)
# The outputs are similar to those of of MacroPCA.

# Make the residual map:
ggpMacroPCApredict = cellMap(predictMacroPCA$stdResid, 
                             showcellvalues="R", sizecellvalues = 0.7,
                             mTitle="MacroPCApredict residual map")
plot(ggpMacroPCApredict) # is very similar to that based on all the data!

# Creating the combined pdf:
# pdf(file="TopGear_MacroPCApredict_residualMap.pdf",width=12,height=10)
# gridExtra::grid.arrange(ggpMacroPCAres,ggpMacroPCApredict,ncol=2)
# dev.off()

## -----------------------------------------------------------------------------
data(data_glass)

## ----results='hide',message=FALSE,warning=FALSE-------------------------------
library(rrcov) # for robust PCA:

## ----fig.height=6,fig.width=7-------------------------------------------------
dim(data_glass) 

# Do not scale the spectra in the glass data:
MacroPCApars$scale = FALSE

# Check data
checkData = checkDataSet(data_glass, silent=TRUE) 
# With checkData = checkDataSet(glass, silent=FALSE) we obtain more information:
#
#  The input data has 180 rows and 750 columns.
#
#  The data contained 11 discrete columns with 3 or fewer values.
#  Their column names are:
#
#  [1] V1  V2  V3  V4  V5  V6  V7  V8  V9  V10 V11
#
#  These columns will be ignored in the analysis.
#  We continue with the remaining 739 columns:
#  
#   [1] V12  V13  V14  V15  V16  V17  V18  V19  V20  V21  V22  V23  V24  V25
#    ......
# [729] V740 V741 V742 V743 V744 V745 V746 V747 V748 V749 V750
#
#  The data contained 2 columns with zero or tiny median absolute deviation.
#  Their column names are:
#
# [1] V12 V13
#
#  These columns will be ignored in the analysis.
#  We continue with the remaining 737 columns:
#
#   [1] V14  V15  V16  V17  V18  V19  V20  V21  V22  V23  V24  V25  V26  V27
#    ......
# [729] V742 V743 V744 V745 V746 V747 V748 V749 V750
#
#  The final data set we will analyze has 180 rows and 737 columns.

remglass = checkData$remX
n <- nrow(remglass); n 
d = ncol(remglass); d # the first 13 variables had scale 0:

# Compare ICPCA and MacroPCA:

ICPCAglass <- ICPCA(remglass,k=4,scale=F,tolProb=0.99)

MacroPCAglass = MacroPCA(remglass,k=4,MacroPCApars=MacroPCApars) # takes 8 seconds

nrowsinblock = 5
rowlabels = rep("",floor(n/nrowsinblock));
rowlabels[1] = "1"
rowlabels[floor(n/nrowsinblock)] = "n";

ncolumnsinblock = 5
columnlabels = rep("",floor(d/ncolumnsinblock));
columnlabels[1] = "1";
columnlabels[floor(d/ncolumnsinblock)] = "d"

ggpICPCA <- cellMap(ICPCAglass$stdResid,
                    rowblocklabels=rowlabels,    
                    columnblocklabels=columnlabels,
                    mTitle="ICPCA residual map",
                    rowtitle="glass samples",
                    columntitle="wavelengths",
                    nrowsinblock=5,
                    ncolumnsinblock=5)

ggpMacroPCA <- cellMap(MacroPCAglass$stdResid,
                       rowblocklabels=rowlabels,
                       columnblocklabels=columnlabels,
                       mTitle="MacroPCA residual map",
                       rowtitle="glass samples",
                       columntitle="wavelengths",
                       nrowsinblock=5,
                       ncolumnsinblock=5)

grid.arrange(ggpMacroPCA,ggpICPCA,nrow=2) 
# pdfName = "Glass_MacroPCA_ICPCA_residualMap.pdf"
# dev.copy(pdf, pdfName, width=12, height=8)
# dev.off()


## ----fig.height=6,fig.width=7-------------------------------------------------

library(rrcov) # only needed for PcaHubert()
ROBPCAglass = PcaHubert(remglass,k=4,alpha=0.5)
# Calculate ROBPCA residuals and standardize them:
Xhat = sweep(ROBPCAglass@scores %*% t(ROBPCAglass@loadings),
             2,ROBPCAglass@center,"+")
Xresid = remglass - Xhat
scaleRes = estLocScale(Xresid,type="1stepM",center=F)$scale
stdResidROBPCA = sweep(Xresid,2,scaleRes,"/")

ggpROBPCA <- cellMap(stdResidROBPCA,
                     rowblocklabels=rowlabels,
                     columnblocklabels=columnlabels,
                     mTitle="ROBPCA residual map",
                     rowtitle="glass samples",
                     columntitle="wavelengths",
                     nrowsinblock=5,
                     ncolumnsinblock=5)

grid.arrange(ggpMacroPCA, ggpROBPCA, nrow=2) 
# pdfName = "Glass_MacroPCA_ROBPCA_residualMap.pdf"
# dev.copy(pdf, pdfName, width=12, height=8)
# dev.off()


## ----fig.height=6,fig.width=7-------------------------------------------------

fastDDCpars = list(fastDDC=TRUE, silent=TRUE)
fastMacroPCApars = list(DDCpars=fastDDCpars, scale=FALSE, silent=TRUE)

fastMacroPCAglass = MacroPCA(data_glass,k=4,MacroPCApars=fastMacroPCApars) # 2 seconds

ggpfastMacroPCA <- cellMap(fastMacroPCAglass$stdResid,
                           columnblocklabels=columnlabels,
                           rowblocklabels=rowlabels,
                           mTitle="MacroPCA with fastDDC=T",
                           columntitle="wavelengths",
                           rowtitle="glass samples",
                           ncolumnsinblock=5,
                           nrowsinblock=5)

grid.arrange(ggpMacroPCA, ggpfastMacroPCA, nrow=2) # The results are similar:
# pdfName = "Glass_MacroPCA_residualMaps.pdf"
# dev.copy(pdf, pdfName, width=12, height=8)
# dev.off()

## -----------------------------------------------------------------------------
data(data_dposs)

## ----fig.height=3,fig.width=7-------------------------------------------------

colnames(data_dposs); dim(data_dposs) 
n = nrow(data_dposs); n

# Count missing cells
missmat = is.na(data_dposs)
sizemat = nrow(missmat)*ncol(missmat); sizemat 
100*sum(as.vector(missmat))/sizemat # 50.2% of the values are missing:

# Count rows with missings
missrow = length(which(rowSums(missmat) > 0))
100*missrow/nrow(missmat) # 84.6% of the rows contain missing values:

# PERFORM ICPCA AND MACROPCA

ICPCAdposs = ICPCA(data_dposs,k=4,scale=TRUE)
names(ICPCAdposs)

# MacroPCA options with fracNA allowing for many NA's:
DDCPars = list(fastDDC=F,fracNA=1.0)
MacroPCAPars = list(DDCpars=DDCPars,scale=TRUE,silent=T)
MacroPCAdposs = MacroPCA(data_dposs,k=4,MacroPCApars=MacroPCAPars) # takes 6 seconds

## SCREE PLOTS

barplot(ICPCAdposs$eigenvalues,
        main="ICPCA scree plot", ylab="eigenvalues",
        names.arg=1:length(ICPCAdposs$eigenvalues))
# dev.copy(pdf,"DPOSS_ICPCA_screeplot.pdf",width=6,height=6)
# dev.off()

barplot(MacroPCAdposs$eigenvalues,
        main="MacroPCA scree plot", ylab="eigenvalues",
        names.arg=1:length(MacroPCAdposs$eigenvalues))
# Not as concentrated in the first eigenvalue.
# dev.copy(pdf,"DPOSS_MacroPCA_screeplot.pdf",width=6,height=6)
# dev.off()

## ----fig.height=4,fig.width=4-------------------------------------------------
## LOADINGS
ICPCAdposs$loadings[,2] = -ICPCAdposs$loadings[,2]

matplot(ICPCAdposs$loadings[,1:2],main="ICPCA loadings",
        xlab="variables",ylab="Loadings",col=c("black","blue"),
        ylim=c(-0.4,0.6),type="l",lty=c(1,2),lwd=2)
abline(v=7.5,col="red")
abline(v=14.5,col="red")
# dev.copy(pdf,"DPOSS_ICPCA_loadings.pdf",width=5,height=5)
# dev.off()

matplot(MacroPCAdposs$loadings[,1:2],main="MacroPCA loadings",
        xlab="variables",ylab="Loadings",col=c("black","blue"),
        ylim=c(-0.3,0.5),type="l",lty=c(1,2),lwd=2)
abline(v=7.5,col="red")
abline(v=14.5,col="red")
# dev.copy(pdf,"DPOSS_MacroPCA_loadings.pdf",width=5,height=5)
# dev.off()

## ----fig.height=4,fig.width=4-------------------------------------------------
# Select the 150 rows with highest OD from MacroPCA
dpossOD = MacroPCAdposs$OD
summary(dpossOD) 
n = length(dpossOD); n
sortOD = order(dpossOD,1:n)
summary(sortOD); sortOD[1:5]
indHighOD = sortOD[n-150+1:n]
indHighOD = na.omit(indHighOD)
summary(indHighOD); length(indHighOD) 
indLowOD = sortOD[1:(n-150)]
indLowOD = na.omit(indLowOD)
summary(indLowOD); length(indLowOD)

## OUTLIER MAPS

dpossColList = list(class1=list(col="black",index=indLowOD),
                    class2=list(col="red",index=indHighOD))

outlierMap(ICPCAdposs,title="ICPCA outlier map", 
           col=dpossColList,labelOut=FALSE)
# dev.copy(pdf,"DPOSS_ICPCA_outlierMap.pdf",width=8,height=8)
# dev.off()

outlierMap(MacroPCAdposs,title="MacroPCA outlier map", 
           col=dpossColList,labelOut=FALSE)
# dev.copy(pdf,"DPOSS_MacroPCA_outlierMap.pdf",width=8,height=8)
# dev.off()



## ----fig.height=4,fig.width=4-------------------------------------------------
## ICPCA SCORE PLOTS
ICPCAdposs$scores[,2] = -ICPCAdposs$scores[,2]

plot(ICPCAdposs$scores[,1:2],main="ICPCA scores",xlab="PC1",ylab="PC2")
points(ICPCAdposs$scores[indHighOD,1:2],pch=16,col="red")
# dev.copy(pdf,"DPOSS_ICPCA_Scores12.pdf",width=5,height=5)
# dev.off()

plot(ICPCAdposs$scores[,c(1,3)],main="ICPCA scores",xlab="PC1",ylab="PC3")
points(ICPCAdposs$scores[indHighOD,c(1,3)],pch=16,col="red")

plot(ICPCAdposs$scores[,2:3],main="ICPCA scores",xlab="PC2",ylab="PC3")
points(ICPCAdposs$scores[indHighOD,2:3],pch=16,col="red")

# MacroPCA SCORE PLOTS
MacroPCAdposs$scores[,2] = -MacroPCAdposs$scores[,2]
MacroPCAdposs$scores[,3] = -MacroPCAdposs$scores[,3]

plot(MacroPCAdposs$scores[,1:2],main="MacroPCA scores",xlab="PC1",ylab="PC2")
points(MacroPCAdposs$scores[indHighOD,1:2],pch=16,col="red")
# dev.copy(pdf,"DPOSS_MacroPCA_Scores12.pdf",width=5,height=5)
# dev.off()

plot(MacroPCAdposs$scores[,c(1,3)],main="MacroPCA scores",xlab="PC1",ylab="PC3")
points(MacroPCAdposs$scores[indHighOD,c(1,3)],pch=16,col="red")

plot(MacroPCAdposs$scores[,2:3],main="MacroPCA scores",xlab="PC2",ylab="PC3")
points(MacroPCAdposs$scores[indHighOD,2:3],pch=16,col="red")


## ----fig.height=5,fig.width=5-------------------------------------------------
# RESIDUAL MAPS

dpossH = na.omit(indHighOD)
summary(dpossH); length(as.vector(dpossH)) 
dpossL = na.omit(indLowOD)
summary(dpossL); length(as.vector(dpossL)) 
set.seed(0)
dpossH = sample(dpossH,150)
dpossL = sample(dpossL,300)
showrowsdposs = c(dpossH,dpossL)
summary(showrowsdposs)
length(showrowsdposs) 

rowlabels = c("OS1","OS2","OS3","OS4","OS5","OS6",
              "S1","S2","S3","S4","S5","S6","S7","S8",
              "S9","S10","S11","S12")


ggpICPCAdposs = cellMap(ICPCAdposs$stdResid,
                        rowblocklabels=rowlabels,
                        mTitle="ICPCA residual map",
                        rowtitle="",
                        showrows=showrowsdposs,
                        nrowsinblock=25,
                        ncolumnsinblock=1,
                        sizetitles=1.5)
plot(ggpICPCAdposs) # not much to see:
# dev.copy(pdf,"DPOSS_ICPCA_residualMap.pdf",width=8,height=6)
# dev.off()


ggpMacroPCAdposs = cellMap(MacroPCAdposs$stdResid, 
                           rowblocklabels=rowlabels,
                           mTitle="MacroPCA residual map",
                           rowtitle="",
                           showrows=showrowsdposs,
                           nrowsinblock=25,
                           ncolumnsinblock=1,
                           sizetitles=1.5)
plot(ggpMacroPCAdposs) # interesting structure:
# dev.copy(pdf,"DPOSS_MacroPCA_residualMap.pdf",width=8,height=6)
# dev.off()


