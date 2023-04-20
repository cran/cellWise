## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  fig.width  = 5 ,
  fig.height = 3.5,
  fig.align  = 'center'
)

## -----------------------------------------------------------------------------
library("cellWise")
library("gridExtra") # has grid.arrange()

# Default options for DDC:
DDCpars = list(fracNA = 0.5, numDiscrete = 3, precScale = 1e-12,
               cleanNAfirst = "automatic", tolProb = 0.99, 
               corrlim = 0.5, combinRule = "wmean",
               returnBigXimp = FALSE, silent = FALSE,
               nLocScale = 25000, fastDDC = FALSE,
               standType = "1stepM", corrType = "gkwls",
               transFun = "wrap", nbngbrs = 100)

# A small list giving the same results:
DDCpars = list(fastDDC = FALSE)


## -----------------------------------------------------------------------------

i = c(1,2,3,4,5,6,7,8,9) 
name = c("aa","bb","cc","dd","ee","ff","gg","hh","ii") 
logic = c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) 
V1 = c(1.3,NaN,4.5,2.7,20.0,4.4,-2.1,1.1,-5)
V2 = c(2.3,NA,5,6,7,8,4,-10,0.5)
V3 = c(2,Inf,3,-4,5,6,7,-2,8)
Vna = c(1,-4,2,NaN,3,-Inf,NA,6,5)
Vdis = c(1,1,2,2,3,3,3,1,2)
V0s = c(1,1.5,2,2,2,2,2,3,2.5) 
datafr = data.frame(i,name,logic,V1,V2,V3,Vna,Vdis,V0s) 
datafr

DDCdatafr = DDC(datafr,DDCpars)

remX = DDCdatafr$remX; dim(remX)
cellMap(DDCdatafr$stdResid)
# Red cells have higher value than predicted, blue cells lower,
# white cells are missing values, all other cells are yellow.

## ----fig.height=6,fig.width=3-------------------------------------------------
set.seed(12345) # for reproducibility
n <- 50; d <- 20
A <- matrix(0.9, d, d); diag(A) = 1
round(A[1:10,1:10],1) # true covariance matrix
library(MASS) # only needed for the following line:
x <- mvrnorm(n, rep(0,d), A)
x[sample(1:(n * d), 50, FALSE)] <- NA
x[sample(1:(n * d), 50, FALSE)] <- 10
x[sample(1:(n * d), 50, FALSE)] <- -10

# When not specifying DDCpars in the call to DDC
# all defaults are used:
DDCx <- DDC(x)

cellMap(DDCx$stdResid)
# Red cells have higher value than predicted, blue cells lower,
# white cells are missing values, all other cells are yellow.

## -----------------------------------------------------------------------------
DDCx$DDCpars # These are the default options:

names(DDCx)

# We will now go through these outputs one by one:

DDCx$colInAnalysis # all columns X1,...,X20 remain:
DDCx$rowInAnalysis # all rows 1,...,50 remain:
DDCx$namesNotNumeric # no non-numeric columns:
DDCx$namesCaseNumber # no column was the case number:
DDCx$namesNAcol # no columns with too many NA's:
DDCx$namesNArow # no columns with too many NA's:
DDCx$namesDiscrete # no discrete columns:
DDCx$namesZeroScale # no columns with scale = 0:
dim(DDCx$remX) # remaining data matrix
round(DDCx$locX,2) # robust location estimates of all 20 columns:
round(DDCx$scaleX,2) # robust scale estimates of all 20 columns:

round(DDCx$Z[1:10,1:10],1)
# Robustly standardized dataset. Due to the high correlations,
# cells in the same row look similar (except for outlying cells).

DDCx$nbngbrs
# For each column the code looked for up to 19 non-self neighbors (highly correlated columns).
# It goes through all of them, unless fastDDC is set to TRUE.

DDCx$ngbrs[1:3,]
# Shows the neighbors, e.g. the nearest non-self neighbor of X1 is X11, then X2,...

round(DDCx$robcors[1:3,],2)
# Robust correlations with these neighbors. In each row the correlations
# are sorted by decreasing absolute value. 

round(DDCx$robslopes[1:3,],2)
# For each column, the slope of each neighbor predicting it.
# For instance, X1 is predicted by its first neighbor with 
# slope 0.97 and by its second neighbor with slope 0.81 .

round(DDCx$deshrinkage,2)
# For each column, the factor by which its prediction is multiplied.

## -----------------------------------------------------------------------------
round(DDCx$Xest[1:12,1:10],2) # the estimated cells of remX:

round(DDCx$stdResid[1:12,1:10],1)
# The standardized residuals of the cells. Note the NA's and some
# large positive and negative cell residuals.

qqnorm(as.vector(DDCx$stdResid)) # Note the far outliers on both sides:

as.vector(DDCx$indcells) # indices of the cells that were flagged by DDC:

plot(DDCx$Ti) # the Ti values of the rows. None are high.
qqnorm(DDCx$Ti) 
DDCx$medTi # median of the raw Ti (used to standardize Ti):
DDCx$madTi # median absolute deviation of the raw Ti (used to standardize Ti):

DDCx$indrows # numeric(0) means no rows are flagged:

as.vector(DDCx$indall) # indices of the flagged cells, including those in flagged rows:

as.vector(DDCx$indNAs) # indices of the missing cells:

round(DDCx$Ximp[1:10,1:10],2)
# The imputed matrix. Both the cellwise outliers and the missing values
# are replaced by their predicted values.

round((DDCx$Ximp - DDCx$remX)[1:10,1:10],2)
# The nonzero values and the NA's correspond to imputed cells.


## ----results='hide',message=FALSE,warning=FALSE-------------------------------
library(robustHD)
data(TopGear)

## -----------------------------------------------------------------------------
dim(TopGear)
rownames(TopGear)[1:13] # "1" to "297" are not useful names
rownames(TopGear) = paste(TopGear[,1],TopGear[,2]) 
# Now the rownames are make and model of the cars.
rownames(TopGear)[165] = "Mercedes-Benz G" # name was too long
myTopGear = TopGear[,-31] # removes the subjective variable `Verdict'

# Transform some variables to get roughly gaussianity in the center:
transTG = myTopGear
transTG$Price = log(myTopGear$Price)
transTG$Displacement = log(myTopGear$Displacement)
transTG$BHP = log(myTopGear$BHP)
transTG$Torque = log(myTopGear$Torque)
transTG$TopSpeed = log(myTopGear$TopSpeed)

# Run the DDC method:
DDCpars = list(fastDDC = FALSE, silent = TRUE)
DDCtransTG = DDC(transTG,DDCpars)
# With DDCpars = list(fastDDC = FALSE, silent = FALSE) we obtain more information:
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



## ----fig.height=8,fig.width=6-------------------------------------------------
remX = DDCtransTG$remX # the remaining part of the dataset
dim(remX)

colSums(is.na(remX)) # There are still NAs, mainly in `Weight':

# Analyze the data by column:
standX = scale(remX,apply(remX,2,median,na.rm = TRUE), 
               apply(remX,2,mad,na.rm = TRUE))
dim(standX)
round(standX[1:5,],1) # has NAs where remX does
transTGcol = remX
transTGcol[abs(standX) > sqrt(qchisq(0.99,1))] = NA
round(transTGcol[1:5,],1) # has NAs in outlying cells as well:

# Make untransformed submatrix of X for labeling the cells in the plot:
tempX = myTopGear[DDCtransTG$rowInAnalysis,DDCtransTG$colInAnalysis]
tempX$Price = tempX$Price/1000 # to avoid printing long numbers
dim(tempX)

# Show the following 17 cars in the cellmap:
showrows = c(12,42,56,73,81,94,99,135,150,164,176,198,209,215,234,241,277)

# Make two ggplot2 objects:
ggpcol = cellMap(standX,
                 mTitle="By column",
                 showrows=showrows,
                 showVals="D",
                 D=tempX,
                 adjustrowlabels=0.5)  
plot(ggpcol)

## ----fig.height=10,fig.width=8------------------------------------------------
ggpDDC = cellMap(DDCtransTG$stdResid, 
                 mTitle="DetectDeviatingCells",
                 showrows=showrows,
                 showVals="D",
                 D=tempX,
                 adjustrowlabels=0.5)
plot(ggpDDC)

# Creating the pdf:
# pdf("cellMap_TopGear.pdf", width = 20, height = 15 )
# gridExtra::grid.arrange(ggpcol,ggpDDC,nrow=1) # combines 2 plots in a figure
# dev.off()


## ----fig.height=8,fig.width=6-------------------------------------------------
# Top Gear dataset: prediction of "new" data
############################################
# For comparison we first remake the cell map of the entire dataset, but now 
# showing the values of the residuals instead of the data values:

dim(remX) # 296 11

ggpDDC = cellMap(DDCtransTG$stdResid, 
                 mTitle="DetectDeviatingCells",
                 showrows=showrows,
                 showVals="R",
                 adjustrowlabels=0.5)
plot(ggpDDC)

## -----------------------------------------------------------------------------
initX = remX[-showrows,]
dim(initX) # 279 11

# Fit initX:
DDCinitX = DDC(initX,DDCpars=DDCpars) 

## -----------------------------------------------------------------------------
newX = remX[showrows,]
dim(newX) # 17  11 

# Make predictions by DDCpredict. 
# Its inputs are:
#   Xnew       : the new data (test data)
#   InitialDDC : Must be provided.
#   DDCpars    : the input options to be used for the prediction.
#                By default the options of InitialDDC are used. 

predictDDC = DDCpredict(newX,DDCinitX)

names(DDCinitX)
# For comparison with:

names(predictDDC) # Fewer, since DDCpredict does not call checkDataSet:

# If you specify the parameters the result is the same:
predictDDC2 = DDCpredict(newX,DDCinitX,DDCpars=DDCpars)
all.equal(predictDDC,predictDDC2) # TRUE

## ----fig.height=10,fig.width=8------------------------------------------------

ggpnew = cellMap(predictDDC$stdResid,
                 showVals="R",
                 mTitle="DDCpredict",
                 adjustrowlabels=0.5) 
plot(ggpnew) # Looks quite similar to the result using the entire dataset:

# Creating the pdf:
# pdf("TopGear_DDCpredict.pdf",width=20,height=15)
# gridExtra::grid.arrange(ggpDDC,ggpnew,nrow=1) 
# dev.off()

## -----------------------------------------------------------------------------
data(data_philips)
dim(data_philips) 
colnames(data_philips) = c("X1","X2","X3","X4","X5","X6","X7","X8","X9")

DDCphilips = DDC(data_philips)

qqnorm(as.vector(DDCphilips$Z)) # rather gaussian, here we only see 2 outliers:

round(DDCphilips$stdResid[1:12,],1) # the standardized residuals:

DDCphilips$indcells # indices of the cells that were flagged:

DDCphilips$indrows # flagged rows:


## ----results='hide',message=FALSE,warning=FALSE-------------------------------
library(robustbase) # for covMcd

## ----fig.height=4,fig.width=8-------------------------------------------------
MCDphilips = robustbase::covMcd(data_philips)
indrowsMCD = which(mahalanobis(data_philips,MCDphilips$center,
                               MCDphilips$cov) > qchisq(0.975,df=9))

plot(sqrt(mahalanobis(data_philips,MCDphilips$center,MCDphilips$cov)),
     main="Philips data",ylab="Robust distances",xlab="",pch=20)
abline(h=sqrt(qchisq(0.975,df=9))) # this horizontal line is the cutoff.
# dev.copy(pdf,"Figure_philips_left.pdf",width=10,height=4)
# dev.off()

## ----fig.height=8,fig.width=6-------------------------------------------------
# cellMaps with rectangular blocks:

ggpMCDphilips = cellMap(data_philips,
                        indrows=indrowsMCD,
                        mTitle="MCD",
                        nrowsinblock=15,
                        ncolumnsinblock=1,
                        drawCircles = TRUE)
plot(ggpMCDphilips)

ggpDDCphilips = cellMap(DDCphilips$stdResid,
                        indrows=DDCphilips$indrows,
                        mTitle="DetectDeviatingCells",
                        nrowsinblock=15,
                        ncolumnsinblock=1,
                        drawCircles = TRUE)
plot(ggpDDCphilips)
# dev.copy(pdf,"Figure_philips_right.pdf",width=6,height=12)
# dev.off()

## -----------------------------------------------------------------------------
data(data_mortality)
dim(data_mortality)
# 198  91
rownames(data_mortality)[1:5] 
colnames(data_mortality)[1:5] 

DDCpars = list(fastDDC = FALSE, silent = TRUE)
DDCmortality = DDC(data_mortality,DDCpars) # 1 second

remX = DDCmortality$remX
dim(remX)

## ----results='hide',message=FALSE,warning=FALSE-------------------------------
library(rrcov) # contains ROBPCA

## ----fig.height=8,fig.width=6-------------------------------------------------
PCAmortality = rrcov::PcaHubert(data_mortality,alpha=0.75,scale=FALSE)

ggpROBPCA = cellMap(remX,
                    indrows=which(PCAmortality@flag==FALSE),
                    mTitle="By row",
                    nrowsinblock=5,
                    ncolumnsinblock=5, 
                    rowtitle = "Years",
                    columntitle = "Age",
                    sizetitles = 2.0,
                    drawCircles = TRUE)

plot(ggpROBPCA)

ggpDDC = cellMap(DDCmortality$stdResid,
                 indrows=DDCmortality$indrows,
                 mTitle="DetectDeviatingCells",
                 nrowsinblock=5,
                 ncolumnsinblock=5,
                 rowtitle = "Years",
                 columntitle = "Age",
                 sizetitles = 2.0,
                 drawCircles = TRUE)
plot(ggpDDC) # Leads to a detailed interpretation:

# pdf("cellmap_mortality.pdf",width=14,height=12)
# gridExtra::grid.arrange(ggpROBPCA,ggpDDC,nrow=1)
# dev.off()

## ----fig.height=8,fig.width=6-------------------------------------------------
rowblocksizes = c(84,14,5,21,6,35,33)
rowlabels = c("19th century","1900-1913","WW1","1919-1939",
"WW2","1946-1980","recent")
colblocksizes = c(5,5,10,10,10,10,10,10,10,10,1)
collabels = c("upto 4","5 to 9","10 to 19","20 to 29","30 to 39","40 to 49","50 to 59","60 to 69","70 to 79","80 to 89","90+")

ggpDDC = cellMap(DDCmortality$stdResid,
                 mTitle = "Cellmap with manual blocks",
                 manualrowblocksizes = rowblocksizes,
                 rowblocklabels = rowlabels,
                 manualcolumnblocksizes = colblocksizes,
                 columnblocklabels = collabels,
                 rowtitle = "Epochs",
                 columntitle = "Age groups",
                 sizetitles = 2.0)
# pdf("cellmap_mortality_manual_blocks.pdf",width=8,height=5)
plot(ggpDDC)
# dev.off()

## -----------------------------------------------------------------------------
data(data_glass)
DDCpars = list(fastDDC = FALSE, silent = TRUE)
DDCglass = DDC(data_glass,DDCpars) # takes 8 seconds
remX = DDCglass$remX
# With DDCpars$silent = FALSE we obtain more information:
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

dim(remX)

## -----------------------------------------------------------------------------

fastDDCpars = list(fastDDC = TRUE, silent = TRUE)
fastDDCglass = DDC(data_glass, fastDDCpars) # takes 2 seconds
remXfast = fastDDCglass$remX
all.equal(remX,remXfast) # The remaining data is the same:

## ----results='hide',message=FALSE,warning=FALSE-------------------------------
library(rrcov) # contains ROBPCA

## ----fig.height=4,fig.width=8-------------------------------------------------
PCAglass = rrcov::PcaHubert(remX,alpha=0.75,scale=FALSE)

n = nrow(remX)
nrowsinblock = 5
rowtitle = "glass samples"
rowlabels = rep("",floor(n/nrowsinblock));
rowlabels[1] = "1"
rowlabels[floor(n/nrowsinblock)] = "n";

d = ncol(remX)
ncolumnsinblock = 5
columntitle = "wavelengths"
columnlabels = rep("",floor(d/ncolumnsinblock));
columnlabels[1] = "1";
columnlabels[floor(d/ncolumnsinblock)] = "d"

ggpROBPCA = cellMap(matrix(0,n,d),
                    indrows=which(PCAglass@flag==FALSE),
                    rowblocklabels=rowlabels,
                    columnblocklabels=columnlabels,
                    mTitle="By row",
                    nrowsinblock=5,
                    ncolumnsinblock=5,
                    columnangle=0,
                    rowtitle="glass samples",
                    columntitle="wavelengths",
                    sizetitles=1.5,
                    drawCircles = TRUE)
plot(ggpROBPCA)

ggpDDC = cellMap(DDCglass$stdResid,
                 indrows=DDCglass$indrows,
                 rowblocklabels=rowlabels,
                 columnblocklabels=columnlabels,
                 mTitle="DDC",
                 nrowsinblock=5,
                 ncolumnsinblock=5,
                 columnangle=0,
                 rowtitle="glass samples",
                 columntitle="wavelengths",
                 sizetitles=1.5,
                 drawCircles = TRUE)
plot(ggpDDC)
# pdf("cellmap_glass_ROBPCA_DDC.pdf",width=16,height=10)
# gridExtra::grid.arrange(ggpROBPCA,ggpDDC,ncol=1)
# dev.off()

ggpfastDDC = cellMap(fastDDCglass$stdResid,
                     indrows=fastDDCglass$indrows,
                     rowblocklabels=rowlabels,
                     columnblocklabels=columnlabels,
                     mTitle="fast DDC",
                     nrowsinblock=5,
                     ncolumnsinblock=5,
                     columnangle=0,
                     rowtitle="glass samples",
                     columntitle="wavelengths",
                     sizetitles=1.5,
                     drawCircles = TRUE)
plot(ggpfastDDC)
# pdf("cellmap_glass_DDC_fastDDC.pdf",width=16,height=10)
# gridExtra::grid.arrange(ggpDDC,ggpfastDDC,ncol=1)
# dev.off()

