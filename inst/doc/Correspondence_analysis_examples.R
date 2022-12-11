## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 8 ,
  fig.height = 12,
  fig.align ='center'
)

## -----------------------------------------------------------------------------
library(cellWise)

## -----------------------------------------------------------------------------
data("data_clothes")
X <- data_clothes
dim(X)
# create matrix S as in paper:
sum(X) # total count is 4373
P    <- X / sum(X) # relative frequencies that add up to 1
rvec <- rowSums(P) # vector of row totals
cvec <- colSums(P) # vector of column totals
R    <- X / rowSums(X) # row profiles
rowSums(R) # all 1,  OK
S <- diag(sqrt(rvec)) %*% (R - rep(1, nrow(X)) %*% 
                             t(cvec)) %*% diag(1/sqrt(cvec))
dimnames(S) <- dimnames(X)
round(S, 3)
d <- ncol(S)

# We verify that the points of S are on a hyperplane by construction:
(eigvals <- eigen(t(S) %*% S)$values)


## -----------------------------------------------------------------------------
DDC.out <- DDC(S)
ggp <- cellMap(R = DDC.out$stdResid,               
               indcells = DDC.out$indcells,
               indrows = DDC.out$indrows,
               mTitle = "",
               rowtitle = "countries",
               columntitle = "brackets    ",
               sizetitles = 2,
               drawCircles = F)

# pdf("clothes_cellmap.pdf", width = 4, height = 8)
plot(ggp)
# dev.off()

## -----------------------------------------------------------------------------
countryInds <- which(rownames(S) %in% 
                       c("GB", "SK", "MT", "GR", "BG", "LV", "RO"))
matplot(t(S[countryInds, ]), pch = 16, col = 1:7)
lines(apply(S, 2, median), lwd = 3) # median profile
# The profile of these countries deviate in a similar way: 
# they trade a lot of cheap clothes, and fewer expensive ones. 

## -----------------------------------------------------------------------------
svd.out <- svd(S)
svd.out$d# S is indeed singular:
(svd.out$d)^2 - eigen(t(S) %*% S)$values # ~0
diff(c(0,cumsum(svd.out$d^2)/sum(svd.out$d^2)))
# This can be shown in a scree plot.

# For plotting the rows:
rowproj <- diag(1/sqrt(rvec)) %*% svd.out$u %*% diag(svd.out$d)
# For plotting the columns = variables:
colproj <- diag(1/sqrt(cvec)) %*% svd.out$v %*% diag(svd.out$d)

# pdf(file="clothes_ClassCorrespA.pdf", width=7, height=6)
plot(rowproj[, 1:2], col = "white", xlab="", ylab="",
     xlim = c(-1, 0.6), ylim = c(-0.6, 0.6))
title(main="Classical correspondence analysis of clothes data", 
      line=1)
text(rowproj[, 1:2], labels = rownames(S))
abline(v=0); abline(h=0)
text(colproj[, 1:2], labels = colnames(S), col="red")
facs = c(0.93,0.85,0.78,0.84,0.87)
shape::Arrows(0, 0, facs*colproj[, 1], facs*colproj[, 2], 
              col = "red", arr.type="simple", arr.width=0.1,
              arr.length = 0.1) # , arr.adj = 0)
title(xlab="Dimension 1", line=2.3)
title(ylab="Dimension 2", line=2.3)
# dev.off()

## -----------------------------------------------------------------------------
# We apply MacroPCA with center fixed at zero.
# As in classical CA, we do not prescale S.
MacroPCApar0 <- list(scale = FALSE, center = rep(0,d))
MacroPCA.out <- MacroPCA(S, k=0, MacroPCApars = MacroPCApar0)

MacroPCA.out <- MacroPCA(S, k=2, MacroPCApars = MacroPCApar0)

(eigvals <- MacroPCA.out$eigenvalues)
diff(c(0,cumsum(eigvals)/sum(eigvals)))

# Compute the principal coordinates for the biplot.
# To make the plot match the orientation in Riani et al:
V <- -MacroPCA.out$loadings
Tscores <- -MacroPCA.out$scores
sVals <- sqrt(nrow(S)*MacroPCA.out$eigenvalues)
U <- Tscores %*% diag(1/sVals)
rowproj <- diag(1/sqrt(rvec)) %*% U %*% diag(sVals)
colproj <- diag(1/sqrt(cvec)) %*% V %*% diag(sVals)

# pdf(file="clothes_MacroCA_biplot.pdf", width=7, height=6)
plot(rowproj[, 1:2], col = "white", xlim=c(-0.95,0.65),
     ylim=c(-0.6,0.6), xlab="", ylab="")
title(main="Cellwise robust correspondence analysis of clothes data", 
      line=1)
text(rowproj[,1:2], labels = rownames(S))
abline(h=0); abline(v=0)
text(colproj[, 1:2], labels = colnames(S), col="red")
facs = c(0.9,0.8,0.5,0.75,0.85)
shape::Arrows(0, 0, facs*colproj[, 1], facs*colproj[, 2], 
              col = "red", arr.type="simple", arr.width=0.1,
              arr.length = 0.1)
title(xlab="Dimension 1", line=2.3)
title(ylab="Dimension 2", line=2.3)
# dev.off()
# Matches Fig 4 of Riani quite well.

## -----------------------------------------------------------------------------
data("data_brands")
X <- data_brands
dim(X) 
sum(X) # total count is 11713
P    <- X/sum(X) # relative frequencies that add up to 1
rvec <- rowSums(P) # vector of row totals
hist(rvec) # Right tail: Chevrolet, Ford, Honda, Toyota
# These brands are well known and sold a lot in the US.
cvec <- colSums(P) # vector of column totals
R    <- X / rowSums(X) # row profiles
S    <- diag(sqrt(rvec)) %*% (R - rep(1, nrow(X)) %*% 
                                t(cvec)) %*% diag(1/sqrt(cvec))
dimnames(S) <- dimnames(X)
d <- ncol(S)

## -----------------------------------------------------------------------------
DDC.out <- DDC(S)
ggp     <- cellMap(R = DDC.out$stdResid,                        
                    indcells = DDC.out$indcells,
                    indrows = DDC.out$indrows,
                    mTitle = "",
                    rowtitle = "brands",
                    columntitle = "perceptions   ",
                    sizetitles = 2.3,
                    drawCircles = F)
#pdf("brands_cellmap.pdf", width = 6, height = 13)
plot(ggp)
#dev.off()

# Volvo is most deviating (3 cells), followed by Hyundai
# (2 cells) and Maserati (2 cells).

## -----------------------------------------------------------------------------
svd.out <- svd(S)
svd.out$d # S is singular:
diff(c(0,cumsum(svd.out$d^2)/sum(svd.out$d^2)))
# Can be plotted in a scree plot.

# To match the plot in Riani at al:
svd.out$v[, 2] = -svd.out$v[, 2]
svd.out$u[, 2] = -svd.out$u[, 2]
rowproj <- diag(1/sqrt(rvec)) %*% svd.out$u %*% diag(svd.out$d)
colproj <- diag(1/sqrt(cvec)) %*% svd.out$v %*% diag(svd.out$d)


# pdf(file="brands_ClassCorrespA.pdf", width=7, height=6)
plot(rowproj[, 1:2], col = "white", xlim=c(-1.1,1.1),
     ylim=c(-0.8,1.5), xlab="", ylab="")
title(main="Classical correspondence analysis of brands data", 
      line=1)
text(rowproj[,1:2], labels = rownames(S))
abline(v=0); abline(h=0)
text(colproj[, 1:2], labels = colnames(S), col="red")
facs = c(0.8,0.9,0.8,0.65,0.92,0.8,0.65)
shape::Arrows(0, 0, facs*colproj[, 1], facs*colproj[, 2], 
              col = "red", arr.type="simple", arr.width=0.1,
              arr.length = 0.1)
title(xlab="Dimension 1", line=2.3)
title(ylab="Dimension 2", line=2.3)
# dev.off()

## -----------------------------------------------------------------------------
MacroPCApar0 <- list(scale = FALSE, center = rep(0, d))
MacroPCA.out <- MacroPCA(S, k = 0, MacroPCApars = MacroPCApar0)

MacroPCA.out <- MacroPCA(S, k = 3, MacroPCApars = MacroPCApar0)


# Computations for the biplot.
V <- MacroPCA.out$loadings
scores <- MacroPCA.out$scores
# To make the plot match the orientation in Riani et al:
V[,2]      <- -V[,2]
scores[,2] <- -scores[,2]
sVals      <- sqrt(nrow(S)*MacroPCA.out$eigenvalues)
U          <- scores %*% diag(1/sVals)
rowproj    <- diag(1/sqrt(rvec)) %*% U %*% diag(sVals)
colproj    <- diag(1/sqrt(cvec)) %*% V %*% diag(sVals)

# pdf(file="brands_MacroCA_biplot.pdf", width=7, height=6)
plot(rowproj[, 1:2], col = "white", xlim=c(-1.1,1.1),
     ylim=c(-0.6,0.6), xlab="", ylab="")
title(main="Cellwise robust correspondence analysis of brands data",
      line=1)
text(rowproj[,1:2], labels = rownames(S))
abline(h=0); abline(v=0)
text(colproj[, 1:2], labels = colnames(S), col="red")
facs = c(0.75,0.76,0.9,0.88,0.65,0.74,0.52)
shape::Arrows(0, 0, facs*colproj[, 1], facs*colproj[, 2], 
              col = "red", arr.type="simple", arr.width=0.1,
              arr.length = 0.1)
title(xlab="Dimension 1", line=2.3)
title(ylab="Dimension 2", line=2.3)
# dev.off()
# Roughly matches Figure 7 of Riani et al (2022).

