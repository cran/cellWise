
outlierMap <- function(res,title="Robust PCA",col="black",
                             pch=16,labelOut=TRUE,id=3){
  
  
  SD <- res$SD
  OD <- res$OD
  cutoffSD <- res$cutoffSD
  cutoffOD <- res$cutoffOD
  if (is.null(pch)) {
    pch <- ifelse(is.null(col), 1, 16)
  }
  if (is.null(col)) {
    col <- "black"
  }
  if (is.list(col)) {
    for (i in 1:length(col)) {
      if (i == 1) {
        plot(SD[col[[i]]$index],OD[col[[i]]$index],
             xlab="Score distance",ylab="Orthogonal distance", 
             main=title,pch=pch,col=col[[i]]$col,
             xlim = c(0, max(SD)*1.1),ylim=c(0,max(OD)*1.1))
      }
      points(SD[col[[i]]$index], OD[col[[i]]$index], pch = pch, 
             col = col[[i]]$col)  
    }
  } else {
    plot(SD,OD,xlab="Score distance",ylab="Orthogonal distance", 
         main=title,pch=pch,col=col,xlim=c(0, max(SD)*1.1), 
         ylim=c(0,max(OD)*1.1))
  }
  abline(v = cutoffSD)
  abline(h = cutoffOD)
  if (labelOut) {
    labelDD(SD,OD,id.n.SD=id,id.n.OD=id)
  }
}


labelDD <- function(x,y,id.n.SD=3,id.n.OD=3,off=0.02) 
{ # used in OutlierMap
  xrange <- graphics::par("usr")
  xrange <- xrange[2] - xrange[1]
  if (id.n.SD > 0 && id.n.OD > 0) {
    n <- length(x)
    indexSD <- sort(x, index.return = TRUE)$ix
    indexSD <- indexSD[(n - id.n.SD + 1):n]
    indexOD <- sort(y, index.return = TRUE)$ix
    indexOD <- indexOD[(n - id.n.OD + 1):n]
    lab <- indexOD
    if (is.character(names(y))) {
      lab <- names(y[indexOD])
    }
    graphics::text(x[indexOD] - off * xrange, y[indexOD], lab)
    lab <- indexSD
    if (is.character(names(x))) {
      lab <- names(x[indexSD])
    }
    graphics::text(x[indexSD] - off * xrange, y[indexSD], lab)
  }
}