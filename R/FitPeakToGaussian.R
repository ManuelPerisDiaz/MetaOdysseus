#' FitpeakToGaussian
#'
#' Using the output from EstimateGaussianParameters, do a nonlinear
#' fit to a Gaussian model
#'
#' @param x vector with the massess
#' @param y vector with the intensity
#' @param centroid mean of the peak
#' @param doPlot - A boolean (default FALSE) to plot the data, estimate, and the fit.
#' @param pTitle - An optional plot title
#' @param xTitle - An optional X-axis label
#' @param yTitle  - An optional Y-axis label
#'
#' @return sum - A summary of the fit
#' @export
#'
#' @examples
#' # Not run
FitPeakToGaussian <- function(x,y, centroid, doPlot=TRUE, pTitle='Peak',xTitle='x', yTitle='y'){

  Mu <- function(x,y){
    l <- length(x)
    sy  <- 0.
    sxy <- 0.
    for (i in 1:l){
      sy <- sy + y[i]
      sxy <- sxy +  x[i]*y[i]
    }
    muEst <- sxy /sy
    muEst
  }
Sigma <- function(x,y, mu){
    l <- length(x)
    sy  <- 0.
    sd2 <- 0.
    for (i in 1:l){
      sy <- sy + y[i]
      sd2 <- sd2 +  y[i] * (x[i]-mu)^2
    }
    varEst <- sd2/ (sy - 1.)
    sdEst = sqrt(varEst)
    sdEst
  }

Mu<-Mu(x,y)
Sigma<-Sigma(x,y,Mu)
# extract the data we need
  l <- length(x)
  i <- which(x==centroid)
  htEst=y[i]


  muEst    <- Mu
  sigmaEst <- Sigma
  htEst    <- htEst




  nls.control <-nls.control(maxiter = 500, tol = 1e-01, minFactor = 1/4.9152e+11,
              printEval = FALSE, warnOnly = FALSE)
 res <- nls(y ~ scale*exp(-0.5*(x-mu)^2/sigma^2),control=nls.control,
             start=c(mu=muEst, sigma=mean(sigmaEst), scale=mean(htEst)))

 # gstart <- data.frame(mu=muEst,scale=mean(htEst),sigma= seq(from=0.1, to=max(lData$sigmaEst),by=0.2))


  # res<- nls2(y ~ scale*exp(-0.5*(x-mu)^2/sigma^2),
  # start=gstart)



  sum <- summary(res)
  ht <- sum$coefficients[3]
  mu <- sum$coefficients[1]
  sigma <- sum$coefficients[2]

  if(doPlot==TRUE){
    yc <- mean(htEst)*exp(-0.5*(x-muEst)^2/sigmaEst^2)
    yc2 <- ht*exp(-0.5*(x-mu)^2/sigma^2)
    xMin <- min(x)
    xMax <- max(x)
    yMax <- max(y)

    tiff("fitteed.tiff",width=4000, height=2500,res=500)
      par(mar=c(5,6,4,1)+.1)


    plot(c(xMin, xMax), c(0, yMax), type='n',
         xlab=xTitle, ylab=yTitle,  cex.lab=2,cex.axis=1.5)
         points(x,y,pch=19)
         lines(x,yc, col='red', lw=2)
         lines(x,yc2, col='blue', lw=2)
         legend("topright", bty="n",
         legend=c("points", "estimate", "fit"),
         col=c('black', 'red', 'blue'),
         lty=c(0, 1, 1),
         lwd=c(0, 1,1),
         pch=c(19, NA, NA))

         dev.off()

  }
  return(sum)
  }
