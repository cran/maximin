\name{maximin}
\alias{maximin}
\alias{maximin.cand}

\title{
  Sequential space-filling design under the criterion of maximin distance
}

\description{
  Generates sequential space-filling design under the criterion of maximin distance;
  both discrete and continuous searches are provided.
}

\usage{
  maximin.cand(n, Xcand, Tmax, Xorig=NULL, init=NULL, verb=FALSE) 
  maximin(n, p, T, Xorig=NULL, Xinit=NULL, verb=FALSE, plot=FALSE, boundary=FALSE)
}

\arguments{
  \item{n}{ the number of space-filling locations }
  \item{Xcand}{ the candidate set, from which each space-filling location is selected }
  \item{Tmax}{ the number of iterations; \code{Tmax <= nrow(Xcand)}; to be safe, set \code{Tmax = nrow(Xcand)}. }
  \item{Xorig}{ the existing design; \code{ncol(Xorig) = ncol(Xcand)} }
  \item{init}{ the initial indices of \code{X}; it can be randomly selected from \code{Xcand} or introduced from a previous experiment. }
  \item{verb}{ progress indicator --- every tenth iteration is printed out; by default \code{verb = FALSE}. }
  \item{p}{ the dimensionality of input space }
  \item{T}{ the number of iterations; \eqn{T > n}{T > n}; setting \eqn{T = 10 * n}{T = 10 * n} is a good starting point. }
  \item{Xinit}{ the (initial) design introduced from a previous experiment }
  \item{plot}{ if \code{plot = TRUE}, then the search space and the "start location --> new location" with each iteration is plotted; 
               if \eqn{p > 2}{p > 2}, then TWO input coordinates are RANDOMLY chosen for plotting; it is worth noticing that the search space 
               only VISUALLY makes sense when \eqn{p = 2}{p = 2}. }
  \item{boundary}{ if \code{boundary = TRUE}, then for each iteration, the "to-be-swapped-in" location will be away from the 
                   boundary in addition to being away from other \code{X} locations and \code{Xorig}; how far is it? \eqn{min(d, 4*d.bound)}{min(d, 4*d.bound)}, 
                   where \code{d} is the Euclidean distance between the "to-be-swapped-in" location and other \code{X} locations as well as \code{Xorig}, while 
                   \code{d.bound} is the minimum Euclidean distance between the "to-be-swapped-in" location and the boundaries.}
}

\details{

  For a visualization of applying \code{maximin.cand} in a real-life problem, see Sun et al. (2018).
  
  \code{maximin.cand} returns the space-filling design, the indices of \code{Xcand}, the minimum pairwise Euclidean distance as well as the corresponding number 
  
  \code{maximin} returns the combined existing design and the space-filling design as well as the minimum pairwise Euclidean distance through all the iterations
}

\value{
\code{maximin.cand} returns
  \item{\code{X}}{ the space-filling design; \code{dim(X) = n * ncol(Xcand)} }
  \item{\code{inds}}{ the indices of \code{Xcand}, which makes \code{X} }
  \item{\code{mis}}{ the minimum distance with each iteration; \code{length(mis) = T + 1} }
  \item{\code{mislen}}{ the number of minimum distance with each iteration; \code{length(mislen) = T + 1} }
  
\code{maximin} returns
  \item{\code{Xf}}{ \code{dim(Xf) = (nrow(Xorig) + n) * p} }
  \item{\code{mi}}{ the minimum distance with each iteration; \code{length(mi) = T + 1} }
}

\references{
  
  F. Sun, R.B. Gramacy, B. Haaland, S.Y. Lu, and Y. Hwang (2018)
  \emph{Synthesizing simulation and field data of solar irradiance};
  preprint on arXiv:1806.05131. \url{http://arxiv.org/abs/1806.05131}
  
  M.H.Y. Tan (2013)
  \emph{Minimax Designs for Finite Design Regions}, Technometrics, 55(3), 346-358. 
  
  M.E. Johnson, L.M. Moore, and D. Yivisaker (1990)
  \emph{Minimax and maximin distance designs}, Journal of Statistical Planning and Inference, 26(2), 131-148. 
}

\author{
  Furong Sun \email{furongs@vt.edu} and Robert B. Gramacy \email{rbg@vt.edu} 
}

\note{
  Regarding \code{maximin.cand}, it is unlikely that \code{Tmax} can be reached, and therefore, 
  loading the \code{temp.RData} from the working directory is usually necessary. 
}

\examples{
   
 \dontrun{
   ## maximin.cand
   # a generic function to expand grids at higher dimension
   expand.grids <- function(x,d) expand.grid(replicate(d, x, simplify=FALSE))
   
   # generate the design
   library(lhs)
   n <- 10
   p <- 2
   Xorig <- as.data.frame(randomLHS(10, p))
   x1 <- seq(0, 1, length.out=100)
   Xcand <- expand.grids(x1, p)
   T <- nrow(Xcand)
   Xsparse <- maximin.cand(n=n, Xcand=Xcand, Tmax=T, Xorig=Xorig, init=NULL, verb=FALSE)
   
   # load the saved data from the working directory
   load("temp.RData")
   maxmd <- as.numeric(format(round(max(na.omit(mind)), 5), nsmall=5))
   
   # visualization
   par(mfrow=c(1,2))
   plot(X[,1], X[,2], xlab=expression(x[1]), ylab=expression(x[2]), xlim=c(0,1), ylim=c(0,1), 
        main=paste0("n=", nrow(X), "_p=", ncol(X), "_T=", t, "_maximin=", maxmd))
   points(Xorig, col=2, pch=20)
   abline(h=c(0,1), v=c(0,1), lty=2, col=2)
   if(!is.null(Xorig)) legend("topright", "Xorig", xpd=TRUE, horiz=TRUE, 
      inset=c(-0.03, -0.05), pch=20, col=2, bty="n")
   plot(log(na.omit(mind)), type="b", xlab="iteration", ylab="log(minimum distance)", 
        main="progress on minimum distance")
   abline(v=nrow(X), lty=2)
   mtext(paste0("design size=", nrow(X)), at=nrow(X), cex=0.6)
  }
  
  ## maximin
  # generate the design
  library(lhs)
  n <- 10
  p <- 2
  T <- 10*n
  Xorig <- randomLHS(10,p)
  Xsparse <- maximin(n=n, p=p, T=T, Xorig=Xorig, Xinit=NULL, 
                     verb=FALSE, plot=FALSE, boundary=FALSE)
  maxmd <- as.numeric(format(round(Xsparse$mi[T+1], 5), nsmall=5))
  
  # visualization
  par(mfrow=c(1,2))
  plot(Xsparse$Xf[,1], Xsparse$Xf[,2], xlab=expression(x[1]), ylab=expression(x[2]), 
       xlim=c(0,1), ylim=c(0,1), 
       main=paste0("n=", n, "_p=", p, "_T=", T, "_maximin=", maxmd))
  points(Xorig, col=2, pch=20)
  abline(h=c(0,1), v=c(0,1), lty=2, col=2)
  if(!is.null(Xorig)) legend("topright", "Xorig", xpd=TRUE, horiz=TRUE, 
     inset=c(-0.03, -0.05), pch=20, col=2, bty="n")
  plot(log(Xsparse$mi), type="b", xlab="iteration", ylab="log(minimum distance)", 
       main="progress on minimum distance")
  abline(v=n, lty=2)
  mtext(paste0("design size=", n), at=n, cex=0.6)
  abline(v=T, lty=2)
  mtext(paste0("max.md=", maxmd), at=T, cex=0.6)
}

\keyword{space-filling design}
\keyword{sequential design}
\keyword{maximin}
