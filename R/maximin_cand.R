#*******************************************************************************
#
# Space-filling Design under Maximin Distance
# Copyright (C) 2018, Virginia Tech
#
# This library is a free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#
# Questions? Contact Furong Sun (furongs@vt.edu) and Robert B. Gramacy (rbg@vt.edu)
#
#*******************************************************************************

## maximin.cand:
##
## generates a space-filling design in a FINITE DESIGN SPACE 
## under the criterion of maximin distance; the candidate pool is pre-specified.

maximin.cand <- function(n, Xcand, Tmax=nrow(Xcand), Xorig=NULL, init=NULL, verb=FALSE, tempfile=NULL) 
{
  ## sanity checks
  if(class(Xcand) != "matrix"){
     Xcand <- as.matrix(Xcand)
  }
  
  if(!is.null(Xorig) && class(Xorig) != "matrix"){
     Xorig <- as.matrix(Xorig)
  }
  
  if(n==1){
     stop("n must be bigger than 1.")
  }
  
  if(!is.null(Xorig))
     if(ncol(Xcand) != ncol(Xorig)) 
        stop("column dimension mismatch between Xcand and Xorig :-(")
  
  if(Tmax <= n){
    warning("Tmax had better be bigger than n.")
  }
  
  ## the number of runs in the candidate set
  ncand <- nrow(Xcand)
  
  ## the indices (Xcand) of the initial design
  if(!is.null(init)){   
     xi <- init  
  }else{
     xi <- sample(1:ncand, n)
  }
  
  X <- Xcand[xi,]            ## the initial design, which is a subset of the candidate set
  
  xo <- setdiff(1:ncand, xi) ## the remaining indices (Xcand) in the candidate set, after the initial design is extracted.
  
  ## (Euclidean) distances among runs in X
  D <- distance(X)                         ## an n * n matrix
  md <- min(as.numeric(D[upper.tri(D)]))   ## the 'md' among X: no diagonal elements (since the diagonal elements of D are ZERO)
  md.ind <- which(D==md, arr.ind=TRUE)[,1] ## the indices in X with the 'md' among runs in X
  #mdlen <- length(md.ind)/2               ## the number of 'md'
  
  ## the "true" candidate set, after removing the initial design
  Xun <- Xcand[-xi,]
  ## obtain the distance matrix between Xun and X
  Du <- distance(Xun, X)
  ## obtain the indices of the runs in Xun with distances to X bigger than "md": "to-be-swapped-in" runs
  uw.ind <- which(rowSums(Du>md) == ncol(Du))
  
  ## distances to the original design given it exists
  if(!is.null(Xorig)){
     D2 <- distance(X, Xorig)  ## an nrow(X) * nrow(Xorig) matrix
     md2 <- min(D2)            ## the `md' between X and Xorig
     
     ## make a decision: use the smaller `md': from local minimum to global minimum
     if(md2 < md){
        md <- md2
        md.ind <- which(D2==md2, arr.ind=TRUE)[,1] ## the indices in X with the `md' between X and Xorig
        #mdlen <- length(md.ind)
     }else if(md2 == md){
        md.ind <- unique(c(md.ind, which(D2==md2, arr.ind=TRUE)[,1]))
        #mdlen <- mdlen + sum(D2==md2)
     }
     
     ## the distance matrix between Xun and Xorig: dimensionality of nrow(Xun) * nrow(Xorig)
     # Du2 <- cbind(Du, distance(Xun, Xorig))
     # uw.ind <- which(rowSums(Du2 > md)==ncol(Du2))
     Du2.newcols <- distance(Xun, Xorig) ## to avoid a duplicate of Du---saving memory!!!
     uw.ind <- which(rowSums(Du > md)==ncol(Du) & rowSums(Du2.newcols > md)==ncol(Du2.newcols))
  }
  
  ## allocate space for 'md' and 'mdlen'
  mind <- rep(NA, Tmax+1)
  #mindlen <- rep(NA, Tmax+1)
  mind[1] <- md 
  #mindlen[1] <- mdlen
  ## there should be improvement with each iteration regarding design criterion: mdprime > md
  for(t in 1:Tmax){
      
      if(length(uw.ind) == 0){
         warning("terminated early since the maximum progress has been achieved :-)")
         return(list(inds=xi, mis=mind))#, mislen=mindlen))
      }
      
      ## randomly select a location (index) with md: to-be-swapped-out
      row.in.ind <- ceiling(runif(1)*length(md.ind)) 
      row.in <- md.ind[row.in.ind] 
      
      ## the "to-be-swapped-out" location: corresponding to row.in
      xold <- matrix(X[row.in,], nrow=1) 
      
      ## randomly select a location (index) from the subsetted candidate set, uw.ind, which guarantees progress
      row.out.ind <- ceiling(runif(1)*length(uw.ind)) 
      row.out <- uw.ind[row.out.ind]
      X[row.in,] <- Xcand[xo[row.out],] ## == Xun[row.out,], the "to-be-swapped-in" location: corresponding to row.out
      
      if(class(X[-row.in,]) != "matrix"){
          X.rrow <- matrix(X[-row.in,], nrow=1)
      }
      
      Xr <- matrix(X[row.in,], nrow=1)
      
      if(class(X[-row.in,]) == "matrix"){
         dr <- distance(Xr, X[-row.in,])
      }else{
         dr <- distance(Xr, X.rrow)
      }
        
      
      ## update D
      D[row.in, -row.in] <- D[-row.in, row.in] <- as.numeric(dr)
      
      dprime <- as.numeric(D[upper.tri(D)])
      mdprime <- min(dprime)
      mdprime.ind <- which(D==mdprime, arr.ind=TRUE)[,1]
      #mdprimelen <- length(mdprime.ind)/2
      
      ## update Xun
      Xun[row.out,] <- Xcand[xi[row.in],] ## == as.numeric(xold)
      ## update Du
      Du[,row.in] <- as.numeric(distance(Xr, Xun))
      Du[row.out,] <- as.numeric(distance(xold, X))
    
    ## distances to the original design given it exists
    if(!is.null(Xorig)){
      
       ## update D2
       D2[row.in,] <- as.numeric(distance(Xr, Xorig))
       md2prime <- min(D2)
      
       ## make a decision: use the smaller "new" md: from local "new" minimum to global "new" minimum
       if(md2prime < mdprime){
          mdprime <- md2prime
          mdprime.ind <- which(D2 == md2prime, arr.ind=TRUE)[,1]
          #mdprimelen <- length(mdprime.ind)
       }else if(md2prime == mdprime){
          mdprime.ind <- unique(c(mdprime.ind, which(D2==md2prime, arr.ind=TRUE)[,1]))
          #mdprimelen <- mdprimelen + sum(D2==md2prime)
       }
      
      ## update Du2
      Du2.newcols[row.out,] <- as.numeric(distance(xold, Xorig))
      uwprime.ind <- which(rowSums(Du > mdprime)==ncol(Du) & rowSums(Du2.newcols > mdprime)==ncol(Du2.newcols))
      
    }else{
      uwprime.ind <- which(rowSums(Du > mdprime)==ncol(Du))
    }
    
    ## sanity check
    if((mdprime < md))# || ((mdprime == md) && (mdprimelen >= mdlen))) 
        stop("md should increase with each iteration :-|")
    
    ## update indices
    xiold <- xi[row.in]
    xi[row.in] <- xo[row.out]
    xo[row.out] <- xiold
    
    md <- mdprime             
    #mdlen <- mdprimelen       
    md.ind <- mdprime.ind     
    uw.ind <- uwprime.ind
    
    mind[t+1] <- md
    #mindlen[t+1] <- mdlen
    
    if(!is.null(tempfile)) 
       save(xi, mind, t, file=tempfile)#mindlen, t, file=tempfile) ## no need to save "X" since "xi" is enough to construct the final design
    
    # progress indicator
    if(verb == TRUE)
       if(t%%10 == 0)
          cat("t=", t, "/Tmax=", Tmax, " is done.\n", sep="")
    
  }
  
  # if all the iterations are done ...
  return(list(inds=xi, 
               mis=mind))#, 
              #mislen=mindlen))
}
