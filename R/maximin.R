## Questions? Contact Furong Sun (furongs@vt.edu) and Robert B. Gramacy (rbg@vt.edu)

## maximin:
##
## generates a space-filling design sequentially in a CONTINUOUS way under the criterion of maximin distance

maximin <- function(n, p, T, Xorig=NULL, Xinit=NULL, verb=FALSE, plot=FALSE, boundary=FALSE)
{ 
  ## sanity check
  if(!is.null(Xorig) && ncol(Xorig) != p) stop("column mismatch between X and Xorig :-(")
  
  if(T <= n) warning("T had better be bigger than n.")
  
  ## a generic function to calculate the distance between the "to-be-swapped-in" location and the boundary
  d.bound <- function(xnew){
    xnew <- as.numeric(xnew)
    ## the distance to one boundary: should be multiplied by 4, but I want to do it later
    d.k1 <- xnew^2
    ## the distance to the other boundary: ...
    d.k2 <- (1-xnew)^2
    return(4*c(d.k1, d.k2)) ## return all of them: we will calculate the 'md' inside the maximin.obj
  }
  
  ## the objective function to be maximized
  maximin.obj <- function(X, row.out, xrow.in, Xorig=NULL){
    ## the distance between the "to-be-swapped-in" location and other X locations & the existing design
    d <- as.numeric(distance(matrix(xrow.in, nrow=1), rbind(X[-row.out,], Xorig)))
    if(boundary == TRUE) d <- c(d, d.bound(xrow.in))
    return(sqrt(min(d))) ## use the square root of the minimum distance as the objective since such a way provides much less extreme values
  }
  
  ## the search space when optimizing
  sp <- function(xs, xd=NULL, md){
    if(!is.null(xd)){
      ds <- as.numeric(distance(matrix(xs, nrow=1), matrix(xd, nrow=1)))
    }else ds <- md
    ## treat 'xs' as the center and 'sqrt(ds)' as the "radius" to justify a hypercube
    lo <- xs - sqrt(ds) ## lower bound
    up <- xs + sqrt(ds) ## upper bound
    ## guarantee the search space to be within the global input space
    if(length(lo[lo<0]) != 0) lo[lo<0] <- 0
    if(length(up[up>1]) != 0) up[up>1] <- 1
    return(list(lo=lo, up=up, ds=ds)) ## track the distance between 'xs' and 'xd' for future use
  }
  
  ## The following function is to be used when "jumping out of the well":
  ##    to justify the start location, which can not be the same as 'xs', 
  ##                  however, it should be as close to 'xs' as possible;
  ##    'xd' is the closest X or Xorig location to 'xs';
  ##    The start location will be in-between 'xs' and 'xd' and its distance 
  ##    to 'xs' is 5% of 'sqrt(ds)', the Euclidean distance between 'xs' and 'xd':
  st <- function(xs, xd, ds){
    st.prime <- rep(NA, length(xs))
    for(i in 1:length(xs)){
      st.prime[i] <- ifelse(xs[i] > xd[i], xs[i] - 0.05*sqrt(ds), xs[i] + 0.05*sqrt(ds))
      if(xs[i] == xd[i]) st.prime[i] <- xs[i]
    }
    return(st.prime)
  }
  
  ## to check whether each element of a vector is the same
  cvs <- function(x) diff(range(x)) < sqrt(.Machine$double.eps)
  
  ## the initial design: an n*p matrix
  if(!is.null(Xinit)){
    X <- Xinit ## from a previous experiment
  }else X <- matrix(runif(n*p), ncol=p) ## or randomly generated within the global input space
  
  D <- distance(X)
  md <- min(as.numeric(D[upper.tri(D)]))
  md.ind <- which(D==md, arr.ind=TRUE)[,1] ## the list of locations among X with md
  
  if(!is.null(Xorig)){
    D2 <- distance(X, Xorig)  ## an nrow(X) * nrow(Xorig) matrix
    md2 <- min(D2)            ## the 'md' between X and Xorig
    if(md2 < md){
      md <- md2
      md.ind <- which(D2==md2, arr.ind=TRUE)[,1]
    }else if(md2 == md) md.ind <- unique(c(md.ind, which(D2==md2, arr.ind=TRUE)[,1]))
  }
  
  ## allocate the 'md' and 'design' for each iteration
  mind <- rep(NA, T+1); mind[1] <- md
  Xind <- vector(mode="list", T+1); Xind[[1]] <- X
  ## build up the design iteratively
  for(t in 1:T){
    
    ## track the locations with 'non-md'
    cand.ind <- setdiff(1:n, md.ind) 
    
    ## only swap out the location with 'md'
    row.out.ind <- ceiling(runif(1)*length(md.ind))
    row.out <- md.ind[row.out.ind]
    xold <- X[row.out,]
    
    ## justify the search space
    sw.in <- sp(xs=xold, md=md)
    low <- sw.in$lo  ## lower bound
    upp <- sw.in$up  ## upper bound
    
    ## the start location for 'optim' must be in the matrix form
    start <- matrix(xold, nrow=1)
    ## try to maximize the 'md' between the "to-be-swapped-in" location and other X locations and Xorig, 
    ## constrained by the search space
    out <- optim(par=start, fn=maximin.obj, method="L-BFGS-B", 
                 lower=low, upper=upp, control=list(fnscale=-1), 
                 X=X, row.out=row.out, Xorig=Xorig)
    md.ind <- md.ind[!md.ind==row.out] ## update the 'md' list
    
    if(plot==TRUE){
      d <- sample(1:p, 2, replace=FALSE) ## randomly choose 2 coordinates to plot
      Xv <- X[,d]
      plot(Xv, xlab=paste0("x", d[1]), ylab=paste0("x", d[2]), xlim=c(0,1), ylim=c(0,1), main=paste0("md.loc at it=", t))
      points(Xorig[,d], pch=20, col="forestgreen")
      points(xold[d[1]], xold[d[2]], pch=20, col="red")
      rect(low[d[1]], low[d[2]], upp[d[1]], upp[d[2]], col=rgb(0,1,0,alpha=0.2)) ## the search space
      if(!is.null(Xorig)){
        legend("topright", c("xold", "Xorig"), xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col=c("red","forestgreen"), bty="n")
      }else legend("topright", "xold", xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col="red", bty="n")
    }
    
    while(all.equal(as.numeric(start), as.numeric(out$par))==TRUE){
      if(length(md.ind) >= 1){ ## if there are more than one location with 'md' ...    ## 1
        row.out.ind <- ceiling(runif(1)*length(md.ind))
        row.out <- md.ind[row.out.ind]
        xold <- X[row.out,]
        
        sw.in <- sp(xs=xold, md=md)
        low <- sw.in$lo
        upp <- sw.in$up
        
        start <- matrix(xold, nrow=1)
        out <- optim(par=start, fn=maximin.obj, method="L-BFGS-B",
                     lower=low, upper=upp, control=list(fnscale=-1),
                     X=X, row.out=row.out, Xorig=Xorig)
        md.ind <- md.ind[!md.ind==row.out] ## update the 'md' list
        
        if(plot==TRUE){
          Xv <- X[,d]
          plot(Xv, xlab=paste0("x", d[1]), ylab=paste0("x", d[2]), xlim=c(0,1), ylim=c(0,1), main=paste0("md.loc at it=", t))
          points(Xorig[,d], pch=20, col="forestgreen")
          points(xold[d[1]], xold[d[2]], pch=20, col="red")
          rect(low[d[1]], low[d[2]], upp[d[1]], upp[d[2]], col=rgb(0,1,0,alpha=0.2)) # the search space
          if(!is.null(Xorig)){
            legend("topright", c("xold", "Xorig"), xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col=c("red","forestgreen"), bty="n")
          }else legend("topright", "xold", xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col="red", bty="n")
        }
      }else if(length(cand.ind) >= 1){ ## if there is no location with 'md', then randomly go to a location with 'non-md' ...     ## 2
        row.out.ind <- ceiling(runif(1)*length(cand.ind))
        row.out <- cand.ind[row.out.ind]
        xold <- X[row.out,]
        
        ## combine the distance matrices
        if(is.null(Xorig)==TRUE) D2 <- NULL
        DD2 <- cbind(D, D2)
        
        ## justify 'xd'
        Dor <- order(DD2[row.out,], decreasing=FALSE)
        xd <- ifelse(rep(Dor[2] <= n, p), X[Dor[2],], Xorig[Dor[2]-n,])
        
        ## justify the search space
        sw.in <- sp(xs=xold, xd=xd, md=NULL)
        low <- sw.in$lo
        upp <- sw.in$up
        
        start <- matrix(xold, nrow=1)
        out <- optim(par=start, fn=maximin.obj, method="L-BFGS-B",
                     lower=low, upper=upp, control=list(fnscale=-1),
                     X=X, row.out=row.out, Xorig=Xorig)
        cand.ind <- cand.ind[!cand.ind==row.out] ## update the 'non-md' list
        
        if(plot==TRUE){
          Xv <- X[,d]
          plot(Xv, xlab=paste0("x", d[1]), ylab=paste0("x", d[2]), xlim=c(0,1), ylim=c(0,1), main=paste0("md.loc at it=", t))
          points(Xorig[,d], pch=20, col="forestgreen")
          points(xold[d[1]], xold[d[2]], pch=20, col="red")
          points(xd[d[1]], xd[d[2]], pch=20, col="blue")
          rect(low[d[1]], low[d[2]], upp[d[1]], upp[d[2]], col=rgb(0,1,0,alpha=0.2))
          if(!is.null(Xorig)){
            legend("topright", c("xold", expression(x[d]), "Xorig"), xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col=c("red","blue", "forestgreen"), bty="n")
          }else legend("topright", c("xold", expression(x[d])), xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col=c("red", "blue"), bty="n")
        }
      }else{   ## after going through each location with 'md' as well as with 'non-md': start "jumping out of the well", i.e.,
        ## ask the stuck location to jump to a place close to the X location which has the maximum md to other X locations and Xorig. 
        ## Such a way guarantees the start location no longer to be stuck, and thus, there is no necessity to be iterative.
        
        ## justify 'xs'
        ## for each location in X, calculate its 'md' to other X locations and Xorig
        me <- apply(DD2, 1, function(x) min(x[x>0]))
        ## justify the index with the maximum md 
        max.ind <- which(rowSums(DD2 < max(me))==1)
        inds <- ceiling(runif(1)*length(max.ind))
        max.inds <- max.ind[inds]
        xs <- X[max.inds,]
        
        ## justify 'xd'
        Dor <- order(DD2[max.inds,], decreasing=FALSE)
        xd <- ifelse(rep(Dor[2] <= n, p), X[Dor[2],], Xorig[Dor[2]-n,])
        
        ## justify the search space
        sw.in <- sp(xs=xs, xd=xd, md=NULL)
        low <- sw.in$lo
        upp <- sw.in$up
        
        ## the new start location
        st.new <- st(xs=xs, xd=xd, ds=sw.in$ds)
        
        ## update X
        X[row.out,] <- st.new ## use the "new" location to replace the stuck location
        
        start <- matrix(X[row.out,], nrow=1)
        out <- optim(par=start, fn=maximin.obj, method="L-BFGS-B", 
                     lower=low, upper=upp, control=list(fnscale=-1), 
                     X=X, row.out=row.out, Xorig=Xorig)
        
        ## plot 'xs' and 'xd'
        if(plot==TRUE){
          Xv <- X[,d]
          plot(Xv, xlab=paste0("x", d[1]), ylab=paste0("x", d[2]), xlim=c(0,1), ylim=c(0,1), main=paste0("md.loc at it=", t))
          points(Xorig[,d], pch=20, col="forestgreen")
          points(xold[[1]], xold[d[2]], pch=20, col="cyan")
          points(xs[d[1]], xs[d[2]], pch=20, col="blue")
          points(st.new[d[1]], st.new[d[2]], pch=20, col="red")
          points(xd[d[1]], xd[d[2]], pch=20, col="orange")
          rect(low[d[1]], low[d[2]], upp[d[1]], upp[d[2]], col=rgb(0,1,0,alpha=0.2))
          arrows(xold[d[[1]]], xold[d[[2]]], st.new[d[[1]]], st.new[d[[2]]], lty=2, col="gray")
          if(!is.null(Xorig)){
            legend("topright", c("st.new", "Xorig"), xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col=c("red","forestgreen"), bty="n")
          }else legend("topright", "st.new", xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col="red", bty="n")
        }
      }
    }
    
    if(plot==TRUE){
      xnew <- as.numeric(out$par)
      Xv <- X[,d]
      plot(Xv, xlab=paste0("x", d[1]), ylab=paste0("x", d[2]), xlim=c(0,1), ylim=c(0,1), main=paste0("md.loc at it=", t))
      points(Xorig[,d], pch=20, col="forestgreen")
      points(X[row.out,][d[1]], X[row.out,][d[2]], pch=20, col="red")
      points(xnew[d[1]], xnew[d[2]], pch=20, col="blue")
      rect(low[d[1]], low[d[2]], upp[d[1]], upp[d[2]], col=rgb(0,1,0,alpha=0.2))
      arrows(X[row.out,][d[1]], X[row.out,][d[2]], xnew[d[1]], xnew[d[2]], length=0.1, lty=2, col="black") ## 'start' --> 'xnew'
      if(!is.null(Xorig)){
        legend("topright", c("xnew", "Xorig"), xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col=c("blue","forestgreen"), bty="n")
      }else legend("topright", "xnew", xpd=TRUE, horiz=TRUE, inset=c(-0.015, -0.045), pch=20, col="blue", bty="n")
    }
    
    X[row.out,] <- out$par
    D[row.out, -row.out] <- D[-row.out, row.out] <- as.numeric(distance(matrix(X[row.out,], nrow=1), X[-row.out,])) ## update D
    md <- min(as.numeric(D[upper.tri(D)]))
    md.ind <- which(D==md, arr.ind=TRUE)[,1]
    
    ## distances to the existing design
    if(!is.null(Xorig)){
      D2[row.out,] <- as.numeric(distance(matrix(X[row.out,], nrow=1), Xorig)) ## update D2
      md2 <- min(D2)
      if(md2 < md){
        md <- md2
        md.ind <- which(D2==md2, arr.ind=TRUE)[,1]
      }else if(md2==md) md.ind <- unique(c(md.ind, which(D2==md2, arr.ind=TRUE)[,1]))
    }
    
    if(plot==TRUE) mtext(paste("md =", format(round(as.numeric(md), 5), nsmall=5)))
    
    mind[t+1] <- md    ## track 'md' at each iteration
    Xind[[t+1]] <- X   ## track 'X' at each iteration
    
    ## When t >= 1.3 * n, if there is no improvement on 'md' during the previous 
    ## (0.3 * n) iterations, then I shall push the algorithm to make progress ...
    ## In addition, I shall not do any push during the last (0.2 * T - 1) iterations 
    ## since such push is only effective in the relatively long term ...
    if((t >= 1.3 * n) && (t <= 0.8 * T) && (cvs(mind[(t+1-0.3*n):(t+1)]) == TRUE)){
      
      if(is.null(Xorig)==TRUE) D2 <- NULL
      DD2 <- cbind(D, D2)
      
      ## justify 'xs'
      me <- apply(DD2, 1, function(x) min(x[x>0]))
      max.ind <- which(rowSums(DD2 < max(me))==1)
      inds <- ceiling(runif(1)*length(max.ind))
      max.inds <- max.ind[inds]
      xs <- X[max.inds,]
      
      ## justify 'xd'
      Dor <- order(DD2[max.inds,], decreasing=FALSE)
      xd <- ifelse(rep(Dor[2] <= n, p), X[Dor[2],], Xorig[Dor[2]-n,])
      
      ## justify the search space
      sw.in <- sp(xs=xs, xd=xd, md=NULL)
      low <- sw.in$lo
      upp <- sw.in$up
      
      ## the new start location
      st.new <- st(xs=xs, xd=xd, ds=sw.in$ds)
      
      ## update X
      X[row.out,] <- st.new
      
      start <- matrix(X[row.out,], nrow=1)
      out <- optim(par=start, fn=maximin.obj, method="L-BFGS-B", 
                   lower=low, upper=upp, control=list(fnscale=-1), 
                   X=X, row.out=row.out, Xorig=Xorig)
      
      X[row.out,] <- out$par
      D[row.out, -row.out] <- D[-row.out, row.out] <- as.numeric(distance(matrix(X[row.out,], nrow=1), X[-row.out,])) ## update D
      md <- min(as.numeric(D[upper.tri(D)]))
      md.ind <- which(D==md, arr.ind=TRUE)[,1]
      
      ## distances to the existing design
      if(!is.null(Xorig)){
        D2[row.out,] <- as.numeric(distance(matrix(X[row.out,], nrow=1), Xorig)) ## update D2
        md2 <- min(D2)
        if(md2 < md){
          md <- md2
          md.ind <- which(D2==md2, arr.ind=TRUE)[,1]
        }else if(md2==md) md.ind <- unique(c(md.ind, which(D2==md2, arr.ind=TRUE)[,1]))
      }
      mind[t+1] <- md    ## track 'md' at each iteration
      Xind[[t+1]] <- X   ## track 'X' at each iteration
    }
    
    if((verb==TRUE) && (t%%10==0)) cat("t=", t, "/T=", T, " is done.\n", sep="")
  }
  
  ## at the last iteration, recover the maximum minimum distance as well as the corresponding X
  mind[t+1] <- max(mind)
  X <- Xind[[which.max(mind)]]
  
  XXorig <- rbind(Xorig, X)
  return(list(Xf=XXorig, mi=mind))
}
