
#make footnote for ggplot2
#png("mwba_gdp.png"); print(g); makeFootnote("xxx");dev.off()
makeFootnote <- function(footnoteText=
                         format(Sys.time(), "%d %b %Y"),
                         size= 1, color= grey(.1))
{
   require(grid)
   pushViewport(viewport())
   grid.text(label= footnoteText ,
             x = unit(1,"npc") - unit(2, "mm"),
             y= unit(2, "mm"),
             just=c("right", "bottom"),
             gp=gpar(cex= size, col=color))
   popViewport()
}

#handeling missing value using last varialbe
imp <- function (a){
  missing <- is.na(a)
  anymis <- any(missing)
  imputed <- a
  while(anymis){
    imputed[tail(which(missing==TRUE),1)] <- imputed[tail(which(missing==TRUE),1)+1]
    missing <- is.na(imputed)
    anymis <- any( missing )
  }
  return (imputed)
}

# transform list 2 df
list2df <- function(ll) {
  return(ldply(ll,function(l){ return(data.frame(rbind(unlist(l))))}))
}

#tauchen method for AR1
# =====================================================================
#tauchen <- function(rho,sigma,mu=0,n,m=3){
#    N=n
#    Z<- array(0, N)
#    Zprob <- matrix(0, nrow=N, ncol=N)
#    a     <- (1-rho)*mu
#
#    Z[N]  <- m * sqrt( sigma^2 / (1 - rho^2) )
#    Z[1]  <- -Z[N];
#    zstep <- (Z[N] - Z[1]) / (N - 1)
#    for ( i in 2:(N-1) ){
#        Z[i] <- Z[1] + zstep * (i - 1)
#    }    
#    
#    Z <- Z + a / (1-rho)
#
#    for (j in 1:N){
#        for (k in 1:N){
#            if (k == 1) {
#                Zprob[j,k] <-     pnorm( (Z[1] - a - rho * Z[j] + zstep / 2) / sigma )
#            }else if(k == N){
#                Zprob[j,k] <- 1 - pnorm( (Z[N] - a - rho * Z[j] - zstep / 2) / sigma )
#            }else{
#                Zprob[j,k] <- pnorm( (Z[k] - a - rho * Z[j] + zstep / 2) / sigma ) -
#                              pnorm( (Z[k] - a - rho * Z[j] - zstep / 2) / sigma )
#            }
#        }
#    }
#    return(list(Pmat=Zprob,zgrid=Z))
#}

tauchen <- function(rho,sigma,mu=0,n,m=3){
    N=n
    Z<- array(0, N)
    Zprob <- matrix(0, nrow=N, ncol=N)

    Z[N]  <- m * sqrt( sigma^2 / (1 - rho^2) )
    Z[1]  <- -Z[N];
    zstep <- (Z[N] - Z[1]) / (N - 1)
    for ( i in 2:(N-1) ){
        Z[i] <- Z[1] + zstep * (i - 1)
    }    
    
    Z <- Z + mu / (1-rho)

    for (j in 1:N){
        for (k in 1:N){
            if (k == 1) {
                Zprob[j,k] <-     pnorm( (Z[1] - mu - rho * Z[j] + zstep / 2) / sigma )
            }else if(k == N){
                Zprob[j,k] <- 1 - pnorm( (Z[N] - mu - rho * Z[j] - zstep / 2) / sigma )
            }else{
                Zprob[j,k] <- pnorm( (Z[k] - mu - rho * Z[j] + zstep / 2) / sigma ) -
                              pnorm( (Z[k] - mu - rho * Z[j] - zstep / 2) / sigma )
            }
        }
    }
    return(list(Pmat=Zprob,zgrid=Z))
}
#rw <- rouwenhorst(rho=rhow, sigma=stdw, mu=0, n=nw)
#rw <- tauchen(rho=rhow, sigma=stdw, mu=0, n=nw,m=3) 
#xw <- rw$zgrid # grid points for income distribution
#income <- exp(xw) #actual income
#wprob  <- rw$Pmat # markov transition matrix, row is today, column is tomorrow

# rouwenhorst discretization for AR1
# 
# translation of \url{http://www.karenkopecky.net/rouwenhorst.m}
# @references \url{http://www.karenkopecky.net/RouwenhorstPaperFinal.pdf}
# @param rho first order autocorrelation 
# @param sigma standard deviation of error term
# @param mu mean of error term
# @param n number of points to use in approximation
# @return list with Pmat (transition matrix) and zgrid (grid points)
# @export
# @examples
# R <- rouwenhorst(rho=0.9,sigma=1.1,mu=0,n=5)
# print(R$zgrid) # support points
# print(R$Pmat)  # transition matrix
# print(rowSums(R$Pmat))
rouwenhorst <- function(rho,sigma,mu=0,n){
  stopifnot(n>1)
  qu <- (rho+1)/2
  nu <- ((n-1)/(1-rho^2))^(1/2) * sigma
  P  <- matrix(c(qu,1-qu,1-qu,qu),nrow=2,ncol=2)
  if (n>2){
    for (i in 2:(n-1)){
      zeros    <- rep(0,i)
      zzeros   <- rep(0,i+1)
      P        <- qu * rbind(cbind(P,zeros,deparse.level=0),zzeros,deparse.level=0) + 
              (1-qu) * rbind(cbind(zeros,P,deparse.level=0),zzeros,deparse.level=0) + 
              (1-qu) * rbind(zzeros,cbind(P,zeros,deparse.level=0),deparse.level=0) + 
                  qu * rbind(zzeros,cbind(zeros,P,deparse.level=0),deparse.level=0)
      P[2:i, ] <- P[2:i, ]/2
    }
  }
  zgrid <- seq(from=mu/(1-rho)-nu,to=mu/(1-rho)+nu,length=n)
  return(list(Pmat=P,zgrid=zgrid))
}

# -----------------------------------
# Compute Transition Matrix
# -----------------------------------
#require(copula)
#computeCopula <- function(rho, grid, nn) {
#  if (rho < 1) {
#    norm.cop  <- normalCopula(rho)
#    vals_grid <- data.matrix( expand.grid(z1=grid, z2=grid) )   
#    vals      <- dCopula(vals_grid, norm.cop)
#    QQ        <- array( vals, dim=c(nn, nn) )
#    QQ        <- QQ / array( rowSums(QQ), dim=c(nn, nn) )
#  }
#  else { QQ <- diag(nn) }
#}

#######################################
############## function:
  merge.with.order <- function(x,y, ..., sort = T, keep_order)
  {
    # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
    add.id.column.to.data <- function(DATA)
    {
      data.frame(DATA, id... = seq_len(nrow(DATA)))
    }
    # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
    order.by.id...and.remove.it <- function(DATA)
    {
      # gets in a data.frame with the "id..." column.  Orders by it and returns it
      if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
 
      ss_r <- order(DATA$id...)
      ss_c <- colnames(DATA) != "id..."
      DATA[ss_r, ss_c]
    }
 
    # tmp <- function(x) x==1; 1  # why we must check what to do if it is missing or not...
    # tmp()
 
    if(!missing(keep_order))
    {
      if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
      if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
      # if you didn't get "return" by now - issue a warning.
      warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
    } else {return(merge(x=x,y=y,..., sort = sort))}
  }

  