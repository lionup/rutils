\name{rouwenhorst}
\alias{rouwenhorst}
\title{rouwenhorst discretization for AR1}
\usage{
  rouwenhorst(rho, sigma, mu = 0, n)
}
\arguments{
  \item{rho}{first order autocorrelation}

  \item{sigma}{standard deviation of error term}

  \item{mu}{mean of error term}

  \item{n}{number of points to use in approximation}
}
\value{
  list with Pmat (transition matrix) and zgrid (grid
  points)
}
\description{
  translation of
  \url{http://www.karenkopecky.net/rouwenhorst.m}
}
\examples{
R <- rouwenhorst(rho=0.9,sigma=1.1,mu=0,n=5)
print(R$zgrid)	# support points
print(R$Pmat)	# transition matrix
print(rowSums(R$Pmat))
}
\references{
  \url{http://www.karenkopecky.net/RouwenhorstPaperFinal.pdf}
}

