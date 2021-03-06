% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linreg.r
\name{linreg}
\alias{linreg}
\alias{lr}
\title{Linear regression}
\usage{
linreg(formula, data)
}
\arguments{
\item{formula}{The formula object defining the right hand side 
and left hand side of your linear regression}

\item{data}{A data frame.}
}
\value{
An S3 class linreg object, which is a list that contanins:
\describe{
  \item{Coefficients}{A vector of estimated coefficients}
  \item{X}{Design matrix}
  \item{vcov}{Variance-covariance matrix}
  \item{formula}{The model formula}}
}
\description{
Solves linear regression with the normal equation
}
\details{
This fgunction takes a data frame and a formula objects
and creates a design matrix, inverts it with the \link{ginv}
function from the MASS package.

It returns an object from class linreg.
}
\references{
http://mathworld.wolfram.com/NormalEquation.html
}
\seealso{
\link{lm} \link{ginv} \link{formula}
}
\author{
Sur Herrera Paredes
}
