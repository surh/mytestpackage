#' Linear regression
#' 
#' Solves linear regression with the normal equation
#' 
#' This fgunction takes a data frame and a formula objects
#' and creates a design matrix, inverts it with the \link{ginv}
#' function from the MASS package.
#' 
#' It returns an object from class linreg.
#' 
#' @param formula The formula object defining the right hand side 
#' and left hand side of your linear regression
#' @param data A data frame.
#' 
#' @return An S3 class linreg object, which is a list that contanins:
#' \describe{
#'   \item{Coefficients}{A vector of estimated coefficients}
#'   \item{X}{Design matrix}
#'   \item{vcov}{Variance-covariance matrix}
#'   \item{formula}{The model formula}}
#'   
#' @references http://mathworld.wolfram.com/NormalEquation.html
#' 
#' @author Sur Herrera Paredes
#' 
#' @seealso \link{lm} \link{ginv} \link{formula}
#' 
#' @importFrom MASS ginv
#' 
#' @aliases lr
#' 
#' @export
linreg <- function(formula, data){
  # Get design matrix
  X <- model.matrix(formula,data = data)
  
  # Get dependent variable
  lhs <- all.vars(update(formula, . ~ 0))
  Y <- matrix(data[,lhs],ncol = 1)
  
  # Solve normal equation
  B_hat <- MASS::ginv(t(X) %*% X) %*% t(X) %*% Y
  
  # Get residuals
  Y_hat <- X %*% B_hat
  e <- Y - Y_hat
  sigma_2 <- sum((e - mean(e))^2) / (length(e) - 1)
  vcov <- sigma_2 * MASS::ginv(t(X)%*%X)
  
  # Format results
  B_hat <- as.vector(B_hat)
  names(B_hat) <- colnames(X)
  
  colnames(vcov) <- colnames(X)
  row.names(vcov) <- colnames(X)
  
  res <- list(Coefficients = B_hat,
              X = X,
              vcov = vcov,
              data = data,
              formula = formula)
  class(res) <- "linreg"
  
  return(res)
}



#' Print linreg object
#' 
#' Pretty display of linreg object
#' 
#' This is an S3 method for the generic \link{print} function.
#' 
#' It builds a matrix of coefficients and standard
#' error and it prints it
#' 
#' @param x Linreg object
#' @param digits Number of digits to print
#' 
#' @author Sur Herrera Paredes
#' 
#' @return Invisibly a matrix of coefficient estimates and
#' standard errors
#' 
#' @seealso \link{print} \link{linreg}
print.linreg <- function(x, digits = 2){
  # Build matrix
  res <- cbind(x$Coefficients, sqrt(diag(x$vcov)))
  row.names(res) <- names(x$Coefficients)
  colnames(res) <- c("Estimate", "SE")
  
  # Round
  res <- round(res, digits)
  
  # Print
  cat("This is a linreg model object:\n")
  printCoefmat(res)
}


#' Summarize linreg model fit
#' 
#' Provides a full summary of the lireg
#' fit and performs t-tests on the coefficients
#' 
#' This function is a n S3 method of the generic \link{summary}
#' function.
#' 
#' Takes a linreg object, calculates the degrees of freedom
#' and performs a two-tailed t-test on each coefficient where
#' the null hypohtesis is that its real value is zero.
#' 
#' @return Invisibly return a matrix of coefficient and
#' test results
#' 
#' @author Sur Herrera Paredes
#' 
#' @seealso \link{summary} \link{linreg}
#' 
#' @param object A linreg object
#' 
#' @export
summary.linreg <- function(object){
  # Calculate degrees of freedom
  df <- nrow(object$X) - ncol(object$X) + 1
  
  # Build matrix
  res <- cbind(object$Coefficients, sqrt(diag(object$vcov)))
  row.names(res) <- names(object$Coefficients)
  colnames(res) <- c("Estimate", "SE")
  
  # Calculate t-value
  res <- cbind(res, t.value = res[,1] / res[,2])
  
  # Calculate p-value
  res <- cbind(res,
               p.value = 2*(1 - pt(q = abs(res[,3]), df = df)))
  
  printCoefmat(res,P.values = TRUE,signif.stars = TRUE, signif.legend = TRUE, has.Pvalue = TRUE)
}

