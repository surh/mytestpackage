#' Perform a permutation test on a linreg object
#' 
#' Permutes the observations and calculates a p-value
#' for every coefficient based on those permutations.
#' 
#' Takes a linreg object, permutes the observations
#' and uses the estimated null distributions to obtain
#' a two-sided p-value of each coefficient in the
#' linreg model.
#' 
#' It ignores cases of undefined values in permutations.
#' 
#' @param model A linreg model
#' @param nperm Number of permutations
#' 
#' @author Sur Herrera Paredes
#' 
#' @seealso \link{linreg}
#' 
#' @return Silently returns a list that includes the
#' following elements:
#' \describe{
#'   \item{perms}{A matrix containing the values for all coefficients
#'   on each permutation}
#'   \item{coefficients}{A matrix of coefficients and their permutation
#'   based p-values}}
#'  
#' @export
permtest <- function(model, nperm = 99){
  if(class(model) != "linreg")
    stop("ERROR: model must be of class linreg", call. = TRUE)
  
  # Get dependent variable
  lhs <- all.vars(update(model$formula, . ~ 0))
  
  # Perform permutations
  perms <- matrix(model$Coefficients, ncol = 1)
  row.names(perms) <- names(model$Coefficients)
  for(i in 1:nperm){
    
    # Permute data
    dat_i <- model$data
    dat_i[,lhs] <- sample(dat_i[,lhs],size = nrow(dat_i), replace = FALSE)
    
    # Fit permuted data
    mi <- linreg(formula = model$formula, data = dat_i)
    perms <- cbind(perms, mi$Coefficients)
  }
  
  # Calculate permutation statistics
  res <- perms[ , 1, drop = FALSE]
  colnames(res) <- "Estimate"
  perm.quantiles <- t(apply(perms, 1, quantile, probs = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  res <- cbind(res, perm.quantiles)
  
  # Calculate permutation p-value
  pvals <- apply(perms, 1, function(x){
    sum(abs(x) >= abs(x[1])) / length(x)
  })
  res <- cbind(res, p.value = pvals )
  
  # Print coefficient matrix
  printCoefmat(res,P.values = TRUE,signif.stars = TRUE, signif.legend = TRUE, has.Pvalue = TRUE)
  
  res <- list(perms = perms, coefficients = res)
  invisible(res)
}