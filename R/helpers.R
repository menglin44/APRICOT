#' Ancestry Generator
#'
#' Helper function to generate 2-way admixture under beta distribution.
#' @param mean mean of ancestry
#' @param sd sd of ancestry
#' @param n number of individuals
#' @return A matrix of global ancestry percentages for n individuals
#' @examples
#' AncBetaOut(0.75, 0.01, 100)
#' @export

AncBetaOut <- function(mean, var, n){
  alpha <- mean^2*((1-mean)/var - 1/mean)
  beta <- (1-mean)/mean*alpha
  ganc1 <- rbeta(n, alpha, beta)
  ganc2 <- 1-ganc1
  anc_table <- matrix(c(ganc1, ganc2), nrow=n,byrow=F)
  return(anc_table)
}

