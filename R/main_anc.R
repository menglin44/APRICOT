#' Power Estimate for Trait Ancestry Associations
#'
#' Estimates power for significant associations between a binary trait and global ancestries.
#' @param incidence.rates A vector of incidence rates for the ancestral populations.
#' @param N.adm Size of simulated individuals from which cases and controls are drawn from. Default is 3e5.
#' @param anc.matrix A sample matrix of global ancestries, where rows represent samples and each column is one ancestry. See the format from the output of \code{AncBetaOut} function. If an ancestry matrix from empirical estimates is not available, \code{AncBetaOut} can be used to supply ancestry draws from simulations.
#' @param n Iteration number of phenotype-ancestry regression tests. Each iteration is an independent simulation of phenotypes and ancestries for case and controls. Default 1000.
#' @param case.size Number of cases. Default 750.
#' @param control.size Number of controls. Default 750.
#' @param h2 Heritability. Default 1.
#' @param liability Liability of the binary trait, used for defining cases and controls. If left as NA (default), then a proxy phenotype based on ancestry and incidence rates is used to reflect the probability of becoming cases.
#' @param anc.names A vector of names for each ancestry, e.g. South Asian, West African, etc. It can be left unspecified (NA as default).
#' @return A list of the following components:\tabular{ll}{
#' \code{ancestry} \tab Name of ancestries\cr
#' \tab \cr
#' \code{liability_model} \tab If a liability treshold model has been used\cr
#' \tab \cr
#' \code{powers} \tab Power per ancestry\cr
#' \tab \cr}
#' @examples
#' # sample ancestries for a size of 100
#' aa_anc <- AncBetaOut(0.75, 0.01, 100)\cr
#' # simulate phenotype-ancestry associations
#' aa_ancp <- powphenanc(incidence.rates=c(0.01, 0.05), anc.matrix=aa_anc, h2=0.5, anc.names=c("eur","waf"))\cr
#' # check powers
#' anc_powers <- aa_ancp$powers\cr
#' @seealso \code{AncBetaOut}
#' @export

# this is the main function for phenotype-ancestry association power
powphenanc <- function(incidence.rates,
                       N.adm = 3e5,
                       anc.matrix,
                       n=1000,
                       case.size=750,
                       control.size=750,
                       h2=1,
                       liability=NA,
                       anc.names=NA) { # the main function

  #start_time <- proc.time()


  ##@ Checking errors @##
  n.pop <- try(length(incidence.rates),silent=T)
  anc.total <- try(rowSums(anc.matrix))
  if(!exists("n.pop")) stop("Please input incidence rates for ancestral populations.")
  if(!exists("anc.total")) stop("Ancestry input is expected to be a matrix, where each column represents one ancestry.")
  if(n.pop != ncol(anc.matrix)) stop("Number of populations are not equal in incidence.rates and anc.matrix." )
  if(length(which(anc.total< 0.95))>0) stop("Some individuals' ancestries don't sum to ~1." )
  if(n.pop==1) stop("There is only one ancestry component." )
  if(!is.na(liability) & !(liability>0 & liability<1)) stop("Liability is between 0 and 1.")
  rm(anc.total)
  ##@ Checking errors ends @##

  ##@ Spitting out parameters to screen and a log file @##
  print(paste("Admixed individuals are: ",n.pop, "-way admixed.", sep=""))
  print(paste("Proposed sample size for case / control:", case.size, control.size))
  print(paste("Size of the admixed population to be simulated:", N.adm))
  print(paste("Number of correlations to assess power:", n))
  ##@ Spitting out parameters ends @##

  ##@ reading in and prepare @##
  anc.matrix <- as.matrix(anc.matrix)


  #initiation of tables for records
  pval.table <- matrix(NA, nrow=n, ncol=n.pop)
  powers <- matrix(NA, nrow=1, ncol=n.pop)

  #converting incidence rates to proxy phenotypes for ancestral individuals
  proxy.anc <- sapply(incidence.rates, FUN=p2logit)

  #naming columns to populations#
  if(is.na(anc.names[1])) {
    colnames(powers)<-paste("Ancestry", c(1:n.pop))
  } else {
    colnames(powers)<-anc.names
  }


  #iterations of associations
  pb <- progress::progres_bar$new(format="Association iterations [:bar] :percent", total=n, clear=F, width=60)

  for(i in c(1:n)){

    index.ancdraw <- sample(c(1:dim(anc.matrix)[1]), size=N.adm, replace=T)
    sim.proxy.pheno <- proxy_pheno(index.ancestry=index.ancdraw, anc.matrix=anc.matrix, proxy.anc=proxy.anc, h2=h2)

    #convert logit to binary phenotypes under two alternative methods
    ## Method 1. If liability threshold is given
    if(!is.na(liability) & liability >0 & liability <1){
      #make sure w/ liability, the pool of N.adm is big enough to get enough cases
      if((1-liability)*N.adm < case.size | liability*N.adm < control.size){
        stop("Unable to get enough cases/controls given the liability. Consider increasing N.adm.")
      }

      sim.bin.pheno <- bin_pheno_liab(sim.proxy.pheno, liability)

      ## Method 2. Use logit pheno to draw binomial disease outcome, when liability not provided
    } else if(is.na(liability)){
      sim.p <- logit2p(sim.proxy.pheno)
      sim.bin.pheno <- bin_pheno(sim.p=sim.p, case.size=case.size, control.size=control.size)
    }

    pval <- sample_pval(pheno=sim.bin.pheno, case.size=case.size, control.size=control.size,  sim.ancestry=anc.matrix[index.ancdraw,])
    pval.table[i,] <- pval

    pb$tick()
  }

  powers[1,] <- apply(pval.table, MARGIN=2, FUN=power_calc)

  ##@ return list @##
  populations = colnames(powers)
  liability_model = c(T, F)[as.numeric(is.na(liability)) + 1]
  powers = powers
  return(list(ancestry=populations, liability_model=liability_model, powers=powers))
}




