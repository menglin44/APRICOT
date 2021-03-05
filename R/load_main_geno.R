#global anc proprotions from beta distribution
GAncBeta <- function(mean, var, samplesize){
  alpha <- mean^2*((1-mean)/var - 1/mean)
  beta <- (1-mean)/mean*alpha
  ganc <- rbeta(samplesize, alpha, beta)
  return(ganc)
}

#set genotypes for an admixed person from binomial sampling based on freq
LAnc2Geno <- function(index,lanc_p,lanc_m, frq_pop1, frq_pop2){
  hap_p <- rep(NA, length(lanc_p[index,]))
  hap_m <- rep(NA, length(lanc_m[index,]))
  #paternal hap
  hap_p[lanc_p[index,]==1] <- rbinom(sum(lanc_p[index,]==1), 1, frq_pop1[lanc_p[index,]==1])
  hap_p[lanc_p[index,]==0] <- rbinom(sum(lanc_p[index,]==0), 1, frq_pop2[lanc_p[index,]==0])
  #maternal hap
  hap_m[lanc_m[index,]==1] <- rbinom(sum(lanc_m[index,]==1), 1, frq_pop1[lanc_m[index,]==1])
  hap_m[lanc_m[index,]==0] <- rbinom(sum(lanc_m[index,]==0), 1, frq_pop2[lanc_m[index,]==0])
  
  hap <- hap_p + hap_m
  return(hap)
}

# prs on causal variants
indiv_prs <- function(genos,weights=weights) return (sum(genos*weights))

# sampling freq from 1/x function, min(x) = 1/2N
# https://rpubs.com/a_pear_9/weird_distributions
rinverse <- function(n=1, minMAF=1e-4) { 
  u = runif(n) # u = CDF value
  return(exp(u)*minMAF)
}


# calc af
af <- function(genos) sum(genos)/(2*length(genos))

# calc power
powers <- function(p, alpha) return(length(which(p<=alpha))/ length(p))


# convert to 0, 1 based on liability
q2b <- function(quant, liability){
  cutoff <- quantile(quant, 1-liability)
  binary <- as.numeric(quant>=cutoff)
  return(binary)
}

# calc maf
maf <- function(genos) 0.5-abs(0.5-af(genos))

