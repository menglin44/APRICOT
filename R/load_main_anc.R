glm_pval<-function(y,x){
  GLM<-glm(y ~ x, family="binomial")
  return(as.numeric(summary(GLM)$coefficients[2,4]))
}

sample_pval<-function(pheno, case.size, control.size,  sim.ancestry){
  index.cases <- which(pheno==1)
  index.ctrl <- which(pheno==0)
  index.samples<- c(sample(index.cases, size=case.size, replace=F), sample(index.ctrl, size=control.size, replace=F))
  p.pops <- apply(sim.ancestry[index.samples,], MARGIN=2, FUN=glm_pval, y=pheno[index.samples])
  
  return(p.pops)
}


power_calc <- function(pvals){
  count<-length(which(pvals<=0.05))
  return(count/length(pvals))
}

logit2p<-function(logit){
  p=1/(1+1/exp(logit))
  return(p)
}

p2logit<-function(p){
  logit<- log(p/(1-p))
  return(logit)
}

bin_pheno <- function(sim.p, case.size, control.size){
  temp.b.pheno <-NA
  iter <-1
  while(length(which(temp.b.pheno==1)) < case.size | length(which(temp.b.pheno==0)) < control.size){
    if(iter > 500) stop("Unable to get enough samples for case/control. Consider increasing N.adm. ")
    temp.b.pheno <- rbinom(length(sim.p), 1, sim.p) # convert proxy to binary
    iter <- iter +1
  }
  return(temp.b.pheno)
}



proxy_pheno <- function(index.ancestry, anc.matrix, proxy.anc, h2){
  sim.proxy.pheno <- c()
  sim.ancestry <- anc.matrix[index.ancestry,]
  
  ## simulating admixed proxy phenotypes based on ancestry and ancestral proxy phenotypes ##
  # genetic / additive
  sim.proxy.gen <- colSums(proxy.anc*t(sim.ancestry)) # multiplied to each indiv's vector of ancestries; additive, no h2
  
  # add nurture noise, based on h2
  var_gen = var(sim.proxy.gen) # phenotypic var attributed to genetics
  var_env = var_gen / h2 * (1-h2) # phenotypic var attibuted to environment
  sig_env = sqrt(var_env) # sigma of normal distribution to model the environmental noise
  sim.proxy.pheno = sim.proxy.gen + rnorm(length(sim.proxy.gen), 0, sig_env)
  
  return(sim.proxy.pheno)
}


bin_pheno_liab <- function(sim.proxy.pheno, liability){
  cutoff <- quantile(sim.proxy.pheno, probs = 1-liability)
  sim.bin.pheno <- as.numeric(sim.proxy.pheno>=cutoff)
  
  return(sim.bin.pheno)
}