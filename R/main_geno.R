#' Power Inference for GWAS Discovery in Admixed Populations
#'
#' Infers statistical power for discovery of variants in a 2-way admixed group, and the corresponding populations of ancestry origin.
#' @param n_snps Total number of SNPs (causal and non-causal). SNPs are assumed independent from each other. Default is 1000.
#' @param popsizes A vector of sample sizes for the admixed group, Pop 1 (ancestral origin from which the major component of ancestry comes), Pop 2 (ancestral origin from which the minor component of ancestry comes). Default is 1000 each.
#' @param mean_ganc Mean global ancestry from Pop 1 in the admixed group. Default is 0.75.
#' @param sd_ganc Standard deviation of the global ancestry distribution from Pop 1 in the admixed group. Default is 0.1.
#' @param prop_causal Proportion of causal variant among total SNPs. Default is 0.1.
#' @param hg2 Heritability of the trait. Default is 0.5.
#' @param fst The baseline fixation index FST among causal variants between Pop 1 and Pop2. Default is 0.2.
#' @param fst_constant True or False: Is FST among causal variants a constant value equaling the baseline (True) or have an increment with the rarity of ancestral minor allele frequency (False)? Default is False.
#' @param freqdist Frequency distribution of ancestral variants (prior to the divergence of Pop 1 and 2). "unif" for a uniform distribution modeling (default), "inverse" for an inverse distribution modeling.
#' @param liability A vector of liability threshold for dichotomous traits in the admixed group, Pop 1, Pop2, respectively. If set as NA (default), the trait is assumed to be quantitative.
#' @param aaf Limits (min, max) on the frequency distribution of ancestral variants. Default is c(0.001, 0.999).
#' @param envanc Options to model environment by ancestry interactions: 0 means no interaction. 1 means interation is modeled as the sum of Gaussian sampling. 2 means the interaction is modeled as linear dependence on the ancestry proportion. Default is 1.
#' @param envvar percentage of the environment by ancestry interaction over the total phenotypic variance. Only effective when \code{envanc} is 2.
#' @param alpha Type I error rate. Default is 0.05.
#' @return A list of powers and the re-estimated heritability from simulation.\tabular{ll}{
#' \code{admixed_power} \tab Power for discovery in the admixed group\cr
#' \tab \cr
#' \code{pop1_power} \tab Power for discovery in Pop 1\cr
#' \tab \cr
#' \code{pop2_power} \tab Power for discovery in Pop 2\cr
#' \tab \cr
#' \code{h2} \tab Estimated heritability from simulated phenotypes and genotypes \cr
#' \tab \cr}
#'
#' @export
# this is the main function for genotype mediated simulation
SimGenoPower <-function(n_snps=1000,
                        popsizes=c(1e3, 1e3, 1e3),
                        mean_ganc=0.75,
                        sd_ganc=0.1,
                        prop_causal=0.1,
                        hg2=0.5,
                        fst=0.2,
                        fst_constant=F,
                        freqdist="unif",
                        liability=NA,
                        aaf=c(1e-3, 0.999),
                        envanc=1,
                        envvar=0,
                        alpha=0.05){

  ptm <- proc.time()

  ###########################
  #       check input       #
  ###########################

  n.pop <- try(length(popsizes),silent=T)
  if(n.pop != 3) stop("Need sample size of two ancestral populations and the admixed population." )
  if (freqdist!="unif" & freqdist!="inverse") stop("Ancestral frequencies are drawn from either uniform distribution or inverse distribution.")
  if (!is.na(liability[1]) & sum(liability>0 & liability<1)!=3) stop("Liability for 3 populations should be a vector of 3 values (adx, pop1, pop2) between 0 and 1.
                                                                   Leave it as NA if for quantitative traits.")
  if (!envanc%in%c(0,1,2)) stop("Specify if environment by ancestry is modeled. 0: no. 1: gaussian. 2: linear.")
  if (envanc==2 & !(envvar>=0 & envvar<1)) stop("Specify the proportion of phenotypic variance explained by environment by ancestry (linear).")
  ######################################
  #     draw ancestry and genotypes    #
  ######################################
  print("Simulating ancestry and genotypes.")

  pop1_idvs <- popsizes[1]
  pop2_idvs <- popsizes[2]
  adx_idvs <- popsizes[3]

  #global
  ganc <- GAncBeta(mean_ganc, sd_ganc^2, adx_idvs)

  #local ancestry per marker: 0/1/2 copies of pop1 ancestry
  lanc_p <- matrix(rbinom(adx_idvs*n_snps, 1, ganc), nrow=adx_idvs, ncol=n_snps)
  lanc_m <- matrix(rbinom(adx_idvs*n_snps, 1, ganc), nrow=adx_idvs, ncol=n_snps)
  lanc_pop1 <- lanc_p + lanc_m

  #set freq per ancestral population at each locus, give fst and human ancestral allele maf
  pop1.frqs <- rep(NA, n_snps)
  pop2.frqs <- rep(NA, n_snps)
  for (s in 1:n_snps) {
    # set allele frequencies for ancestries
    p.pop1 <- 0
    p.pop2 <- 0
    if(freqdist=="unif"){
      anc.frq <- runif(1, min=aaf[1], max=aaf[2])
    }else if(freqdist=="inverse"){
      anc.frq <- rinverse(1, minMAF=1e-3)
    }

    anc.maf <- 0.5 - abs(0.5-anc.frq)

    if(fst_constant){
      Fst <- fst
    }else{
      Fst <- fst + (1-anc.maf)*0.3 # the smaller maf, the bigger fst
    }

    #Balding-Nichol model
    while((p.pop1 <= 0 || p.pop1 >= 1) || p.pop2 <= 0 || p.pop2 >= 1){
      p.pop1 <- rbeta(1, anc.frq*(1-Fst)/Fst,(1-anc.frq)*(1-Fst)/Fst)
      p.pop2 <- rbeta(1, anc.frq*(1-Fst)/Fst,(1-anc.frq)*(1-Fst)/Fst)
    }
    pop1.frqs[s] <- p.pop1
    pop2.frqs[s] <- p.pop2
  }


  # set genotypes based on frequency
  #1. pop1 and pop2
  pop1_genos <- matrix(rbinom(n_snps*pop1_idvs, 1, pop1.frqs), byrow=T, nrow=pop1_idvs, ncol=n_snps) +
    matrix(rbinom(n_snps*pop1_idvs, 1, pop1.frqs), byrow=T, nrow=adx_idvs, ncol=n_snps)

  pop2_genos <- matrix(rbinom(n_snps*pop2_idvs, 1, pop2.frqs), byrow=T, nrow=pop2_idvs, ncol=n_snps) +
    matrix(rbinom(n_snps*pop2_idvs, 1, pop2.frqs), byrow=T, nrow=pop2_idvs, ncol=n_snps)

  #2. admixed population
  adx_genos <- t(sapply(c(1:adx_idvs), FUN=LAnc2Geno, lanc_p=lanc_p,lanc_m=lanc_m, frq_pop1=pop1.frqs, frq_pop2=pop2.frqs))
  adx.frqs <- apply(adx_genos, MARGIN=2, FUN=af)

  ##############################################
  #  simulating phenotype based on prs and h2  #
  ##############################################
  print("Simulating phenotypes.")

  # effect size is tied to fst+ancestral maf
  # only causal variants are assigned with weights; assuming trait is inverse normalized
  #weights <- rbinom(n_snps,1,prop_causal) * rnorm(n_snps, 0, 1) * ((rbinom(n_snps, 1, pop1.frqs/(pop1.frqs+pop2.frqs))-0.5)*2)
  causal_index <- sample(c(1:n_snps), size=round(n_snps*prop_causal), replace=F)
  if_causal <- rep(0, n_snps)
  if_causal[causal_index]<-1
  weights <- if_causal * rnorm(n_snps, 0, 1) * ((rbinom(n_snps, 1, pop1.frqs/(pop1.frqs+pop2.frqs))-0.5)*2)

  adx_prs <- apply(adx_genos,MARGIN=1,FUN=indiv_prs,weights=weights)
  pop1_prs <- apply(pop1_genos, MARGIN=1,FUN=indiv_prs,weights=weights)
  pop2_prs <- apply(pop2_genos, MARGIN=1,FUN=indiv_prs,weights=weights)

  # add environmental noise # including gaussian and ancestry associated
  pop1_phenos <- pop1_prs + rnorm(pop1_idvs, 0, 1) * sqrt(var(pop1_prs)*(1/hg2-1))
  pop2_phenos <- pop2_prs + rnorm(pop2_idvs, 0, 1) * sqrt(var(pop2_prs)*(1/hg2-1))
  if(envanc==1){# including gaussian and ancestry associated
    # make sure the h2 is not changed
    #scaling <- sqrt(var(adx_prs) / (var(adx_prs) + var(pop1_prs)*ganc^2 + var(pop2_prs)*((1-ganc)^2)))
    m_ganc <- mean(ganc)
    v_ganc <- var(ganc)
    scaling <- sqrt(var(adx_prs)/(var(adx_prs) + var(pop1_prs)*(m_ganc^2 + v_ganc) + var(pop2_prs)*(1-2*m_ganc + m_ganc^2 + v_ganc)))
    adx_phenos <- adx_prs +
      (rnorm(adx_idvs, 0, 1) * sqrt(var(adx_prs)*(1/hg2-1)) +
       rnorm(pop1_idvs, 0, 1) * sqrt(var(pop1_prs)*(1/hg2-1))*ganc +
       rnorm(pop2_idvs, 0, 1) * sqrt(var(pop2_prs)*(1/hg2-1))*(1-ganc))*scaling
  }else if(envanc==0){
    adx_phenos <- adx_prs + rnorm(adx_idvs, 0, 1) * sqrt(var(adx_prs)*(1/hg2-1))
  }else if(envanc==2){
    scaling <- sqrt((1-(hg2+envvar))/(1-hg2))
    if(scaling==0){# when p = 1-hg2; the original equation doesn't stand
      beta_theta <- sqrt(var(adx_prs)*(1-hg2)/(var(ganc)*hg2))
      adx_phenos <- adx_prs + beta_theta*ganc
    }else{
      beta_theta <- sqrt((envvar*(1-hg2)*var(adx_prs))/(hg2*(1-(hg2+envvar))*var(ganc)))
      adx_phenos <- adx_prs + scaling * (rnorm(adx_idvs, 0, 1) * sqrt(var(adx_prs)*(1/hg2-1)) +
                                           beta_theta*ganc)
    }
  }
  h2 <- var(adx_prs)/var(adx_phenos)

  # if the trait is dichotomous, converting from
  if(!is.na(liability[1])){
    adx_phenos <- q2b(adx_phenos, liability[1])
    pop1_phenos <- q2b(pop1_phenos, liability[2])
    pop2_phenos <- q2b(pop2_phenos, liability[3])
  }


  #############################################
  #  perform correlations on causal variants  #
  #############################################

  p_adx <-c()
  p_1<-c()
  p_2<-c()
  index_causal <- which(weights!=0)

  pb<- progress::progress_bar$new(format="Performing associations [:bar] :percent", total=length(index_causal),clear=F, width=80)
  for(i in index_causal){

    if(!is.na(liability[1])){#dichotomous

      if(var(adx_genos[,i])!=0){
        glm_temp <- summary(glm(adx_phenos~ganc + adx_genos[,i], family="binomial"))
        p_adx <- c(p_adx, glm_temp$coefficients[3,4])
      }else{
          p_adx <- c(p_adx, 1)
      }


      if(var(pop1_genos[,i])!=0){
          glm_temp <-summary(glm(pop1_phenos~ pop1_genos[,i], family="binomial"))
          p_1 <- c(p_1, glm_temp$coefficients[2, 4])
      }else{
          p_1 <- c(p_1, 1)
      }


      if(var(pop2_genos[,i])!=0){
        glm_temp <-summary(glm(pop2_phenos~ pop2_genos[,i], family="binomial"))
      }else{
        p_2 <- c(p_2, 1)
      }

    }else{# quantitative
      if(var(adx_genos[,i])==0){
        p_adx <- c(p_adx, 1)
      }else{
        lm_temp <- anova(lm(adx_phenos~ganc + adx_genos[,i]))
        p_adx <- c(p_adx, lm_temp$`Pr(>F)`[2])
      }

      if(var(pop1_genos[,i])==0){
        p_1<- c(p_1,1)
      }else{
        lm_temp <-anova(lm(pop1_phenos~ pop1_genos[,i]))
        p_1 <- c(p_1, lm_temp$`Pr(>F)`[1])
      }

      if(var(pop2_genos[,i])==0){
        p_2<- c(p_2,1)
      }else{
        lm_temp <-anova(lm(pop2_phenos~ pop2_genos[,i]))
        p_2 <- c(p_2, lm_temp$`Pr(>F)`[1])
      }
    }
    pb$tick()
  }

  admixed_power <- powers(p_adx, alpha)
  pop1_power <- powers(p_1, alpha)
  pop2_power <- powers(p_2, alpha)

  print(paste("Elapsed time: ", as.numeric((proc.time() - ptm)[3]), "s", sep=""))
  return(list(admixed_power=admixed_power, pop1_power=pop1_power, pop2_power=pop2_power, h2_adx=h2))

}
