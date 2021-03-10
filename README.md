# APRICOT
**A**dmixed **P**opulation powe**R** **I**nference **C**omputed for phen**O**typic **T**raits

This package provides a mechanistic simulation framework to infer statistical power for discovery of causal variants in an admixed population and the correpsonding ancestral origin populations. It sets in a 2-way admixture scenario with assumption of independent markers, and allows for flexibility in adjusting environmental effect, trait divergence, polygenicity, heritability etc. *APRICOT* also provides a function to estimate power for phenotype-ancestry association for dichotomous traits based on incidence rates and ancestry distributions. See [*Lin et al. bioRxiv, 2021*](https://www.biorxiv.org/content/10.1101/2021.03.09.434643v1) for more details on the model.

## Setup

In R, you can install *APRICOT* with the following commands:
```
devtools::install_github("menglin44/APRICOT")
library(APRICOT)
```

## Functions and Usage

- Genotype-mediated simulation framework to test for GWAS power in an admixed group
  - See ```?SimGenoPower```
- Power for phenotype-ancestry associations
  - See ```?powphenanc```
  - This function requires an example ancestry matrix as part of input. If no emprical estimates are available, one can generate an approximate ancestry under a beta distribtuion. A helper function ```AncBetaOut``` is provided.

## Reference

[*Lin et al. bioRxiv, 2021*](https://www.biorxiv.org/content/10.1101/2021.03.09.434643v1) 
