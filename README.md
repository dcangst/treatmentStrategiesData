[![DOI](https://zenodo.org/badge/257538463.svg)](https://zenodo.org/badge/latestdoi/257538463)
# Treatment Strategies Data Files & Scripts

This repository contains all experimental data for the paper LINK AND CITE ONCE PUBLISHED as well as the scripts used for analysis.

# Description of files

## Folder scripts

Contains scripts used for data analysis, plot creation and statistics. All plots were generated in R 3.6.2 in April 2020.

`R sessionInfo()`

```R
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] kableExtra_1.2.1   knitr_1.30         multcompView_0.1-8 multcomp_1.4-14   
 [5] TH.data_1.0-10     MASS_7.3-51.6      survival_3.1-12    mvtnorm_1.1-1     
 [9] forcats_0.5.0      stringr_1.4.0      dplyr_1.0.2        purrr_0.3.4       
[13] readr_1.3.1        tidyr_1.1.2        tibble_3.0.3       ggplot2_3.3.2     
[17] tidyverse_1.3.0    here_0.1          

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5        lubridate_1.7.9   lattice_0.20-41   zoo_1.8-8         utf8_1.1.4       
 [6] assertthat_0.2.1  rprojroot_1.3-2   digest_0.6.25     R6_2.4.1          cellranger_1.1.0 
[11] backports_1.1.10  reprex_0.3.0      evaluate_0.14     httr_1.4.2        pillar_1.4.6     
[16] rlang_0.4.7       readxl_1.3.1      rstudioapi_0.11   blob_1.2.1        Matrix_1.2-18    
[21] rmarkdown_2.3     labeling_0.3      splines_4.0.2     webshot_0.5.2     munsell_0.5.0    
[26] broom_0.7.0       compiler_4.0.2    modelr_0.1.8      xfun_0.17         pkgconfig_2.0.3  
[31] htmltools_0.5.0   tidyselect_1.1.0  codetools_0.2-16  fansi_0.4.1       viridisLite_0.3.0
[36] crayon_1.3.4      dbplyr_1.4.4      withr_2.3.0       grid_4.0.2        jsonlite_1.7.1   
[41] gtable_0.3.0      lifecycle_0.2.0   DBI_1.1.0         magrittr_1.5      scales_1.1.1     
[46] cli_2.0.2         stringi_1.5.3     farver_2.0.3      fs_1.5.0          xml2_1.3.2       
[51] ellipsis_0.3.1    generics_0.0.2    vctrs_0.3.4       sandwich_3.0-0    tools_4.0.2      
[56] glue_1.4.2        hms_0.5.3         colorspace_1.4-1  rvest_0.3.6       haven_2.3.1
```

## 1-phenotypeDefinition.csv

Phenotypes were defined using the following set of rules:

- phenotype U (uninfected): no growth, empty well
- phenotype WT: growth only on agarN and in liquid without treatment
- phenotype A: growth in liquid with A or on plate with A, not in other conditions
- phenotype B: growth in liquid with B or on plate with B, not in other conditions
- phenotype A/B, a mixed culture of A and B, no growth on AB but growth in A and B (liquid or solid media)
- phenotype AB, double Resistance, growth on A and B (liquid or solid media)
- E: "error", i.e. not any of the above.

### Columns (`col_types = "llllllllic"`)
  - **agarN**: boolean, visible growth on MS Agar containing no drugs.
  - **agarA**: boolean, visible growth on MS Agar containing 40µg/ml nalidixic acid.
  - **agarB**: boolean, visible growth on MS Agar containing 100µg/ml streptomycin.
  - **agarAB**: boolean 40µg/ml nalidixic acid (agarA) and 100µg/ml streptomycin (agarAB)
  - **growth**: boolean, OD(595nm) > 0.1 after incubation in any treatment.
  - **growthA**: boolean, OD(595nm) > 0.1 after incubation in media containing 20µg/ml nalidixic acid.
  - **growthB**: boolean, OD(595nm) > 0.1 after incubation in media containing 12.5µg/ml streptomycin.
  - **growthAB**:  boolean, OD(595nm) > 0.1 after incubation in media containing 20µg/ml nalidixic acid and 12.5µg/ml streptomycin.
  - **p256**: integer 0-255, representation of boolean growth variables as an 8-bit integer.
  - **p**: phenotype name. A corresponds to nalidix acid, B to streptomycin.

## 20180424-sensitiveCommunity.csv

  Results from experiment with a completely sensitive community.
  created from 'raw' data files using the script `prepareData.R` script.

### Community
  |p  | p256|         f|
  |:--|----:|---------:|
  |AB |    3|         0|
  |A  |    5|         0|
  |B  |   15|         0|
  |WT |    1| 0.8510638|
  |U  |    0| 0.1489362|

### Treatments
| plate|conditions_name     |drugs       |concentration | period|probability |description                                            |
|-----:|:-------------------|:-----------|:-------------|------:|:-----------|:------------------------------------------------------|
|     1|no treatment        |none        |1             |      1|1           |no treatment                                           |
|     2|monotherapyA        |A           |1             |      1|1           |1-40: Nal 20ug/ml (2xMIC) /40+: Sm 12.5ug/ml (2xMIC)   |
|     3|monotherapyB        |B           |1             |      1|1           |1-40: Sm 12.5ug/ml (2xMIC) /40+: Nal 20ug/ml (2xMIC)   |
|     4|combination therapy |AB          |1             |      1|1           |Combination therapy: Nal20 + Sm12.5                      |
|     5|cycling             |A, B        |1,1           |      2|1           |Cycling: Nal / Sm, period = 2                            |
|     6|mixing              |A, B        |1,1           |      1|0.5, 0.5    |Mixing: Nal / Sm, probability = 0.5                    |

### Rates

| turnover| infection|
|--------:|---------:|
|      0.2|       0.3|

### Columns (`col_types = "iiciiicicdic"`)
 - **plate**: integer, number of plate
 - **well**: integer, number of well on plate
 - **row**: character [A-P], row on plate
 - **col**: integer 1-24, column on plate
 - **rep**: integer 1-4, number of replicate
 - **transfer**: integer, transfer number
 - **turnoverStrain**: character, phenotype of strain replacing this culture (see above for phenotype definition)
 - **infectionToWell**: integer 1-384, well that is infected by this culture
 - **treatment**: character, drug used to treat this well (see above)
 - **OD**: double, OD(595) after incubation
 - **p256**: integer 0-255, phenotype after incubation (see above)
 - **p**: character, phenotype after incubation (see above)

## 20190627-doubleResistantCommunity.csv

  same as above except community is modelled to contain both single resistant and double resistant bacteria.

### Community
|p  | p256|         f|
|:--|----:|---------:|
|AB |    3| 0.0531915|
|A  |    5| 0.1063830|
|B  |   15| 0.1063830|
|WT |    1| 0.5212766|
|U  |    0| 0.2127660|


## 20190816-singleResistantCommunity.csv

  same as above except community is modelled to contain only single resistant bacteria.

### Community
|p  | p256|         f|
|:--|----:|---------:|
|AB |    3|         0|
|A  |    5| 0.1063830|
|B  |   15| 0.1063830|
|WT |    1| 0.5744681|
|U  |    0| 0.2127660|


## 20191124-combinationTreatments.csv

  same experiment parameters as 20180424-SensitiveCommunity.csv, but only combination treatments using different concentrations (different multiples of MIC)

| plate|conditions_name |drugs  | concentration| period| probability|description                         |
|-----:|:---------------|:------|-------------:|------:|-----------:|:-----------------------------------|
|     1|no treatment    |none   |             1|      1|           1|no treatment                        |
|     2|combo: 2/2      |AB22   |             1|      1|           1|Combination therapy: Nal20 + Sm12.5  |
|     3|combo: 1/2      |AB12   |             1|      1|           1|Combination therapy: Nal10 + Sm12.5  |
|     4|combo: 2/1      |AB21   |             1|      1|           1|Combination therapy: Nal20 + Sm6.25  |
|     5|combo: 1/1      |AB11   |             1|      1|           1|Combination therapy: Nal10 + Sm6.25  |
|     6|combo: 0.5/1    |AB051  |             1|      1|           1|Combination therapy: Nal5 + Sm6.25   |
|     7|combo: 1/0.5    |AB105  |             1|      1|           1|Combination therapy: Nal10 + Sm3.125 |
|     8|combo: 0.5/0.5  |AB0505 |             1|      1|           1|Combination therapy: Nal5 + Sm3.125  |

## 'F' files (e.g. 20180424-F-sensitiveCommunity.csv)

Tables of phenotype frequencies in the experiments

### Columns (col_types = "iiicid")
- **plate**: integer, number of plate
- **transfer**: integer, transfer number
- **rep**: integer 1-4, number of replicate
- **p**: character, phenotype after incubation (see above)
- **n**: integer, number of cultures with phenotype
- **f**: double, frequency of phenotype = n/94

## 'meanF' files (e.g. 20180424-meanF-sensitiveCommunity.csv)

Tables of mean phenotype frequencies in the experiments

### Columns (col_types = "iicdidd")
- **plate**: integer, number of plate
- **transfer**: integer, transfer number
- **p**: character, phenotype after incubation (see above)
- **mean**: double, mean frequency of phenotype
- **n**: integer, number of replicates
- **se**: double,standard error of the mean = sd/sqrt(n)
- **max**: double, maximal frequency of phenotyp across n replicates
- **min**: double, minimal frequency of phenotyp across n replicates

## raw data files (in Folder /raw/)

### Columns (col_types = "iiciiicicdllllc")
- **plate**: integer, number of plate
- **well**: integer, number of well on plate
- **row**: character [A-P], row on plate
- **col**: integer 1-24, column on plate
- **rep**: integer 1-4, number of replicate
- **transfer**: integer, transfer number
- **turnoverStrain**: character, phenotype of strain replacing this culture (see above for phenotype definition)
- **infectionToWell**: integer 1-384, well that is infected by this culture
- **treatment**: character, drug used to treat this well (see above)
- **OD**: double, OD(595) after incubation
- **agarN**: boolean, visible growth on MS Agar containing no drugs.
- **agarA**: boolean, visible growth on MS Agar containing 40µg/ml nalidixic acid.
- **agarB**: boolean, visible growth on MS Agar containing 100µg/ml streptomycin.
- **agarAB**: boolean 40µg/ml nalidixic acid (agarA) and 100µg/ml streptomycin (agarAB)
- **type**: character, either 'exp' or 'blank', indicating blank wells.

## simulation results (in Folder /simulation/)

### Columns (col_types = "iicd")
- **plate**: integer, number of plate
- **transfer**: integer, transfer number
- **p**: character, phenotype (see above)
- **f**: double, frequency of phenotype
