# Treatment Strategies Data Files

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

### 20180424-sensitiveCommunity.csv

  Results from experiment with a completely sensitive community.
  created from 'raw' data files using the script `prepareData.R` script.

#### Community
  |p  | p256|         f|
  |:--|----:|---------:|
  |AB |    3|         0|
  |A  |    5|         0|
  |B  |   15|         0|
  |WT |    1| 0.8510638|
  |U  |    0| 0.1489362|

#### Treatments
| plate|conditions_name     |drugs       |concentration | period|probability |description                                            |
|-----:|:-------------------|:-----------|:-------------|------:|:-----------|:------------------------------------------------------|
|     1|no treatment        |none        |1             |      1|1           |no treatment                                           |
|     2|monotherapyA        |A           |1             |      1|1           |1-40: Nal 20ug/ml (2xMIC) /40+: Sm 12.5ug/ml (2xMIC)   |
|     3|monotherapyB        |B           |1             |      1|1           |1-40: Sm 12.5ug/ml (2xMIC) /40+: Nal 20ug/ml (2xMIC)   |
|     4|combination therapy |AB          |1             |      1|1           |Combination therapy: Nal20 + Sm12.5                      |
|     5|cycling             |A, B        |1,1           |      2|1           |Cycling: Nal / Sm, period = 2                            |
|     6|mixing              |A, B        |1,1           |      1|0.5, 0.5    |Mixing: Nal / Sm, probability = 0.5                    |

#### Rates

| turnover| infection|
|--------:|---------:|
|      0.2|       0.3|

#### Columns (`col_types = "iiciiicicdic"`)
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

### 20190627-doubleResistantCommunity.csv

  same as above except community is modelled to contain both single resistant and double resistant bacteria.

#### Community
|p  | p256|         f|
|:--|----:|---------:|
|AB |    3| 0.0531915|
|A  |    5| 0.1063830|
|B  |   15| 0.1063830|
|WT |    1| 0.5212766|
|U  |    0| 0.2127660|


### 20190816-singleResistantCommunity.csv

  same as above except community is modelled to contain only single resistant bacteria.

#### Community
|p  | p256|         f|
|:--|----:|---------:|
|AB |    3|         0|
|A  |    5| 0.1063830|
|B  |   15| 0.1063830|
|WT |    1| 0.5744681|
|U  |    0| 0.2127660|


### 20191124-combinationTreatments.csv

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

### 'F' files (e.g. 20180424-F-sensitiveCommunity.csv)

Tables of phenotype frequencies in the experiments

#### Columns (col_types = "iiicid")
- **plate**: integer, number of plate
- **transfer**: integer, transfer number
- **rep**: integer 1-4, number of replicate
- **p**: character, phenotype after incubation (see above)
- **n**: integer, number of cultures with phenotype
- **f**: double, frequency of phenotype = n/94

### 'meanF' files (e.g. 20180424-meanF-sensitiveCommunity.csv)

Tables of mean phenotype frequencies in the experiments

#### Columns (col_types = "iicdidd")
- **plate**: integer, number of plate
- **transfer**: integer, transfer number
- **p**: character, phenotype after incubation (see above)
- **mean**: double, mean frequency of phenotype
- **n**: integer, number of replicates
- **se**: double,standard error of the mean = sd/sqrt(n)
- **max**: double, maximal frequency of phenotyp across n replicates
- **min**: double, minimal frequency of phenotyp across n replicates

### raw data files (in Folder /raw/)

#### Columns (col_types = "iiciiicicdllllc")
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
