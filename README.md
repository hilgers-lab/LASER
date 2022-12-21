# SaiLoR
## Splicing influence estimation by Long Reads

SaiLoR determines the regulatory connections between exons, 5' ends, and 3' ends by analyzing every read as a complete transcript and using multinomial testing to evaluate the frequency of co-occurrence among these features. 


### Installation 

```
install.packages("devtools")
devtools::install_github("hilgers-lab/SaiLoR", build = TRUE, build_vignettes = TRUE)
```
### Usage 
The vignette contains a step by step guide for data processing and identification of exon-couplings.
```
library(SaiLoR)
browseVignettes("SaiLoR")
```

## Contact

Developer Carlos Alfonso-Gonzalez. For questions or feedback you can contact:

alfonso@ie-freiburg.mpg.de


