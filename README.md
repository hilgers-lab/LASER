# SaiLoR
## Splicing influence estimation by Long Reads
<!-- badges: start -->

![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/hilgers-lab/SaiLoR)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-brightgreen)](https://github.com/hilgers-lab/SaiLoR/graphs/contributors)
[![Install](https://img.shields.io/badge/Install-Github-brightgreen)](#installation)
[![Downloads](https://img.shields.io/github/downloads/hilgers-lab/SaiLoR/total)]()
![GitHub](https://img.shields.io/github/license/hilgers-lab/SaiLoR)
<!-- badges: end -->

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
vignette("SaiLoR")
```

# Release 

Initial Release 0.1.0

Release date: 20th Dec 2022
This release corresponds to the SaiLoR version used by Alfonso-Gonzalez et al. manuscript

## Contact

Developer Carlos Alfonso-Gonzalez. For questions or feedback you can contact:

alfonso@ie-freiburg.mpg.de
