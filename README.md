# Dimensionality repository
This repository contains all the [functions](functions), [data](data) and scripts needed to reproduce the analysis presented in [Nakamura's et al paper](https://www.biorxiv.org/content/10.1101/508002v3). To run the entire code the user can download the zip file.

#[Data](data)

This folder contains data of South America small mammals (cricetids and marsupials) used to illustrate the application of dimensionality framework. The following archives contains:

#[Functions](functions)

This folder contains all R functions needed to run the analysis of dimensionality:

- dimensionality.R: R function to calculate the Evenness of Eigenvalues metric proposed by [Stevens and Tello (2014)](https://www.researchgate.net/publication/262605747_On_the_measurement_of_dimensionality_of_biodiversity)

- ImportanceVal_V2.R: function used to calculate Importance Values following the method presented in [Nakamura et al (2017)](https://onlinelibrary.wiley.com/doi/full/10.1111/aec.12529)

- Camargo_function.R: function used to calculate Camargo evenness index.

- simm_comm.R: function used to simulate metacommunities accordingly to one of the scenarios presented in [Nakamura et al (2017)](https://www.biorxiv.org/content/10.1101/508002v3)

- test_TryCatch_metricsFunc_27-7-2019.R: R function used to overcome an aleatory error presented in the simulation of null metacommunities and diversity metrics.
