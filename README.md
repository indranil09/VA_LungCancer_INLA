# Lung Cancer Incidence in Virginia: A spatial zipcode-level analysis via INLA
### Authors: Indranil Sahoo, Jinlei Zhao, Xiaoyan Deng, Myles Cockburn, Katherine Tossas, Robert Winn, Dipankar Bandyopadhyay

In this paper, we analyze overdispersed lung cancer counts in Virginia during the five-year period of 2014 - 2018, obtained at the zip code level using Poisson and negative binomial spatial regression models. The spatial correlation across zip codes is modeled using the popular Conditional Autoregressive (CAR) model. Since covariates historically related to lung cancer have missing values in our data, we model the missing values as latent Gaussian Markov Random Fields (GMRFs), and use a Bayesian hierarchical modeling framework to fit the model within the `INLA` framework. The latent effects required to impute the missing covariates are implemented using the `MIINLA` package in `R` (available from Github repository at https://github.com/becarioprecario/MIINLA). Finally, a comparative study among models has been performed to identify the best model in the presence of overdispersion, spatial correlation and missing covariate information. 

The dataset used in this study has been compiled in `VAlungcancerzip.sas7bdat`. To load the dataset into your current workspace, use the following code: 
```{r}
require(sas7bdat)
valung = read.sas7bdat(file = "VAlungcancerzip.sas7bdat", encoding="", debug=FALSE)
```
The code for this study has been organized as follows: 
* Data loading and pre-processing
* Exploratory analysis of lung cancer counts
* Implementing GLM in INLA with simultaneous missing covariates imputation
* Plots from the manuscript
* Comparison with no-imputation models
