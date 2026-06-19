# XGamma-RightCensures
R functions for estimation of the parameter of the XGamma distribution 

**"Maximum likelihood, Bayesian inference, and prediction for the X-Gamma distribution with right-censored survival data"**

This repository provides functions for estimation, Bayesian inference, interval estimation, prediction, and Monte Carlo simulation under the X-Gamma distribution in the presence of right-censored survival data.

## Features
# Functions used to get:
* Maximum likelihood estimation (MLE)
* Importance Sampling (IS)
* Independent Metropolis-Hastings (IMH)
* Random Walk Metropolis (RWM)
* Credible intervals and HPD intervals
* Predictive inference

## Repository Structure

```text
R/
├── MLE.XGamma.R
├── Bayes.XGamma.R
├── Predict.XGamma.R
├── Intervals.XGamma.R
├── MonteCarlo.XGamma.R

data/
├── lung.csv

examples/
├── Example1.R

paper/
├── manuscript.pdf
```

## Requirements

Required R packages:

```r
install.packages(c(
  "survival",
  "coda",
  "progress",
  "parallel"
))
```


## Reproducing the Results

## Reproducing the Results

The results reported in the paper can be reproduced by following the steps below:

1. Compile the file `fit.functionsXGamma.R`, which contains all functions required for maximum likelihood estimation, Bayesian inference, interval estimation, and predictive analysis under the X-Gamma distribution.

```r
source("fit.functionsXGamma.R")
```

2. Run the scripts corresponding to the analyses presented in the manuscript:

* `simulation_study.R`: reproduces the Monte Carlo simulation study.
* `real_data_analysis.R`: reproduces the analyses of the real survival datasets.
* `prediction_analysis.R`: reproduces the predictive performance study.

3. All tables and figures reported in the manuscript are generated automatically by these scripts.

The computations were performed using R (version 4.x or later). Required packages can be installed using

```r
install.packages(c(
  "survival",
  "coda",
  "parallel",
  "progress"
))
```

Expected running times depend on the number of Monte Carlo replications and MCMC iterations specified in the scripts.

## Authors

Erlandson Ferreira Saraiva
Instituto de Matemática, Federal University of Mato Grosso do Sul (UFMS)

Valdemiro Piedade Vigas
Departamento de Ciências Contabeis da Faculdade de Ciências Contábeis

Natan Hilário da Silva 
Instituto de Ciência Matemática e Computação, Universidade de São Paulo Paulista (USP)

Adriano Kamimura Suzuki 
Instituto de Ciência Matemática e Computação, Universidade de São Paulo Paulista (USP)

## Citation

If you use this code, please cite:

Saraiva, E. F. (2026), Vigas, V. P., da Silva, N. H. and Suzuki, A. K.
Maximum likelihood, Bayesian inference, and prediction for the X-Gamma distribution with right-censored survival data.

## License

MIT License
