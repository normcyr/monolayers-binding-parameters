# R script to calculate binding parameters in monolayer experiments

## Description

This repository offer a complementory R script to a scientific paper which will soon be published (currently pending review - reference will be added accordingly). The goal of this script is to facilitate the calculation of binding parameters (and corresponding errors) associated with Langmuir monolayers experiments.

## How to use it

You can either clone the git repository locally or simply download the R script file named `binding-parameters-calculations.R`.

Then, using your preferred R interpreter (*eg* [RStudio](https://github.com/rstudio/rstudio); a more exhaustive list can be found [here](https://en.wikipedia.org/wiki/R_programming_language#Interfaces)), open the R script.

We provide you with an example data set (`example-data.csv`) which can also be downloaded and used as trial data. The first column has to be Π<sub>i</sub> data and the second column ΔΠ data.

Your data has to be provided as a `csv` file where the first line corresponds to headers. Alternatively, you may modify the R script to suit your data format.

The script will read your data and generate a ΔΠ vs Π<sub>i</sub> plot like the following:

![alt text][plot]

The authors do not pretend to great artistic talent. You are therefore encouraged to modify the resulting plot to your needs and artistic capabilities. R offers various packages for plot generation.

## References and further readings

References for the calculations:

- Calvez *et al*. (2009) *Biochimie* 91:718. ([Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/19345719))
- Calvez *et al*. (2011) *Langmuir* 27:1373. ([Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/21210634))
- Boisselier *et al*. (2017) *Advances in Colloid and Interface Science* 243:60. ([Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/28372794))

A web server offering similar features is available [here](http://www.crchudequebec.ulaval.ca/BindingParametersCalculator/).

[plot]: https://normcyr.github.io/img/example-plot.svg "Example plot generated by the R script"
