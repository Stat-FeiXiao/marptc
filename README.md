[KM.pdf](https://github.com/user-attachments/files/28768403/KM.pdf)# marptc
An R package for fitting marginal semiparametric promotion time cure model under clustered survival data.
We consider two marginal methods: *generalized estimating equations*  (**GEE**)  and *quadratic inference functions*   (**QIF**).

## Package description and included main functions

Installation of this package can be done locally after downloading the package manually from this github website. We will also upload this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package. Currently, it can be loaded using R command
```R
devtools::install_github("Stat-FeiXiao/marptc")
library(marptc)
```

The main function included in our R package is *marptc()* and there is also a function *print.marptc()* for printing fitted results with a better presentation. To sum up, they can be called via:
- **marptc**: fit the models in various ways with synopsis
```R
marptc(formula, data, id, Var = TRUE, Ibeta = NULL, stad = TRUE, boots = FALSE, 
       nboot = 100, method = "GEE", corstr = "independence", IC = FALSE, itermax = 100, eps = 1e-06) 
```
- **print.marptc**: print outputted results from the previous function *marptc()* with syntax
```R
print.marptc(fit)
```

## An example using a periodontal disease data is shown below:

```R
## library
library(survival)
library(survminer)

#### Data preparation
```R
data(teeth)
Data <- teeth
Data$DF <- Data$x14 # the number of decayed and filled sites, per tooth
Data$Crown <- Data$x15 # tooth has crown (0 = yes, 1 = no)
Data$Gender <- Data$x49 # Gender (0=female, 1=male)
Data$Tobacco <- Data$x51 # use of tobacco (0=yes, 1=no)
```

## plot a figure to show the existence of a cure fraction
```R
ggsurvplot(survival::survfit(survival::Surv(time, event) ~ 1, data = teeth), 
           ylim = c(0.6,1),
           ylab = "Survival Probability", xlab = "Survival Time (in Years)", 
           censor.shape="+"
)
```
<img width="794" height="443" alt="KM" src="https://github.com/user-attachments/assets/d3d3679f-dc8e-4aff-a4ba-2903936714cc" />

#### Fit the marginal semi-parametric promotion time cure model using GEE method
- exchangeable correlation
```R
teeth.gee.ex <- marptc(
        formula = Surv(time, event) ~ DF + Crown + Gender + Tobacco, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "GEE", corstr="exchangeable", data = Data
)
print.marptc(teeth.gee.ex)
```
- AR(1) correlation
```R
teeth.gee.ar1 <- marptc(
        formula = Surv(time, event) ~ DF + Crown + Gender + Tobacco, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "GEE", corstr="AR1", data = Data
)
print.marptc(teeth.gee.ar1)
```
- independence correlation
```R
teeth.gee.ind <- marptc(
        formula = Surv(time, event) ~ DF + Crown + Gender + Tobacco, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "GEE", corstr="independence", data = Data
)
print.marptc(teeth.gee.ind)
```
#### Fit the marginal semi-parametric promotion time cure model using QIF method
- exchangeable correlation
```R
teeth.qif.ex <- marptc(
        formula = Surv(time, event) ~ DF + Crown + Gender + Tobacco, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "QIF", corstr="exchangeable", data = Data
)
print.marptc(teeth.qif.ex)
```
- AR(1) correlation
```R
teeth.qif.ar1 <- marptc(
        formula = Surv(time, event) ~ DF + Crown + Gender + Tobacco, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "QIF", corstr="AR1", data = Data
)
print.marptc(teeth.qif.ar1)
```
- independence correlation
```R
teeth.qif.ind <- marptc(
        formula = Surv(time, event) ~ DF + Crown + Gender + Tobacco, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "QIF", corstr="independence", data = Data
)
print.marptc(teeth.qif.ind)
```
