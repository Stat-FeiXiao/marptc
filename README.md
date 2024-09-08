# qifptc
An R package for fitting marginal semiparametric promotion time cure model under clustered survival data.
- We consider two marginal methods: *generalized estimating equations*  (**GEE**)  and *quadratic inference functions*   (**QIF**).

## Package description and included main functions

Installation of this package can be done locally after downloading the package manually from this github website. We will also upload this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package. Currently, it can be loaded using R command
```R
devtools::install_github("Stat-FeiXiao/qifptc")
library(qifptc)
```

The main function included in our R package is *qifptc()* and there is also a function *print.qifptc()* for printing fitted results with a better presentation. To sum up, they can be called via:
- **smgeecure**: fit the models in various ways with synopsis
```R
qifptc(formula, data, id, Var = TRUE,Ibeta=NULL, stad=TRUE,boots=FALSE,nboot=100,
       method = "GEE", corstr="independence",itermax = 100, eps = 1e-06) 
```
- **print.qifptc**: print outputted results from the previous function *qifptc()* with syntax
```R
print.qifptc(fit)
```

## Two numerical illustrations

### An example using a periodontal disease data is shown below:

```R
## library
library(survival)
library(CureAuxSP)

#### Data preparation
```R
data(Teeth)
Teeth$Mg <- Teeth$x10 # 1 for Mucogingival defect
Teeth$Endo <- Teeth$x16 # 1 for endo Therapy
Teeth$Decay <- Teeth$x21 # 1 for decayed tooth
Teeth$Gender <- Teeth$x49 # 1 for female
```

## plot a figure to show the existence of a cure fraction
```R
ggsurvplot(survival::survfit(survival::Surv(time, event) ~ x49, data = Data), 
           ylim = c(0.6,1),
           ylab = "Survival Probability", xlab = "Survival Time (in Years)", 
           censor.shape="+",
           legend.title = "Gender",
           legend.labs = c("Female","Male")
)
```


#### Fit the marginal semi-parametric promotion time cure model using GEE method
- exchangeable correlation
```R
teeth.gee.ex <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "exchangeable", Var = T, nboot = 100
)
print.ptcqif(teeth.gee.ex)
```
- AR(1) correlation
```R
teeth.gee.ar1 <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "exchangeable", Var = T, nboot = 100
)
print.smgeecure(teeth.gee.ar1)
```
- independence correlation
```R
teeth.gee.ind <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "independence", Var = T, nboot = 100
)
print.ptcqif(teeth.gee.ind)
```
#### Fit the marginal semi-parametric promotion time cure model using QIF method
- exchangeable correlation
```R
teeth.qif.ex <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "exchangeable", Var = T, nboot = 100
)
print.ptcqif(teeth.qif.ex)
```
- AR(1) correlation
```R
teeth.qif.ar1 <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "exchangeable", Var = T, nboot = 100
)
print.ptcqif(teeth.qif.ar1)
```
- independence correlation
```R
teeth.qif.ind <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "independence", Var = T, nboot = 100
)
print.ptcqif(teeth.qif.ind)
```
