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
smgeecure(formula, cureform, data, id, model = c("aft", "ph"),
          corstr = c("independence", "exchangeable", "ar1"),
          Var = TRUE, nboot = 100, stdz = TRUE, esmax = 20, eps = 1e-04)
```
- **print.qifptc**: print outputted results from the previous function *qifptc()* with syntax
```R
print.qifptc(fit)
```

## Two numerical illustrations

### An example using a real dataset from TCGA program is shown below:

```R
## library
library(survival)
library(CureAuxSP)

#### Data preparation
```R
data(tonsil)
tonsil <- tonsil[-c(141,136,159),]
tonsil$Sex <- ifelse(tonsil$Sex == 1, 0, 1) # 1="Female"
tonsil$Cond <- ifelse(tonsil$Cond == 1, 0, 1) # 0=no disability
tonsil$T <- ifelse(tonsil$T < 4, 0, 1)
tonsil$Grade2 <- ifelse(tonsil$Grade==2,1,0)
tonsil$Grade3 <- ifelse(tonsil$Grade==3,1,0)
table(tonsil$Inst)
```

## plot a figure to show the existence of a cure fraction
```R
plot(
  survival::survfit(survival::Surv(yobs, delta) ~ 1, data = sdata.TCGA), 
  conf.int = T, mark.time = TRUE, lwd = 2,
  ylab = "Survival Probability", xlab = "Survival Time (in Years)", 
  xlim = c(0,25), ylim = c(0,1)
)
```
![Figure_Tonsil_KM_SexTumorsize](https://github.com/user-attachments/assets/7874b5c8-46ea-4235-af6b-6e6e49c592ac)


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
