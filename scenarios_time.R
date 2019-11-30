library(magrittr)
library(dplyr)
library(stringr)
library('gRain')
source('utils.R')

yn  <- c('yes', 'no')
### First scenario, var_Z*/var_Z**=0.675

# Build DAG
A0  <- cptable( ~ A0, values=c(0.5, 0.5), levels=yn)
H <- cptable( ~ H, values=c(0.5, 0.5), levels=yn)
R <- cptable( ~ R | A0 + H, values=c(0.3, 0.7, 0.7, 0.3, 0.4, 0.6, 0.6, 0.4), levels=yn)
A1 <- cptable( ~ A1 | H, values=c(0.8, 0.2, 0.2, 0.8), levels=yn)
Q <- cptable( ~ Q | R, values=c(0.1, 0.9, 0.9, 0.1), levels=yn)
Y <- cptable( ~ Y | A1 + Q, values=c(0.2, 0.8, 0.3, 0.7, 0.7, 0.3, 0.8, 0.2), levels=yn)
cptlist <- compileCPT(list(A0, R, H, A1, Q, Y))
dag <- grain(cptlist)

# Compute variances
var_1 <- get_var_nonparam_influence(L0=c(), L1=c('Q'), dag=dag)
var_8 <- get_var_nonparam_influence(L0=c('H'), L1=c('Q'), dag=dag)
var_1/var_8

### Second scenario, var_Z*/var_Z**=1.08

# Build DAG
A0  <- cptable( ~ A0, values=c(0.5, 0.5), levels=yn)
H <- cptable( ~ H, values=c(0.5, 0.5), levels=yn)
R <- cptable( ~ R | A0 + H, values=c(0.01, 0.99, 0.01, 0.99, 0.99, 0.01, 0.99, 0.01), levels=yn)
A1 <- cptable( ~ A1 | H, values=c(0.8, 0.2, 0.8, 0.2), levels=yn)
Q <- cptable( ~ Q | R, values=c(0.99, 0.01, 0.01, 0.99), levels=yn)
Y <- cptable( ~ Y | A1 + Q, values=c(0.2, 0.8, 0.2, 0.8, 0.8, 0.2, 0.8, 0.2), levels=yn)
cptlist <- compileCPT(list(A0, R, H, A1, Q, Y))
dag <- grain(cptlist)

# Compute variances
var_1 <- get_var_nonparam_influence(L0=c(), L1=c('Q'), dag=dag)
var_8 <- get_var_nonparam_influence(L0=c('H'), L1=c('Q'), dag=dag)
var_1/var_8