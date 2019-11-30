library(magrittr)
library(dplyr)
library(stringr)
library('gRain')
source('utils.R')

yn  <- c('yes', 'no')

# Build DAG, Z1~A strong, Z1~Z2, U~Z2, Y~U weak, var_Z*/var_Z**=0.0396
Z1  <- cptable( ~ Z1, values=c(0.5, 0.5), levels=yn)
U  <- cptable( ~ U, values=c(0.5, 0.5), levels=yn)
A <- cptable( ~ A | Z1, values=c(0.01, 0.99, 0.99, 0.01), levels=yn)
Z2 <- cptable( ~ Z2 | Z1 + U, values=c(0.6, 0.4, 0.5, 0.5, 0.6, 0.4, 0.4, 0.6), levels=yn)
Y <- cptable( ~ Y | U + A, values=c(0.6, 0.4, 0.5, 0.5, 0.6, 0.4, 0.4, 0.6), levels=yn)
cptlist <- compileCPT(list(Z1, U, A, Z2, Y))
dag <- grain(cptlist)

# Compute variances
var_1 <- get_var_nonparam_influence_point(Z=c(), dag=dag)
var_2 <- get_var_nonparam_influence_point(Z=c('Z1', 'Z2'), dag=dag)
var_3 <- get_var_nonparam_influence_point(Z=c('Z1'), dag=dag)
var_1/var_2
var_1/var_3


# Build DAG, Z1~A, Z1~Z2 weak, U~Z2 strong, var_Z*/var_Z**=1.44
Z1  <- cptable( ~ Z1, values=c(0.5, 0.5), levels=yn)
U  <- cptable( ~ U, values=c(0.5, 0.5), levels=yn)
A <- cptable( ~ A | Z1, values=c(0.5, 0.5, 0.5, 0.5), levels=yn)
Z2 <- cptable( ~ Z2 | Z1 + U, values=c(0.1, 0.9, 0.1, 0.9, 0.9, 0.1, 0.9, 0.1), levels=yn)
Y <- cptable( ~ Y | U + A, values=c(0.99, 0.01, 0.01, 0.99, 0.99, 0.01, 0.01, 0.99), levels=yn)
cptlist <- compileCPT(list(Z1, U, A, Z2, Y))
dag <- grain(cptlist)

# Compute variances
var_1 <- get_var_nonparam_influence_point(Z=c(), dag=dag)
var_2 <- get_var_nonparam_influence_point(Z=c('Z1', 'Z2'), dag=dag)
var_3 <- get_var_nonparam_influence_point(Z=c('Z1'), dag=dag)
var_1/var_2
var_1/var_3

