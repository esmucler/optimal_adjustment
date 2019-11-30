fix_names <- function(df, L0, L1){
  cnames <- colnames(df)
  for (ind in seq_along(cnames)){
    if (cnames[ind] %in% L0){
      cnames[ind] <- paste(cnames[ind], '_L0', sep='')
    } else if (cnames[ind] %in% L1){
      cnames[ind] <- paste(cnames[ind], '_L1', sep='')
    }
  }
  colnames(df) <- cnames
  return(df)
}

yn_to_bin <- function(x){
  if (x=='yes'){
    return(1)
  } else if (x=='no'){
    return(0)
  }
}

###### Time dependent

get_adj_var_names <- function(df){
  var_names_0 <- colnames(df)[which(str_detect(colnames(df), 'L0'))]
  var_names_1 <- colnames(df)[which(str_detect(colnames(df), 'L1'))]
  return(list('var_names_0'=var_names_0, 'var_names_1'=var_names_1))
}

get_sel <- function(df, var_names, values){
  if (length(var_names)==1){
    sel <- apply(t(df[var_names])==values, 2, all)
  } else if(length(var_names)>1) {
    sel <- apply(df[var_names], 1, function(x) {all(x==values)})
  } else {
    sel <- rep(TRUE, nrow(df))
  }
  return(sel)
}

get_pi0 <- function(L0_value='yes', L0, dag){
  if (length(L0)>=1){
    pi0_df <- querygrain(dag, nodes = c('A0', L0), type='conditional', result='data.frame')
    pi0_df <- fix_names(pi0_df, L0=L0, L1=NULL)
    var_names <- get_adj_var_names(df=pi0_df)[['var_names_0']]
    sel <- get_sel(df=pi0_df, var_names=var_names, values=L0_value)
    out <- pi0_df[sel & pi0_df['A0']=='yes', 'Freq']
    out_comp <- pi0_df[sel & pi0_df['A0']=='no', 'Freq']
    stopifnot( abs(out + out_comp - 1)<= 1e-5 )
    return(out)
  } else {
    pi0_df <- querygrain(dag, nodes = 'A0', type='marginal', result='data.frame')
    out <- pi0_df$A0[pi0_df$A0['A0']=='yes', 'Freq']
    out_comp <- pi0_df$A0[pi0_df$A0['A0']=='no', 'Freq']
    stopifnot( abs(out + out_comp - 1)<= 1e-5 )
    return(out)
  }
}

get_pi1 <- function(L0_value='yes', L1_value='yes', L0, L1, dag){
  pi1_df <- querygrain(dag, nodes = c('A1', c(L0, L1, 'A0')), type='conditional', result='data.frame')
  pi1_df <- fix_names(pi1_df, L0=L0, L1=L1)
  var_names_0 <- get_adj_var_names(df=pi1_df)[['var_names_0']]
  var_names_1 <- get_adj_var_names(df=pi1_df)[['var_names_1']]
  sel_0 <- get_sel(df=pi1_df, var_names=var_names_0, values=L0_value)
  sel_1 <- get_sel(df=pi1_df, var_names=var_names_1, values=L1_value)
  out <- pi1_df[sel_0 & sel_1 & pi1_df['A0']=='yes' & pi1_df['A1']=='yes', 'Freq']
  out_comp <- pi1_df[sel_0 & sel_1 & pi1_df['A0']=='yes' & pi1_df['A1']=='no', 'Freq']
  stopifnot( abs(out + out_comp - 1)<= 1e-5 )
  return(out)
}

get_pl1 <- function(L0_value, L1_value, L0, L1, dag){
  # prob of L1 given L0 and A0=1
  if(length(L1)>=1){
    pl1_df <- querygrain(dag, nodes = c(L1, 'A0', L0),type='conditional', result='data.frame')
    pl1_df <- fix_names(pl1_df, L0=L0, L1=L1)
    var_names_0 <- get_adj_var_names(df=pl1_df)[['var_names_0']]
    var_names_1 <- get_adj_var_names(df=pl1_df)[['var_names_1']]
    sel_0 <- get_sel(df=pl1_df, var_names=var_names_0, values=L0_value)
    sel_1 <- get_sel(df=pl1_df, var_names=var_names_1, values=L1_value)
    return(pl1_df[sel_0 & sel_1 & pl1_df['A0']=='yes','Freq'])
  } else {
    return(0)
  }
}

get_b1 <- function(L0_value, L1_value, L0, L1, dag){
  b1_df <- querygrain(dag, nodes = c('Y', c(L0, L1, 'A0', 'A1')),type='conditional', result='data.frame')
  b1_df <- fix_names(b1_df, L0=L0, L1=L1)
  var_names_0 <- get_adj_var_names(df=b1_df)[['var_names_0']]
  var_names_1 <- get_adj_var_names(df=b1_df)[['var_names_1']]
  sel_0 <- get_sel(df=b1_df, var_names=var_names_0, values=L0_value)
  sel_1 <- get_sel(df=b1_df, var_names=var_names_1, values=L1_value)
  out <- b1_df[sel_0 & sel_1 & b1_df['A0']=='yes' & b1_df['A1']=='yes' & b1_df['Y']=='yes', 'Freq']
  out_comp <- b1_df[sel_0 & sel_1 & b1_df['A0']=='yes' & b1_df['A1']=='yes' & b1_df['Y']=='no', 'Freq']
  stopifnot( abs(out + out_comp - 1)<= 1e-5 )
  return(out)
}

get_b0 <- function(L0_value, L0, L1, dag){
  if (length(L1) >= 1){
    all_combs_L1 <- expand.grid(lapply(1:length(L1), function(x)  {c('yes','no')}))
    colnames(all_combs_L1) <- L1
    out <- 0
    for (ind in 1:nrow(all_combs_L1)){
      prob <- get_pl1(L0_value, L1_value=all_combs_L1[ind,], L0=L0, L1=L1, dag=dag)
      val <- get_b1(L0_value=L0_value, L1_value=all_combs_L1[ind,], L0=L0, L1=L1, dag=dag)
      out <- out +  prob * val
    }
    return(out)
  } else {
    return(get_b1(L0_value=L0_value, L1_value='no', L0=L0, L1=L1, dag=dag))
  }
}


target_scalar <- function(x, dag){
  # x = (a0, a1, h, y)
  L0 <- c()
  L1 <- 'H'
  A0_bin <- yn_to_bin(x$A0)
  A1_bin <- yn_to_bin(x$A1)
  Y_bin <- yn_to_bin(x$Y)
  pi0 <- get_pi0(L0_value=c(), L0=L0, dag=dag)
  pi1 <- get_pi1(L0_value=c(), L1_value=x$H, L0=L0, L1=L1, dag=dag)
  val <- (A0_bin * A1_bin * Y_bin) / (pi0 * pi1)
  prob <- x$Freq
  out <- val * prob
  return(out)
}

get_target_chi <- function(dag){
  all_combs_target_probs <- querygrain(dag, nodes = c('A0', 'A1', 'H', 'Y'), type='joint', result='data.frame')
  out <- 0
  for (ind in 1:nrow(all_combs_target_probs)){
    out <- out + target_scalar(all_combs_target_probs[ind,],dag=dag)
  }
  return(out)
}

influence_nonparam_scalar <- function(x, L0, L1, dag, chi, square){
  # x = (a0, a1, l0, l1, y)
  var_names_0 <- get_adj_var_names(df=x)[['var_names_0']]
  var_names_1 <- get_adj_var_names(df=x)[['var_names_1']]
  
  A0_bin <- yn_to_bin(as.character(x$A0))
  A1_bin <- yn_to_bin(as.character(x$A1))
  L0_value <- c(t(x[var_names_0]))
  L1_value <- c(t(x[var_names_1]))
  Y_bin <- yn_to_bin(as.character(x$Y))
  pi0 <- get_pi0(L0_value=L0_value, L0=L0, dag=dag)
  pi1 <- get_pi1(L0_value=L0_value, L1_value=L1_value, L0=L0, L1=L1, dag=dag)
  b1 <- get_b1(L0_value=L0_value, L1_value=L1_value, L0=L0, L1=L1, dag=dag)
  b0 <- get_b0(L0_value=L0_value, L0=L0, L1=L1, dag=dag)
  term1 <- A0_bin * A1_bin * (Y_bin - chi) / (pi0 * pi1)
  term2 <- (A0_bin/pi0) * (A1_bin/pi1 - 1) * (b1 - b0)
  term3 <- (A0_bin/pi0 - 1) * (b0 - chi)
  if (square){
    val <- (term1 - term2 - term3)^{2}
  } else{
    val <- (term1 - term2 - term3)
  }
  prob <- x$Freq
  out <- val * prob
  return(out)
}

get_var_nonparam_influence <- function(L0, L1, dag){
  chi <- get_target_chi(dag=dag)
  all_combs_probs <- querygrain(dag, nodes = c('A0', 'A1', 'Y', L0, L1), type='joint', result='data.frame')
  all_combs_probs <- fix_names(all_combs_probs, L0, L1)
  out <- 0
  for (ind in 1:nrow(all_combs_probs)){
    out <- out + influence_nonparam_scalar(x=all_combs_probs[ind,], L0=L0, L1=L1, chi=chi, dag=dag, square=TRUE)
  }
  return(out)
}

get_mean_nonparam_influence <- function(L0, L1, dag){
  chi <- get_target_chi(dag=dag)
  all_combs_probs <- querygrain(dag, nodes = c('A0', 'A1', 'Y', L0, L1), type='joint', result='data.frame')
  all_combs_probs <- fix_names(all_combs_probs, L0, L1)
  out <- 0
  for (ind in 1:nrow(all_combs_probs)){
    out <- out + influence_nonparam_scalar(x=all_combs_probs[ind,], L0=L0, L1=L1, chi=chi, dag=dag, square=FALSE)
  }
  return(out)
}

###### Time independent

get_pi_point <- function(Z_value, Z, dag){
  # Compute P(A=1 | Z=Z_value)
  if (length(Z)>=1){
    pi_df <- querygrain(dag, nodes = c('A', Z), type='conditional', result='data.frame')
    sel <- get_sel(df=pi_df, var_names=Z, values=Z_value)
    out <- pi_df[sel & pi_df['A']=='yes', 'Freq']
    out_comp <- pi_df[sel & pi_df['A']=='no', 'Freq']
    stopifnot( abs(out + out_comp - 1)<= 1e-5 )
  } else {
    pi_df <- querygrain(dag, nodes = c('A'), type='joint', result='data.frame')
    out <- pi_df[pi_df['A']=='yes', 'Freq']
  }
  
  return(out)
}

get_b_point <- function(Z_value, Z, dag){
  # Compute P(A=1 | Z=Z_value)
  pi_df <- querygrain(dag, nodes = c('Y', 'A', Z), type='conditional', result='data.frame')
  sel <- get_sel(df=pi_df, var_names=Z, values=Z_value)
  out <- pi_df[sel & pi_df['A']=='yes' & pi_df['Y']=='yes', 'Freq']
  return(out)
}

target_scalar_point <- function(x, Z, dag){
  freq <- x$Freq
  x_vals <- x[1:(length(x) - 1)]
  len <- length(x_vals)
  y <- as.character(x_vals$Y)
  a <- as.character(x_vals$A)
  z <- apply(select_at(x, vars(-c('A', 'Y', 'Freq'))), 2, as.character)
  y_bin <- yn_to_bin(y)
  a_bin <- yn_to_bin(a)
  pi_point <- get_pi_point(Z_value = z, Z=Z, dag=dag)
  
  val <- a_bin * y_bin / pi_point
  out <- val * freq
  return(out)
}

get_target_chi_point <- function(Z, dag){
  all_combs_target_probs <- querygrain(dag, nodes = c(Z, 'A', 'Y'), type='joint', result='data.frame')
  out <- 0
  for (ind in 1:nrow(all_combs_target_probs)){
    out <- out + target_scalar_point(all_combs_target_probs[ind,], Z=Z, dag=dag)
  }
  return(out)
}

influence_nonparam_point_scalar <- function(x, Z, dag, chi, square){
  freq <- x$Freq
  x_vals <- x[1:(length(x) - 1)]
  len <- length(x_vals)
  y <- as.character(x_vals$Y)
  a <- as.character(x_vals$A)
  z <- apply(select_at(x, vars(-c('A', 'Y', 'Freq'))), 2, as.character)
  y_bin <- yn_to_bin(y)
  a_bin <- yn_to_bin(a)
  pi_point <- get_pi_point(Z_value = z, Z=Z, dag=dag)
  b_point <- get_b_point(Z_value = z, Z=Z, dag=dag)
  
  val <- a_bin/pi_point * (y_bin - b_point) + b_point - chi
  
  if (square){
    val <- val^2
  } 
  
  out <- val * freq
  return(out)
}

get_var_nonparam_influence_point <- function(Z, dag){
  chi <- get_target_chi_point(Z=Z, dag=dag)
  all_combs_probs <- querygrain(dag, nodes = c(Z, 'A', 'Y'), type='joint', result='data.frame')
  out <- 0
  for (ind in 1:nrow(all_combs_probs)){
    out <- out + influence_nonparam_point_scalar(x=all_combs_probs[ind,], Z=Z, chi=chi, dag=dag, square=TRUE)
  }
  return(out)
}

get_mean_nonparam_influence_point <- function(Z, dag){
  chi <- get_target_chi_point(Z=Z, dag=dag)
  all_combs_probs <- querygrain(dag, nodes = c(Z, 'A', 'Y'), type='joint', result='data.frame')
  out <- 0
  for (ind in 1:nrow(all_combs_probs)){
    out <- out + influence_nonparam_point_scalar(x=all_combs_probs[ind,], Z=Z, chi=chi, dag=dag, square=FALSE)
  }
  return(out)
}
