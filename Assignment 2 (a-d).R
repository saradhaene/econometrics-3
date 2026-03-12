install.packages("rlang")

packageVersion("rlang")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(tidyr)
})

# plot levels, dgp, cpi, ir
plot_levels <- function(data) {
  df_long <- data |>
    select(date, gdp, cpi, ir) |>
    pivot_longer(cols = c(gdp, cpi, ir), names_to = "series", values_to = "value")
  
  ggplot(df_long, aes(x = date, y = value)) +
    geom_line() +
    facet_wrap(~ series, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = NULL, y = NULL)
}

# pllot growth and ir level
plot_growth <- function(data) {
  df_long <- data |>
    select(date, gdp_growth, cpi_growth, ir) |>
    pivot_longer(cols = c(gdp_growth, cpi_growth, ir),
                 names_to = "series", values_to = "value")
  
  ggplot(df_long, aes(x = date, y = value)) +
    geom_line() +
    facet_wrap(~ series, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = NULL, y = NULL)
}

# prepare AR data
prepare_AR_data <- function(data, variable, p, p_max, intercept = TRUE) {
  stopifnot(p_max >= p, p >= 0)
  
  n <- nrow(data)
  start <- p_max + 1  
  
  y <- data[[variable]][start:n]
  y <- as.numeric(y)
  
  X <- NULL
  if (intercept) {
    X <- matrix(1, nrow = length(y), ncol = 1)
    colnames(X) <- "intercept"
  }
  
  if (p > 0) {
    for (lag in 1:p) {
      xlag <- data[[variable]][(start - lag):(n - lag)]
      xlag <- as.numeric(xlag)
      X <- cbind(X, xlag)
      colnames(X)[ncol(X)] <- paste0("lag_", lag)
    }
  }
  
  list(X = X, y = y)
}

# manual OLS
estimate_AR_model <- function(y, X) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  
  XtX_inv <- solve(t(X) %*% X)
  beta_hat <- XtX_inv %*% (t(X) %*% y)
  residuals <- y - X %*% beta_hat
  
  n <- length(y)
  k <- ncol(X)
  sigma_hat <- as.numeric(t(residuals) %*% residuals) / (n - k)
  
  se_beta <- sqrt(diag(sigma_hat * XtX_inv))
  t_stats <- as.numeric(beta_hat) / se_beta
  
  list(
    beta_hat = as.numeric(beta_hat),
    residuals = as.numeric(residuals),
    se_beta = se_beta,
    t_stats = t_stats
  )
}

# AIC and loglikelihood
calculate_AIC_loglik <- function(residuals, p, N) {
  dsigma_hat <- sum(residuals^2) / (N - p)
  dsigmalogdet <- log(dsigma_hat)
  dloglik <- -0.5 * (dsigmalogdet + (N - p) * (1 + log(2 * pi)))
  daic <- dsigmalogdet + (2.0 * p / N)
  list(loglik = dloglik, aic = daic, sigma_hat = dsigma_hat)
}

# run ar models
run_AR_models <- function(data, variable, max_lags = 5, intercept = TRUE) {
  N <- nrow(data)
  
  for (lag_p in 0:max_lags) {
    prep <- prepare_AR_data(data, variable, lag_p, max_lags, intercept)
    X <- prep$X
    y <- prep$y
    
    if (is.null(X) || ncol(as.matrix(X)) == 0) {
      stop("Design matrix X is empty. Use intercept=TRUE or p>0.")
    }
    
    est <- estimate_AR_model(y, X)
    aic <- calculate_AIC_loglik(est$residuals, lag_p, N)
    
    # results table
    var_names <- colnames(as.matrix(X))
    results_df <- data.frame(
      lag = rep(lag_p, length(var_names)),
      variable = ifelse(var_names == "intercept", "const",
                        paste0(variable, "_", gsub("lag_", "", var_names))),
      coef = est$beta_hat,
      se = est$se_beta,
      `t-value` = est$t_stats,
      row.names = NULL
    )
    
    cat("\n============================\n")
    cat(sprintf("AR(%d) for %s\n", lag_p, variable))
    print(results_df, row.names = FALSE)
    cat(sprintf("Log-likelihood: %.4f\n", aic$loglik))
    cat(sprintf("AIC: \t\t%.4f\n", aic$aic))
  }
}

# data
main <- function() {
  file_path <- "C:\\Users\\Sara\\Downloads\\Econometrics III\\data_assignment2.csv"
  data <- read_csv(file_path, show_col_types = FALSE)
  
  if (is.character(data$date) && grepl("Q", data$date[1])) {
    data$date <- as.Date(as.yearqtr(data$date, format = "%YQ%q"))
  } else {
    data$date <- as.Date(data$date)
  }
  
  total_quarters <- nrow(data)
  years_short <- 5
  quarters_short <- total_quarters - years_short * 4 - 1
  
  data <- data |>
    mutate(
      gdp_growth = 100 * (log(gdp) - lag(log(gdp))),
      cpi_growth = 100 * (log(cpi) - lag(log(cpi)))
    ) |>
    slice(-1)
  
  data_out_covid <- data |> slice(1:quarters_short)
  
  # plots
  print(plot_levels(data))
  print(plot_growth(data))
  
  # AR models
  max_lags <- 5
  intercept <- TRUE
  run_AR_models(data_out_covid, variable = "ir", max_lags = max_lags, intercept = intercept)
  
  return(data)
}

# store main in data
data <- main()


## assignment a
# building var design matrices
prepare_VAR_data <- function(data, vars, p, intercept = TRUE) {
  k <- length(vars)
  n <- nrow(data)
  start <- p + 1
  
  Y <- as.matrix(data[start:n, vars, drop = FALSE])
  
  X <- NULL
  if (intercept) {
    X <- matrix(1, nrow = nrow(Y), ncol = 1)
    colnames(X) <- "const"
  }
  
  # add lags of all variables
  for (lag in 1:p) {
    lag_block <- as.matrix(data[(start - lag):(n - lag), vars, drop = FALSE])
    colnames(lag_block) <- paste0(vars, "_L", lag)
    X <- cbind(X, lag_block)
  }
  
  list(Y = Y, X = X)
}

# estimating var
estimate_VAR <- function(Y, X) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  XtX_inv <- solve(t(X) %*% X)
  B_hat <- XtX_inv %*% (t(X) %*% Y)  # (m x k)
  
  U_hat <- Y - X %*% B_hat
  T <- nrow(Y)
  k <- ncol(Y)
  m <- ncol(X)
  
  Sigma_u <- (t(U_hat) %*% U_hat) / (T - m)
  
  list(B_hat = B_hat, U_hat = U_hat, Sigma_u = Sigma_u)
}

get_A_mats <- function(B_hat, vars, p, intercept = TRUE) {
  k <- length(vars)
  idx_start <- if (intercept) 2 else 1 
  A_list <- vector("list", p)
  
  for (lag in 1:p) {
    rows <- idx_start + (lag - 1) * k + (0:(k - 1))
    A_list[[lag]] <- t(B_hat[rows, , drop = FALSE])  
    colnames(A_list[[lag]]) <- vars
    rownames(A_list[[lag]]) <- vars
  }
  A_list
}

# companion matrix and check stability
companion_matrix <- function(A_list) {
  p <- length(A_list)
  k <- nrow(A_list[[1]])
  K <- k * p
  
  Fmat <- matrix(0, nrow = K, ncol = K)
  
  top <- do.call(cbind, A_list)
  Fmat[1:k, 1:(k*p)] <- top
  
  if (p > 1) {
    Fmat[(k + 1):K, 1:(K - k)] <- diag(K - k)
  }
  
  Fmat
}

check_stability <- function(A_list) {
  Fmat <- companion_matrix(A_list)
  eigvals <- eigen(Fmat, only.values = TRUE)$values
  mod <- Mod(eigvals)
  list(
    eigvals = eigvals,
    max_modulus = max(mod),
    stable = all(mod < 1)
  )
}

# run var on gdp and inflation
run_part_a_VAR3 <- function(data) {
  vars <- c("gdp_growth", "cpi_growth")
  p <- 3
  intercept <- TRUE
  
  prep <- prepare_VAR_data(data, vars, p, intercept)
  est <- estimate_VAR(prep$Y, prep$X)
  
  cat("\n============================\n")
  cat("Question a: estimated VAR(3) with constant\n")
  cat("Variables:", paste(vars, collapse = ", "), "\n\n")
  
  B <- est$B_hat
  rownames(B) <- colnames(prep$X)
  colnames(B) <- vars
  print(round(B, 6))
  
  # stability?
  A_list <- get_A_mats(est$B_hat, vars, p, intercept)
  stab <- check_stability(A_list)
  
  cat("\nStability check (companion eigenvalues):\n")
  cat("Max |eigenvalue| =", round(stab$max_modulus, 6), "\n")
  cat("Stable? ", ifelse(stab$stable, "YES (all |λ| < 1)", "NO (some |λ| >= 1)"), "\n")
}

run_part_a_VAR3(data)

#assignment b ----

aic_calc <- function(sigma_u_LS, T_eff, p, vars_num){
  #sigma_u is the unbiased estimated variance matrix of u 
  
  #obtain the MLE residual variance matrix
  sigma_u_MLE <- ((T_eff - (1 + vars_num * p))/T_eff ) * sigma_u_LS 
  
  #Log-L
  logL <- ( -(T_eff/2) * vars_num ) * (1 + log(2*pi)) - (T_eff/2) * log(det(sigma_u_MLE))
  
  #aic
  param_num <- vars_num * (1 + vars_num * p)
  aic <- 2 * param_num - 2 * logL
  
  return(list(logL = logL, aic = aic))
  
} 

run_part_b_VAR_p_selection <- function(data, short = TRUE){
  #specify candidate models to compare
  p_comp <- 5
  #variable specification
  vars <- c("gdp_growth", "cpi_growth")
  #include intercept
  intercept <- TRUE
  
  candidates_list <- c()
  logL_list <- c()
  aic_list <- c()
  stab_list <- c()
  
  for (p in 1:p_comp){
    
    candidates_list <- append(candidates_list, paste0("VAR(", p,")"))
    
    #data alignment to ensure comparison
    if (p < p_comp){
      data_trimmed <- data[-(1:(p_comp-p)),]
    }else{
      data_trimmed <- data 
    }
    
    prep <- prepare_VAR_data(data_trimmed, vars, p, intercept) 
    T_eff <- nrow(prep$X)
    
    #estimation
    est <- estimate_VAR(prep$Y, prep$X)
    
    #parameter estimates
    B <- est$B_hat
    rownames(B) <- colnames(prep$X)
    colnames(B) <- vars
    
    #check stability using 
    A_list <- get_A_mats(B, vars, p, intercept)
    stab <- check_stability(A_list)
    stab_list <- append(stab_list, ifelse(stab$stable, "Stable", "Non-stable"))
    
    ll <- aic_calc(est$Sigma_u, T_eff, p, length(vars))
    logL <- ll$logL
    aic <- ll$aic
    logL_list <- append(logL_list, logL)
    aic_list <- append(aic_list, aic)
    
    aic_adj <- aic/T_eff
    
    if (short == FALSE){
      cat("\n============================\n")
      cat("VAR(", p, ") parameter estimates:\n")
      print(round(B, 6))
      cat("\nStability check (companion eigenvalues):\n")
      cat("Max |eigenvalue| =", round(stab$max_modulus, 6), "\n")
      cat("VAR(", p, ")", "Stable?", ifelse(stab$stable, "YES (all |λ| < 1)", "NO (some |λ| >= 1)"), "\n")
      cat("log likelihood:", logL, "AIC:", aic)
    }
  }
  
  cat("\n============================\n")
  cat("Candidate Comparison\n")
  #concise results table
  table_b <- data.frame(
    Candidates = candidates_list,
    log_likelihood = logL_list,
    AIC = aic_list,
    Stability = stab_list
  )
  print(table_b)
  optimal_p <<- which.min(table_b$AIC)
  cat("The Model with the smallest AIC:", table_b$Candidates[optimal_p])
  
  
  
}

#full results with estimates and stability + results
run_part_b_VAR_p_selection(data, short = FALSE)
#only AIC comparison table
run_part_b_VAR_p_selection(data, short = TRUE)

# assignment b (pckg) ----
install.packages("vars")
library(vars)

data_pkg <- as.data.frame(data)
data_pkg <- data_pkg[, -c(1:4)]
data_ts <- ts(data_pkg, start=1, frequency=1)

p_comp <- 5
for (p in 1:5){
  cat("\n==========\n")
  cat("lag p:", p)
  if (p < p_comp){
    data_trimmed <- data_ts[-(1:(5-p)),]
  }else{
    data_trimmed <- data_ts 
  }
  var_model <- VAR(data_trimmed, p=3, type="const")
  ll <- as.numeric(logLik(var_model))
  cat("\nLL:", ll)
  cat("\nAIC:", AIC(var_model))
  
}




#assignment c ----

irf_function <- function(data, vars, p, irf_periods){
  
  
}

#

#variable specification
vars <- c("gdp_growth", "cpi_growth")
#include intercept
intercept <- TRUE

#selection matrix J
k <- length(vars)
K <- k * optimal_p
J <- matrix(0, nrow = k, ncol = K)
J[1:k, 1:k] <- diag(k)

#Var(3) estimation
prep <- prepare_VAR_data(data, vars, optimal_p, intercept) 
est <- estimate_VAR(prep$Y, prep$X)
B <- est$B_hat

#companion matrix and eigen decomposition
A_list <- get_A_mats(B, vars, optimal_p, intercept) #phi matrices
psi <- companion_matrix(A_list)
eigen_psi <- eigen(psi)
psi_egvc <- eigen_psi$vectors
psi_egvl_diag <- diag(eigen_psi$values)

irf_period <- 10
irf_list <- vector("list", irf_period + 1)
#contains all impulses of all cross-elements of each period

#for each irf_period
for (i in 0:irf_period){
  if (i == 0){
    psi_power_i <- diag(nrow(psi_egvc))
  }else{
    psi_power_i <- psi_egvc %*% (psi_egvl_diag^i) %*% solve(psi_egvc) 
  }
  irf_i <- J %*% psi_power_i %*% t(J)
  irf_list[[i+1]] <- irf_i
}

#to individual irf of each cross-element
irf_decomp <- list()
for (m in 1:k){
  for (n in 1:k){
    responses <- c()
    for (i in 1:(irf_period + 1)){
      responses <- append(responses, irf_list[[i]][m,n])
    }
    irf_decomp <- append(irf_decomp, list(responses))
  }
}

#IRF Plots
period_indicators <- seq(0, irf_period, 1)

par(mfrow = c(2,2), mar = c(4,4,2,1), oma = c(4,4,4,2))
plot(period_indicators, irf_decomp[[1]], type = "l", ylab = "", xlab = "") 
plot(period_indicators, irf_decomp[[2]], type = "l", ylab = "", xlab = "") 
plot(period_indicators, irf_decomp[[3]], type = "l", ylab = "", xlab = "") 
plot(period_indicators, irf_decomp[[4]], type = "l", ylab = "", xlab = "") 
mtext("GDP Growth Impulse", side = 3, line = 2.5, at = 0.25, outer = TRUE)
mtext("CPI Growth Impluse",  side = 3, line = 2.5, at = 0.75, outer = TRUE)
mtext("GDP Growth Response", side = 2, line = 2.5, at = 0.75, outer = TRUE)
mtext("CPI Growth Response", side = 2, line = 2.5, at = 0.25, outer = TRUE)
par(mfrow = c(1,1))


# assignment c (package) ----

var_3 <- VAR(data_ts, p = 3, type="const")

#impluse, gdp_gr on gdp_gr
irf_result <- irf(var_3,
                  impulse = "gdp_growth",
                  response = "gdp_growth",
                  n.ahead = 10,
                  ortho = FALSE)
plot(period_indicators, irf_decomp[[1]], type = "l") 
plot(irf_result, plot.type = "single", ylab = "Response", xlab = "Horizon")
#impluse, cpi_gr on gdp_gr
plot(period_indicators, irf_decomp[[2]], type = "l") 
irf_result <- irf(var_3,
                  impulse = "cpi_growth",
                  response = "gdp_growth",
                  n.ahead = 10,
                  ortho = FALSE)
plot(irf_result, plot.type = "single", ylab = "Response", xlab = "Horizon")
#impluse, gdp_gr on cpi_gr
plot(period_indicators, irf_decomp[[3]], type = "l") 
irf_result <- irf(var_3,
                  impulse = "gdp_growth",
                  response = "cpi_growth",
                  n.ahead = 10,
                  ortho = FALSE)
plot(irf_result, plot.type = "single", ylab = "Response", xlab = "Horizon")
#impluse, cpi_gr on cpi_gr
plot(period_indicators, irf_decomp[[4]], type = "l") 
irf_result <- irf(var_3,
                  impulse = "cpi_growth",
                  response = "cpi_growth",
                  n.ahead = 10,
                  ortho = FALSE)
plot(irf_result, plot.type = "single", ylab = "Response", xlab = "Horizon")





# assignment d ----

data_pre_covid <- data |> slice(1:(nrow(data) - 20))

# redo part b: lag selection on pre-COVID data
run_part_b_VAR_p_selection(data_pre_covid, short = TRUE)

# variable specification
vars <- c("gdp_growth", "cpi_growth")
intercept <- TRUE

# selection matrix J
k <- length(vars)
K <- k * optimal_p
J <- matrix(0, nrow = k, ncol = K)
J[1:k, 1:k] <- diag(k)

# VAR estimation on pre-COVID data
prep <- prepare_VAR_data(data_pre_covid, vars, optimal_p, intercept)
est <- estimate_VAR(prep$Y, prep$X)
B <- est$B_hat

# companion matrix and eigen decomposition
A_list <- get_A_mats(B, vars, optimal_p, intercept)
psi <- companion_matrix(A_list)
eigen_psi <- eigen(psi)
psi_egvc <- eigen_psi$vectors
psi_egvl_diag <- diag(eigen_psi$values)

irf_period <- 10
irf_list <- vector("list", irf_period + 1)

# for each irf period
for (i in 0:irf_period){
  if (i == 0){
    psi_power_i <- diag(nrow(psi_egvc))
  }else{
    psi_power_i <- psi_egvc %*% (psi_egvl_diag^i) %*% solve(psi_egvc)
  }
  irf_i <- J %*% psi_power_i %*% t(J)
  irf_list[[i+1]] <- irf_i
}

# to individual irf of each cross-element
irf_decomp_precovid <- list()
for (m in 1:k){
  for (n in 1:k){
    responses <- c()
    for (i in 1:(irf_period + 1)){
      responses <- append(responses, irf_list[[i]][m,n])
    }
    irf_decomp_precovid <- append(irf_decomp_precovid, list(responses))
  }
}

# IRF plots
period_indicators <- seq(0, irf_period, 1)

par(mfrow = c(2,2), mar = c(4,4,2,1), oma = c(4,4,4,2))
plot(period_indicators, irf_decomp_precovid[[1]], type = "l", ylab = "", xlab = "")
plot(period_indicators, irf_decomp_precovid[[2]], type = "l", ylab = "", xlab = "")
plot(period_indicators, irf_decomp_precovid[[3]], type = "l", ylab = "", xlab = "")
plot(period_indicators, irf_decomp_precovid[[4]], type = "l", ylab = "", xlab = "")
mtext("GDP Growth Impulse", side = 3, line = 2.5, at = 0.25, outer = TRUE)
mtext("CPI Growth Impulse",  side = 3, line = 2.5, at = 0.75, outer = TRUE)
mtext("GDP Growth Response", side = 2, line = 2.5, at = 0.75, outer = TRUE)
mtext("CPI Growth Response", side = 2, line = 2.5, at = 0.25, outer = TRUE)
par(mfrow = c(1,1))





# assignment e ----

run_part_e_VAR2 <- function(data){
  
  # extend the VAR system to three dimensions by adding the interest rate
  vars <- c("gdp_growth", "cpi_growth", "ir")
  
  p <- 2
  intercept <- TRUE
  
  prep <- prepare_VAR_data(data, vars, p, intercept)
  est <- estimate_VAR(prep$Y, prep$X)
  
  cat("\n============================\n")
  cat("Question e: estimated VAR(2) with constant\n")
  cat("Variables:", paste(vars, collapse = ", "), "\n\n")
  
  B <- est$B_hat
  rownames(B) <- colnames(prep$X)
  colnames(B) <- vars
  print(round(B, 6))
}

run_part_e_VAR2(data)

# pre-COVID version ----

data_pre_covid <- data[1:(nrow(data) - 20), ]

run_part_e_VAR2(data_pre_covid)