#' simulate_ontrt_data
#' @description Creates normal random numbers emulating a two-arm trial where 
#' every patient continues to stay on their randomized treatment
#' @param mu_trt numeric; vector of visit means in trt group
#' @param mu_ctl numeric; vector of visit means in ctl group
#' @param covar_trt matrix; covariance matrix of single subjects in trt group
#' @param covar_ctl matrix; covariance matrix of single subjects in ctl group
#' @param n_trt numeric; number of subjects in trt group
#' @param n_ctl numeric; number of subjects in ctl group
#' @param n_sim numeric; number of simulation runs
#' @param check_input logical; if TRUE, input arguments will be checked for consistency
#' @import tidyverse
#' @importFrom mvtnorm
simulate_ontrt_data <- function(mu_trt, mu_ctl,
                                covar_trt, covar_ctl,
                                n_trt, n_ctl,
                                n_sim,
                                check_input) {
  
  
  # Number of visit
  n_visit <- length(mu_trt)
  # Total sample size
  n_total <- n_trt + n_ctl
  
  # Consistency check input arguments 
  if (check_input) {
    if (n_visit != length(mu_ctl)) stop("Number of visits much be identical between groups")
    if (n_visit != dim(covar_trt)[1]) stop("Input for trt is not consistent")
    if (n_visit != dim(covar_ctl)[1]) stop("Input for ctl is not consistent")
    if (n_visit != dim(covar_trt)[1]) stop("Input for trt is not consistent")
    if (dim(covar_trt)[2] != dim(covar_trt)[1]) stop("covar_trt is not quadratic")
    if (dim(covar_ctl)[2] != dim(covar_ctl)[1]) stop("covar_ctl is not quadratic")
  }
  
  
  # Create separate vectors with simulated responses for trt and ctl
  # Responses are simulated for all patients in all simulation runs
  # Responses for a single subject are on consecutive positions
  # baseline: visit=0 
  reponse_trt_mat <- mvtnorm::rmvnorm(n = n_trt*n_sim, mean = mu_trt, sigma = covar_trt)
  data_trt <- expand_grid( sim_run = 1:n_sim,
                           subjid = seq.int(1, n_trt),
                           data.frame(visitn = 0:(n_visit-1), 
                                      mean_trt = mu_trt, 
                                      mean_ctl = mu_ctl)) %>% 
    mutate(response_ontrt = as.vector(t(reponse_trt_mat)), 
           group = "trt")
  rm("reponse_trt_mat")
  
  reponse_ctl_mat <- mvtnorm::rmvnorm(n = n_ctl*n_sim, mean = mu_ctl, sigma = covar_ctl)
  data_ctl <- expand_grid( sim_run = 1:n_sim,
                           subjid = seq.int(n_trt+1, n_total),
                           data.frame(visitn = 0:(n_visit-1), 
                                      mean_trt = mu_trt, 
                                      mean_ctl = mu_ctl)) %>% 
    mutate(response_ontrt = as.vector(t(reponse_ctl_mat)), 
           group = "ctl")
  rm("reponse_ctl_mat")
  
  # Create output data frame
  # First part of the data frame is for subject on trt, second for subjects on ctl
  simulated_data <- bind_rows(data_trt, data_ctl)
  
  # Add separate column with baseline response
  simulated_data <- simulated_data %>% 
    group_by(sim_run, subjid) %>% 
    mutate(baseline_var = response_ontrt[visitn == 0])
  
  return(simulated_data)
}