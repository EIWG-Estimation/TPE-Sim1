# --------------------------------------------------------------
#  Date    : 18-03-23
#  Authors : Thomas Drury (original Tobias Muetze)
#  History : Original code simulate_data.R
#  Purpose : Creates simulations for DNAR_DROP scenario 
#       
# --------------------------------------------------------------
#  Notes:
#  
# 1. Source R code and load libraries
# 2. Define general simulation study parameters
# 3. Simulate on-trt data
# 4. Simulate off-trt data
# 5. Simulate missing data
#
# --------------------------------------------------------------

# --------------------------------------------------------------
# 1. Source R code and load libraries 
# --------------------------------------------------------------

setwd("/shared-scratch/area/tad66240/TPE-Sim1/")
rm(list = ls())

starttime = Sys.time()

set.seed(12345)
library(tidyverse)

r_func_vec = c("add_discontinuation", "add_missingness", "add_response_offtrt", "simulate_ontrt_data")

for(i in seq_along(r_func_vec)) {
  source(paste0("code/functions/", r_func_vec[i], ".R"))
}

# ---------------------------------------------------------------
# 2. Define general simulation study parameters
# ---------------------------------------------------------------

scenario = "dnar_drop"

# Simulations and number of subjects 
n_sim = 10000
n_ctl = 200
n_trt = 200


# Visits
visits  = c(0, 4, 8, 14, 20, 26)
n_visit = length(visits)


# Spatial power correlation matrix - corresponds to AR(1) for equidistant visits
rho      = 0.8
exponent = abs(matrix(visits, nrow = n_visit, ncol = n_visit, byrow = TRUE) - visits) / (visits[2] - visits[1])
corr_mat = rho^exponent


# Variance at baseline is pooled variance from each group
sigma2_bl = 0.48


# Variance for treatment and control group at each visit
sigma2_trt = c(sigma2_bl, 0.75, 0.8, 0.9, 1.06, 1.14)
sigma2_ctl = c(sigma2_bl, 0.8, 1.1, 1.4, 1.23, 1.48)


# Percentage of pts. that discontinued randomized trt
rate_trtdiscont_trt = (175 - c(175, 169, 164, 156, 151, 149)) / 175 
rate_trtdiscont_ctl = (178 - c(178, 172, 165, 158, 140, 133)) / 178 


# On-treatment mean for treatment
mu_trt = c(7.92, 7.55, 7.20, 7.1, 7.05, 7.05)
mu_ctl = c(7.92, 7.82, 7.8, 7.8, 7.78, 7.78)


# covariance matrix for treatment group
mat1 = matrix(sigma2_trt, nrow = n_visit, ncol = n_visit, byrow = TRUE)
mat2 = matrix(sigma2_trt, nrow = n_visit, ncol = n_visit, byrow = FALSE)
covar_trt = sqrt(mat1 * mat2) * corr_mat
rm(list = c("mat1", "mat2"))


# covariance matrix for control group
mat1 = matrix(sigma2_ctl, nrow = n_visit, ncol = n_visit, byrow = TRUE)
mat2 = matrix(sigma2_ctl, nrow = n_visit, ncol = n_visit, byrow = FALSE)
covar_ctl = sqrt(mat1 * mat2) * corr_mat
rm(list = c("mat1", "mat2"))


# Beta parameter for logistic regression modeling of DNAR trt discont
beta_discont_dnar = data.frame(group = rep(c("ctl", "trt"), each = 5),
                               visitn = rep(1:5, times = 2),
                               beta_main = rep(-21, times = 10),
                               beta_prev_resp = c(1.42, 1.14, 1.33, 1.51, 1.46,
                                                  1.42, 1.14, 1.47, 1.48, 1.40),
                               beta_current_resp = c(0.69, 0.69, 0.69, 0.72, 0.72,
                                                     0.72, 0.75, 0.77, 0.77, 0.76),
                               beta_baseline = c(0, 0.3, 0.1, 0.05, 0,
                                                 0, 0.3, 0.1, 0.05, 0),
                               stringsAsFactors = FALSE)


# ---------------------------------------------------------------
# 3. Simulate on-trt data
# ---------------------------------------------------------------

df_data_ontrt = simulate_ontrt_data(mu_trt = mu_trt,
                                    mu_ctl = mu_ctl, 
                                    covar_trt = covar_trt,
                                    covar_ctl = covar_ctl, 
                                    n_trt = n_trt, 
                                    n_ctl = n_ctl, 
                                    n_sim = n_sim, 
                                    check_input = TRUE)


# ---------------------------------------------------------------
# 4. Simulate off-treatment data
# ---------------------------------------------------------------

# Add discontinuation indicators to data frame
df_data_offtrt_dnar = add_discontinuation(data = df_data_ontrt, beta_discont = beta_discont_dnar)
rm("df_data_ontrt")


# Drop scenario means immediate change
df_data_offtrt_dnar = df_data_offtrt_dnar %>%
  ungroup %>%
  mutate(response_offtrt = case_when(ontrt == 1 ~ response_ontrt,
                                     ontrt == 0 & group == "ctl" ~ response_ontrt - 0.6,
                                     ontrt == 0 & group == "trt" ~ response_ontrt - 0.2))


# ---------------------------------------------------------------
# 5. Simulate missing data
# ---------------------------------------------------------------

p_miss_vec = seq(0.1, 0.6, by = 0.1)

for(i in seq_along(p_miss_vec)) {

  p_miss = p_miss_vec[i]
  theta  = data.frame(group = c("ctl", "trt"),
                      theta1 = c(log(p_miss/(1-p_miss)), log(p_miss/(1-p_miss))),
                      theta2 = c(0, 0), 
                      stringsAsFactors = FALSE)
  
  df_final_dnar = add_missingness(data = df_data_offtrt_dnar, theta = theta)
  
  file_name = paste0("/shared/259/arenv/arwork/tad66240/", scenario, "_p", str_replace(string = as.character(p_miss), pattern = "\\.", replacement = ""), ".csv")
  write_csv(df_final_dnar, file = file_name)
  
  rm("df_final_dnar")
  
}

Sys.time() - starttime







