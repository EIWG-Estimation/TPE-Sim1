#' add_missingness
#' @description Add a column with missing data indicator and response with NA 
#' @param data data frame; output from add_discontinuation
#' @param theta data frame; parameters of logistic model for missingness
add_missingness <- function(data, theta) {
  
  # Merge theta parameters with data set
  data_missingness <- full_join(data, theta, 
                                by = "group")
  rm("data")
  
  # Add column with missing data indicator
  # 1. Count off-treatment visits
  # 2. Calculate logit of missingness probability
  # 3. Calculate missingness probability
  #   - missingness prob is zero when patient is on-trt
  # 4. Simulate missingess indicator
  #   - missingness is not intermittent
  data_missingness <- data_missingness %>% 
    group_by(sim_run, subjid) %>% 
    mutate(vists_since_trtdisc = cumsum(ontrt == 0),
           logit_prob = theta1 + theta2 * vists_since_trtdisc,
           prob = 1/(1+exp(-logit_prob)),
           prob = if_else(ontrt == 1, 0, prob),
           is_missing = as.numeric(rbinom(n = n(), size = 1, prob = prob)),
           first_missing_visit = min(visitn[is_missing == 1], max(visitn)),
           is_missing = if_else(visitn > first_missing_visit, 1, is_missing),
           response = if_else(is_missing == 1, NA_real_, response_offtrt)) %>% 
    select(-contains("theta"), -contains("prob"),
           -first_missing_visit, -vists_since_trtdisc)
  
  return(data_missingness)
}