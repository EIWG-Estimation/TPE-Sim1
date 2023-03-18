#' add_discontinuation
#' @description Add a column which indicated whether a patient 
#' is on or off their assigned treatment 
#' @param data data frame; output from simulate_ontrt_data
#' @param beta_discont data frame; parameters of logistric model for discontinuation
add_discontinuation <- function(data, beta_discont){
  
  # Indicator whether inputs include beta coef for current response
  is_dnar <- any(names(beta_discont) == "beta_current_resp")
  
  # 1. Merge beta parameter for logistic discontinuation model with the data
  # 2. create column with response from previous visit
  # 3. Calculate the logit of the distcont probability
  # 4. Calculate discontinuation probability (p=0 for visit = 0)
  # 5. Simulate on-trt indicator
  #     - trt-discont is not intermittent
  data_discont <- full_join(data, beta_discont, 
                            by = c("group", "visitn")) 
  rm("data")
  data_discont <- data_discont %>% 
    group_by(sim_run, subjid) %>% 
    mutate(prev_response = lag(response_ontrt, default = 0),
           logit_prob = beta_main + beta_prev_resp * prev_response + beta_baseline * baseline_var,
           logit_prob = if(is_dnar){logit_prob + beta_current_resp * response_ontrt} else {logit_prob},
           prob = 1/(1+exp(-logit_prob)),
           prob = if_else(is.na(prob), 0, prob),
           ontrt = 1-rbinom(n = n(), size = 1, prob = prob)) %>%  
    mutate(first_offtrt_visit = min(visitn[ontrt == 0], max(visitn)),
           ontrt = if_else(visitn > first_offtrt_visit, 0, ontrt)) %>% 
    select(-contains("beta_"), -contains("prob"),
           -first_offtrt_visit, -prev_response)
  
  return(data_discont)  
}