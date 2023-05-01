#' add_response_offtrt
#' @description Add a column with the off-treatment data 
#' @param data data frame; output from add_discontinuation
#' @param trt_offtrt_effect character; describes how patient in treatment group 
#' loose effect after trt-distcont
#' @param gamma numeric; For linear loss of trt effect, 
#' gamma is the proportion that is lost between visits.
#' For instant drop, gamma is the proportion to which the trt effect drops
add_response_offtrt <- function(data, trt_offtrt_effect = c("linear", "instantdrop"), gamma) {
  
  if (!(trt_offtrt_effect %in% c("linear", "instantdrop"))) {
    stop("trt_offtrt_effect must be either 'linear' or 'instantdrop'.")
  }
  # Check that gamma is defined
  if (is.null(gamma)) {
    stop("gamma must be defined")
  }
  if (!is.numeric(gamma) | (gamma < 0) | (gamma > 1)) {
    stop("gamma must be between 0 and 1")
  }
  
  # Separate data frames for treatment and control
  data_trt <- data %>% 
    filter(group == 'trt')
  data_ctl <- data %>% 
    filter(group == 'ctl')
  rm("data")
  
  # For control off-treatment data is on-treatment data
  data_ctl$response_offtrt <- data_ctl$response_ontrt
  
  # Create response for drop to gamma*(initial effect) for those off-treatment
  if (trt_offtrt_effect == "instantdrop") {
    data_trt <- data_trt %>% 
      mutate(response_offtrt = if_else(ontrt == 0, 
                                       response_ontrt + gamma*(mean_ctl - mean_trt), 
                                       response_ontrt))
  }
  # Create off-treatment data for linear leveling off treatment effect
  if (trt_offtrt_effect == "linear") {
    data_trt <- data_trt %>% 
      group_by(sim_run, subjid) %>% 
      mutate(vists_since_trtdisc = cumsum(ontrt == 0),
             factor_trt_loss = pmin(vists_since_trtdisc*gamma, 1),
             response_offtrt = if_else(ontrt == 0, 
                                       response_ontrt + factor_trt_loss*(mean_ctl - mean_trt), 
                                       response_ontrt)) %>% 
      select(-vists_since_trtdisc, -factor_trt_loss)
  }
  
  return(data_trt %>% bind_rows(data_ctl))
}