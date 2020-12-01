generate_init_condi <- function(r0,
                                Di = 2.9,
                                Dp = 2.3,
                                De = 2.9,
                                Dq = c(7, 7, 7, 7, 7),
                                alpha = 0.55,
                                Dh = 8.6,
                                N = 57779622,
                                flowN = c(0, 0, 0, 0, 0)
) {
  
  stopifnot(r0>=0 & r0<=1 & Di>=0 & Dp>=0 & De>=0 & all(Dp>=0) & alpha>=0 & alpha<=1 & Dh>=0 & N>=0 & all(flowN>=0))
  
  ## N            : population size
  ## H0           : initial number of hospitalized cases based on the reports
  ## R0           : initial number of removed individuals
  ## De           : latent period
  ## r0           : initial ascertainment rate
  ## realData     : real data from the CDC
  load('base.RData')
  start_point=which(SA$Date==as.Date('2020/03/15'))
  propA0 <- (1-r0)/(r0)
  E0 <- (sum(cases_SA[(start_point+2):(start_point+4)])) * (1 + propA0)
  P0 <- (sum(cases_SA[start_point:(start_point+1)])) * (1 + propA0)
  R0 <- SA[start_point,3]
  D0 <- SA[start_point,2]
  H0 <- SA[start_point-7,1]*0.5
  if(R0+D0+H0>=SA[start_point,1]){
    I0=0
  }else{
    I0 <- SA[start_point,1]-H0-R0-D0
  }
  A0 <- I0 * propA0
  S0 <- N - E0 - I0 - A0 - H0 - R0 - D0
  yinit <- c(S = S0, E = E0, P = P0, A = A0,I = I0, H = H0, R = R0,  D = D0)
  time=seq.Date(from=min(SA$Date),to=as.Date("2020/07/31"),by='days')
  data_SA=as.data.frame(time[c(start_point:length(time))])
  colnames(data_SA) <- 'OnsetDate'
  data_SA$CaseNum=cases_SA[c(start_point:length(time))]
  data_SA$CaseSum=SA$Confirmed[c(start_point:length(time))]
  realData_all <- data_SA
  rownames(realData_all) <- data_SA$OnsetDate
  realData_all <- realData_all[,-1]
  realData <- realData_all 
  
  #cases
  daily_new_case <- cases_SA[c(start_point:length(time))]
  daily_new_case_all <- realData[, 1]
  ##
  init_states <- round(c(S = S0, E = E0, P = P0, I = I0, A = A0, H = H0, R = R0+D0), 0)
  
  ## helper function
  # transform variables to a form that SEIRpred can use
  # so that SEIRpred can be re-used as much as possible
  transform_var_main_5stage=function(pars) {
    b_vec <- pars[1:4]
    b_vec <- c(b_vec[1], b_vec[2:4])
    ##
    r12 <- pars[5]
    r3 <- 1 / (1 + (1 - r12) / (r12 * exp(pars[6])))
    r4 <- 1 / (1 + (1 - r3) / (r3 * exp(pars[7])))
    r5 <- 1 / (1 + (1 - r4) / (r4 * exp(pars[8])))
    r_vec <- c(r12,r3,r4,r5)
    
    return(list(b_vec, r_vec))
  }
  
  return(list(Di=Di,
              Dp=Dp,
              De=De,
              Dq=Dq,
              alpha=alpha,
              Dh=Dh,
              N=N,
              flowN=flowN,
              daily_new_case = daily_new_case, 
              daily_new_case_all = daily_new_case_all, 
              init_states = init_states,
              days_to_fit=1:170,
              stage_intervals=list(
                c(start=1, end=12),
                c(start=13, end=47),
                c(start=48, end=78),
                c(start=79, end=170)
              ),
              var_trans_fun=transform_var_main_5stage,
              par_lower = c(b12 = 0, b3 = 0, b4 = 0, b5 = 0, r12 = 0, delta3 = -10, delta4 = -10, delta5 = -10),
              par_upper = c(b12 = 2, b3 = 2, b4 = 2, b5 = 2, r12 = 1, delta3 = 10, delta4 = 10, delta5 = 10)))
  # TODO: please confirm the following:
  # boundaries for delta3-5 will not be used, they are here merely to meet the formality imposed by runMCMC
}

# get_init_sets_list is an alias of generate_init_condi in order not to break exsiting code
get_init_sets_list = generate_init_condi

delta_mean <- 0
delta_sd <- 1
beta_shape1 <- 10
beta_shape2 <- 90
