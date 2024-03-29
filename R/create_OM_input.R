#' @title create_OM_input
#'
#' @description function to create life history parameters for simulation model (without temporal or spatial components)
#'
#' Life history parameters
#' @param Amax maximum age
#' @param Linf von Bertalanffy asymptotic length
#' @param vbk von Bertalanffy growth coefficient
#' @param t0 default = -0.01; von Bertalanffy theoretical age at length = 0
#' @param M natural mortality
#' @param LWa length-weight a parameter
#' @param LWb length-weight b parameter
#' @param Lmat length at 50% maturity (for fecundity calculation)
#' @param agefec default = NULL; fecundity-at-age vector
#' @param h steepness parameter
#' @param cstar default = 1; logical parameter for whether larval are assumed to be capable of locating recruitment habitat (1) or not (0)
#' @param binwidth default = 1; width of length bins
#' @param start_ages default = 0; age to start (either 0 or 1)
#'
#' Fishery parameters
#' @param q catchability
#' @param S1 first selectivity parameter for double normal 
#' @param S2 second selectivity parameter for double normal
#' @param S3 third selectivity parameter for double normal
#' @param S4 fourth selectivity parameter for double normal
#' @param S5 fifth selectivity parameter for double normal
#' @param S6 sixth selectivity parameter for double normal
#' @param sizemin minimum size limit - for retention calculation
#' 
#' Dynamic effort parameters 
#' @param persis default = 1; stickyness parameter for effort dynamics (seen in van Poorten and Camp 2019)
#' @param sig1e shape parameter that describes how strongly fishing effort responds to changes in total utility
#' @param max_eff total maximum effort (obtain from data)
#' @param w_cost power weight for fisher cost function
#'
#' Utility parameters - for recreational fisheries
#' parameter order for "three_param": max PWU, steepness, inflection point
#' @param Ut_cpue utility for CPUE
#' @param Ut_hpue utility for harvest
#' @param Ut_size utility for size
#' @param Ut_dist dis-utility for distance
#' @param Ut_crowd dis-utility for crowding
#'
#' Preference movement parameters
#' @param w_dep power weight for depth preference
#' @param w_hab power weight for habitat preference
#' @param w_dist power weight for fish movement distance
#' 
#' AR parameters
#' @param AR_flag designed to be a flag for whether or not there artificial reefs are implemented
#' @param AR_prop_M percent change at artificial reef sites
#' @param AR_prop_rec percent change of recruitment habitat from artificial reefs (changes ahab and bhab - habitat specific a and b parameters for recruitment)
#' @param AR_prop_q percent change of catchability at artificial reef sites
#' @param AR_prop_U percent change of utility of artificial reefs (how much fisher utilty is improved by artificial reefs inherently)
#' NR parameters
#' @param NR_flag designed to be a flag for whether or not there are natural reefs
#' @param NR_prop_M percent change at natural reef sites
#' @param NR_prop_rec percent change of recruitment habitat from natural reefs (changes ahab and bhab - habitat specific a and b parameters for recruitment)
#' @param NR_prop_q percent change of catchability at natural reef sites
#' @param NR_prop_U percent change of utility of natural reefs (how much fisher utilty is improved by natural reefs inherently)
#'
#' @param nseason specify number of sub-time periods per year; default=1 
#'

#################################################
## Input list for create_OM_input R function ####
#################################################
input_list <- list(Amax = 20, # SEDAR 52
                  Linf = 85.6, # SEDAR 52
                  vbk = 0.1919, # SEDAR 52
                  t0 = -0.3925, # SEDAR 52
                  M = c(2,1.2,0.19,0.15,0.129,0.115,0.106,0.099,0.095,0.091,0.088,0.086,0.085,0.083,0.082,0.081,0.081,0.08,0.08,0.079,0.078), # SEDAR 52
                  # M = 0.08, # SEDAR 52
                  LWa = 1.673e-5, # SEDAR 52
                  LWb = 2.953, # SEDAR 52
                  Lmat = 35, # Glenn et al., 2017
                  agefec = c(0, 0, 0.35E6, 2.62E6, 9.07E6, 20.3E6, 34.71E6, 49.95E6, 64.27E6, 76.76E6, 87.15E6, 95.53E6, 102.15E6, 107.3E6, 111.27E6, 114.3E6, 116.61E6, 118.36E6, 119.68E6, 120.67E6, 123.234591E6)*1e-3, # SEDAR 52
                  # agefec = NULL, # SEDAR 52
                  h = 0.99, # SEDAR 52
                  cstar = 1,
                  binwidth = 1,
                  start_ages = 0,
                  q = 0.00083, # 
                  S1 = 3.19257,
                  S2 = -0.84110, 
                  S3 = -0.07316,
                  S4 = -0.21389,
                  S5 = -12.43890,
                  S6 = -4.45326,
                  persis = 0.5,
                  sig1e = 2,
                  sizemin = 40.64, # 16" regulation
                  max_eff = 1069073.6, # calculated from data
                  w_cost = 1,
                  Ut_cpue = c(1.5,1,3.855415), # tuned
                  Ut_hpue = c(1.5,2,1.435766), # tuned
                  Ut_size = c(1.5,0.15,42.70597), # tuned
                  Ut_dist = c(1.5,-0.006,185.4731), # tuned
                  Ut_crowd =  c(1.5,-0.0003,2333.956), # tuned
                  w_dep = 1,
                  w_hab = 1,
                  w_dist = 1,
                  AR_flag = 0,
                  AR_prop_M = 0,
                  AR_prop_rec = 0,
                  AR_prop_q = 0,
                  AR_prop_U = 0,
                  NR_flag = 1,
                  NR_prop_M = -0.2,
                  NR_prop_rec = 0.2,
                  NR_prop_q = 0.2,
                  NR_prop_U = 0.2,
                  nseason = 1)


create_OM_input <-
function(input_list)
{
  with(input_list, {
    ## Length at age
    # length bins
    mids <- seq((binwidth/2), Linf*1.3, by = binwidth) # midlength bins
    highs <- mids+(binwidth/2) # length bins high
    lows <- mids-(binwidth/2) # length bins low

    # length and weight at age
    Ages <- seq(start_ages, to = (Amax+1-(1/nseason)), by = (1/nseason)) # all ages including age at 0 if start_ages = 0 and seasonality
    L_a <- Linf*(1-exp(-vbk*(Ages-t0))) # length at age
    W_a <- LWa*L_a^LWb # weight at age

    # fecundity at age
    if(is.null(agefec)){
      Wmat <- LWa*Lmat^LWb
      agefec <- ifelse(W_a<Wmat, 0, W_a-Wmat)
    } 



    ## Selectivity
    S_fa <- rep(NA, length(Ages)) # selectivity by age
      # double-normal selectivity
      # peak - endpoint where selectivity = 1.0
      peak <- S1+binwidth+((0.99*Amax-S1-binwidth)/(1+exp(-S2)))
      # joiner functions for asc and desc components
      jf1 <- (1+exp(-20*(Ages-S1)/(1+abs(Ages-S1))))^-1
      jf2 <- (1+exp(-20*(Ages-peak)/(1+abs(Ages-peak))))^-1
      # t1 and t2
      t1 <- exp(-(start_ages-S1)^2/exp(S3))
      t2 <- exp(-(Amax-peak)^2/exp(S4))
      # asc and desc portions
      asc <- (1+exp(-S5))^-1+(1-(1+exp(-S5))^-1)*(exp(-(Ages-S1)^2/exp(S3))-t1)/(1-t1)
      desc <- 1+(((1+exp(-S6))^-1)-1)*(exp(-(Ages-peak)^2/exp(S4))-1)/(t2-1)
      # selectivity by age
      S_fa <- asc*(1-jf1)+jf1*((1-jf2)+jf2*desc)
      S_fa <- S_fa/max(S_fa)

    # retention
    ret_sigma <- 0.1*sizemin # get sd from CV
    ret <- pnorm((L_a-sizemin)/ret_sigma) # normal distribution of size limit
    ret <- ret/max(ret)



    ###################
    ## output list ####
    ###################
    output <- NULL
      output <- list(nseason = nseason,
                    Amax = length(Ages),
                    Ages = Ages,
                    Linf = Linf,
                    vbk = vbk,
                    t0 = t0,
                    length_bins = data.frame(mids = mids,
                                            highs = highs,
                                            lows = lows),
                    M = M,
                    LWa = LWa,
                    LWb = LWb,
                    agefec = agefec,
                    h = h,
                    cstar = cstar,
                    binwidth = binwidth,
                    start_ages = start_ages,
                    q = q,
                    persis = persis,
                    sig1e = sig1e,
                    ret = ret,
                    sel_param = data.frame(S1 = S1,
                                          S2 = S2,
                                          S3 = S3,
                                          S4 = S4,
                                          S5 = S5,
                                          S6 = S6),
                    eff_param = data.frame(max_eff = max_eff,
                                          w_cost = w_cost),
                    Ut_list = list(cpue = Ut_cpue,
                                  hpue = Ut_hpue,
                                  size = Ut_size,
                                  dist = Ut_dist,
                                  crowd = Ut_crowd),
                    L_a = L_a,
                    W_a = W_a,
                    S_fa = S_fa,
                    mvt_param = data.frame(w_dep = w_dep,
                                          w_hab = w_hab,
                                          w_dist = w_dist),
                    reef_flags = list(AR = data.frame(AR_flag = AR_flag,
                                                    AR_prop_M = AR_prop_M,
                                                    AR_prop_rec = AR_prop_rec,
                                                    AR_prop_U = AR_prop_U,
                                                    AR_prop_q = AR_prop_q),
                                      NR = data.frame(NR_flag = NR_flag,
                                                      NR_prop_M = NR_prop_M,
                                                      NR_prop_rec = NR_prop_rec,
                                                      NR_prop_U = NR_prop_U,
                                                      NR_prop_q = NR_prop_q))
      ) # end of output list

    return(output)
  }) # end of with function
} # end of create_OM_input function