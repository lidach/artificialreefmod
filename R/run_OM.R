##########################
## Calculate distance ####
##########################

## PC dist - caculate distance to each site
earth.dist <- function (long1, lat1, long2, lat2)
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.137   #Mean Earth Radius
  d <- R * c
  return(d)
}



#######################################
## Artificial reef operating model ####
#######################################
#' @title run_OM
#'
#' @description function to create initial operating model
#'
#' @param seed for stochasticity
#' @param nyears_forward run model for how many years (excluding burn in years)
#' @param nyears_init how many years to run burn in period
#' @param AR_year
#' @param AR_site_nums
#' @param coord_file data frame with red snapper habitat suitability, NR, depth, longitude, and latitude coordinates for each cell/site
#' @param landing_site_file data frame with longitude and latitude coordinates for landing sites
#' @param RS_tag_file tagging data used for movement parameterization
#' @param RS_data_file numbers at age 0-10 for red snapper, used for movement parameterization
#' @param H_type default = NULL; habitat scenario ("null", "deg", "real")
#' @param lh_list data list for creating operating model, includes data from create_OM_input()

run_OM <- 
function(nyears_forward = 50,
        nyears_init = 50,
        AR_year = 30,
        AR_site_nums = NULL,
        coord_file,
        landing_site_file,
        RS_tag_file,
        RS_data_file,
        H_type = "null",
        lh_list
        )
{
  with(lh_list, {

  #####################
  ## Control Panel ####
  #####################
  ## Temporal parameters
  nyears <- nyears_forward+nyears_init
  fish_strt <- (nyears_init+1):nyears
  AR_strt <- (nyears_init+AR_year+1):nyears

  ## Spatial inputs
  nsites <- nrow(coord_file) # number of cells/fishing sites
  nlandings <- nrow(landing_site_file) # number of landing sites
  long <- coord_file$long
  lat <- coord_file$lat
  depth <- coord_file$depth

  # distance for landing sites
  dist <- matrix(NA, nrow = nsites, ncol = nrow(landing_site_file))
  for(j in 1:ncol(dist)){
    for(i in 1:nsites){
      xone = long[i]
      yone = lat[i]
      # dist[i,j] <- sqrt((landing_site_file$long[j]-xone)^2 + (landing_site_file$lat[j]-yone)^2)
      dist[i,j] <- earth.dist(landing_site_file$long[j], landing_site_file$lat[j], xone, yone)
    }
  }
  # weighted mean of landing sites to each site
  dist.wt <- apply(dist, 1, function(x) weighted.mean(x, landing_site_file$pop))

  # AR and NR sites
  NR <- coord_file$HB
  AR <- rep(0, nsites)
    AR_sites <- AR_site_nums
    AR[AR_sites] <- 1
  NR_sites <- which(NR>0)

  ## H scenarios
  # depth penalty for R0, in meters
  dep_loc <- which(depth < 15 | depth > 75)
  # near shore penalty
  city_pen <- rep(0, nsites)
    city_pen[which(dist.wt < 100)] <- 1-(dist.wt[which(dist.wt < 100)]/max(dist.wt[which(dist.wt < 100)]))
  # different habitat quaility types and penalty for near shore
  if("null" %in% H_type)	H_param <- rep(1, nsites)
  if("deg" %in% H_type)	H_param <- rep(0.75, nsites)
  if("real" %in% H_type){
    H_param <- rep(1,nsites)
    H_param <- 1-city_pen
  }

  ## R0
  Peast <- 0.23 # Proportion of recruits allocated to east of the Mississippi (per RS assessment)
  R0 <- 1.63e8 # unfished recruitment (SEDAR 52)
  FLsubset <- RS_data_file[RS_data_file$Lon > -85.5 & RS_data_file$Lon < -84.7 & RS_data_file$Lat > 28.8 & RS_data_file$Lat < 29.6 & RS_data_file$Depth < 70, ]
  EastMS <- RS_data_file[RS_data_file$Lon > -89 & RS_data_file$Depth < 70, ]
  PropR_FL <- dim(FLsubset)[1]/dim(EastMS)[1]
  R0_FL <- R0*Peast*PropR_FL
  R0_site <- rep(R0_FL/nsites, nsites)
    R0_site <- R0_site*(1-city_pen)
  # for local recruitment
  H <- matrix(rep(H_param, each = nyears), nrow = nyears, ncol = nsites, byrow = FALSE)

  ## cost - for effort allocation
  cost_null <- 50 + dist.wt*0 # this represents low differences in costs between sites so HIGH movement of anglers (all sites have same cost)
  # cost_mult <- 50 + dist.wt*3 # intermediate # not yet exported
  # cost_exp <- 5 + dist.wt^1.5 # this represents great differences in costs between sites, so LOW movement # not yet exported
  cost <- cost_null


  ## Influence of ARs and NRs
  # natural mortality
  M_site <- array(NA, c(nyears, Amax, nsites))
    for(i in 1:nyears) for(j in 1:nsites) M_site[i,,j] <- M
  # catchability
  q_t <- matrix(0, nrow = nyears, ncol = nsites)
    q_t[fish_strt,] <- q
  U_AR <- matrix(1, nrow = nyears, ncol = nsites)
  U_NR <- matrix(1, nrow = nyears, ncol = nsites)

  if(reef_flags$AR$AR_flag == 1){
    M_site[AR_strt,,AR_sites] <- M_site[AR_strt,,AR_sites]*(1+reef_flags$AR$AR_prop_M)
    H[AR_strt,AR_sites] <- H[AR_strt,AR_sites]*(1+reef_flags$AR$AR_prop_rec)
    q_t[AR_strt,AR_sites] <- q_t[AR_strt,AR_sites]*(1+reef_flags$AR$AR_prop_q)
    U_AR[AR_strt,AR_sites] <- U_AR[AR_strt,AR_sites]*(1+reef_flags$AR$AR_prop_U)
  }

  if(reef_flags$NR$NR_flag == 1){
    M_site[,,NR_sites] <- M_site[,,NR_sites]*(1+reef_flags$NR$NR_prop_M)
    H[,NR_sites] <- H[,NR_sites]*(1+reef_flags$NR$NR_prop_rec)
    q_t[,NR_sites] <- q_t[,NR_sites]*(1+reef_flags$NR$NR_prop_q)
    U_NR[,NR_sites] <- U_NR[,NR_sites]*(1+reef_flags$NR$NR_prop_U)
  }



  #############################
  ## Recruitment processes ####
  #############################
  ## Survivorship
  So <- c(1,cumprod(exp(-M))[1:Amax-1]) 

  ## Habitat-based Beverton-Holt
  EPR0 <- sum(agefec*So) # unfished eggs biomass-per-recruit
  E0 <- sum(R0_site*EPR0)
  # BH recruitment
  a_rec <- EPR0*((1-h)/(4*h))
  b_rec <- ifelse(R0_site > 0, (5*h-1)/(4*h*R0_site), 0)
  # BH recruitment by site (H matrix) - habitat specific
  ahab <- a_rec*H
  bhab <- t(t(H)*b_rec)/H



  ################
  ## Movement ####
  ################
  ## movement matrix - dist_mat
  dist_mat <- matrix(NA, nrow = nsites, ncol = nsites)
  for(i in 1:nsites){
    for(j in 1:nsites){
      dist_mat[i,j] <- earth.dist(long[i], lat[i], long[j], lat[j])
    }
  }

  ## Depth
  # aggregate, mean numbers caught at depth and sd
    agg_list <- list()
    for(i in 1:11){
      agg_list[[i]] <- do.call(data.frame, aggregate(RS_data_file[,50+i], by = list(RS_data_file$Depth), FUN = function(x) c(Av = mean(x), sig = sd(x))))
    }
  # fitting mean depth at age and variance using VBF, then assuming those follow a normal
    means <- sapply(agg_list, FUN = function(x) weighted.mean(x[,1], w = x[,2]))
     mod_means <- nls(means ~ C*(1-exp(-r*(0:10-a))), start = list(C = 80, r = 1, a = 1))
      pred_means <- predict(mod_means)[1:11]
    vars <- sapply(agg_list, FUN = function(x) wtd.var(x[,1], w = x[,2]))
      mod_vars <- nls(vars ~ C*(1-exp(-r*(0:10-a))), start = list(C = 1000, r = 1, a = 1))
      pred_vars <- predict(mod_vars)[1:11]
  # calculations
    pref_depth <- matrix(NA, ncol = Amax, nrow = ceiling(max(depth)))
    pref_depth_cumd <- matrix(NA, ncol = Amax, nrow = ceiling(max(depth)))
    for(i in 1:length(pred_means)){
      pref_depth[,i] <- dnorm(seq(1, ceiling(max(depth))), mean = pred_means[i], sd = sqrt(pred_vars[i]))
      for(j in 1:length(pred_means)){
        pref_depth_cumd[j,i] <- pnorm(j+0.5, mean = pred_means[i], sd = sqrt(pred_vars[i])) - pnorm(j-0.5, mean = pred_means[i], sd = sqrt(pred_vars[i]))
      }
    }
    pref_depth[,12:Amax] <- pref_depth[,11]
  # every column has max at 1
    pref_depth_std <- t(t(pref_depth)/apply(pref_depth, 2, max))
    row.names(pref_depth_std) <- seq(1, ceiling(max(depth)))

  ## Substrate
  coord_dat <- data.frame(AR = AR/sum(AR),
                          HB = coord_file$HB/sum(coord_file$HB),
                          RS_HS = coord_file$RS_HS/sum(coord_file$RS_HS))
    coord_dat$AR[is.na(coord_dat$AR)] <- 0
  sub_mod <-  lm(RS_HS ~ AR+HB, data = coord_dat)
  hab_suit_per <- predict(sub_mod)

  ## Distance
  RS_tag <- RS_tag_file[!is.na(RS_tag_file[,1]),]	
  RS_tag$travdist <- earth.dist(RS_tag[,"Lon1"], RS_tag[,"Lat1"], RS_tag[,"Lon2"], RS_tag[,"Lat2"]) # Calculating the distance a red snapper traveled
  RS_tag$travdist <- RS_tag$travdist*(365/RS_tag$days_at_large) # Calculating the distance they would have traveled in a full year based on their days at large
  move_RS <- function(theta){
    lam_dist <- exp(theta[1])
    NLL <- -1*sum(dexp(RS_tag$travdist[RS_tag$days_at_large > 200], lam_dist, log = TRUE))
    return(NLL)
  }
  move_distcost <- optim(log(0.05), fn = move_RS, method = "BFGS")
  lam <- exp(move_distcost$par)


  ## movement matrix without threshold preference (initial year)
  pmove_init <- array(0, dim = c(nsites, nsites, Amax))
  for(a in 1:Amax){
    pmove_init[,,a] <- t(t(exp(-lam*dist_mat)^mvt_param$w_dist) * (pref_depth_std[round(depth), a]^mvt_param$w_dep) * (hab_suit_per^mvt_param$w_hab))
    pmove_init[,,a] <- pmove_init[,,a]/rowSums(pmove_init[,,a])
  }



  ################################################
  ## Creating empty vectors, matrices, arrays ####
  ################################################
  Nage <- array(0, dim = c(nyears, Amax, nsites)) # Actual population-at-age in each cell for each year before mortality
  Recruit <- matrix(0, nrow = nyears, ncol = nsites) # recruitment
  Eggs <- matrix(0, nrow = nyears, ncol = nsites) # eggs
  Larv <- matrix(0, nrow = nyears, ncol = nsites) # larvae
  VB <- matrix(0, nrow = nyears, ncol = nsites) # vulnerable/exploitable biomass (for now used in gravity model)
  Catch_bio <- matrix(0, nrow = nyears, ncol = nsites) # total catch (biomass)
  Catch_num <- matrix(0, nrow = nyears, ncol = nsites) # catch at age
  Catch_numage <- array(0, dim = c(nyears, Amax, nsites)) # total catch
  FM <- array(0, dim = c(nyears, Amax, nsites)) # fishing mortality
  pmove <- array(0, dim = c(nsites, nsites, Amax)) # movement
  Eff_space <- matrix(0, nrow = nyears, ncol = nsites) # effort expended
  Sel <- matrix(S_fa, nrow = Amax, ncol = nsites) # selectivity by space

  # Utiltiy metrics for effort dynamic
  CPUE_Ut <- HPUE_Ut <- Size_Ut <- matrix(NA, nrow = nyears, ncol = nsites)
  U_cpue <- U_hpue <- U_size <- U_dist <- U_crd <- Tot_Ut <- matrix(NA, nrow = nyears, ncol = nsites) 
  Tot_Ut_sum <- NULL
  pmax_eff <- persis_pmax_eff <- NULL
  grav_wt <- matrix(NA, nrow = nyears, ncol = nsites)



  ############################
  ## Model Initialization ####
  ############################
  # effort dynamic
  tot_eff <- eff_param$max_eff
  pmax_eff[1] <- sum(Eff_space[1,])/tot_eff
  persis_pmax_eff[1] <- pmax_eff[1]

  # numbers  at age
  Recruit[1,] <- Nage[1,1,] <- R0_FL*((pref_depth_std[round(depth),1]*hab_suit_per)/sum(pref_depth_std[round(depth),1]*hab_suit_per))
  Nage[1,,] <- R0_FL*((pref_depth_std[round(depth),1]*hab_suit_per)/sum(pref_depth_std[round(depth),1]*hab_suit_per))

  # recruitment
  Eggs[1,] <- colSums(agefec*Nage[1,,])
  Larv[1,] <- Eggs[1,]*((pref_depth_std[round(depth),1]*hab_suit_per)/sum(pref_depth_std[round(depth),1]*hab_suit_per)) 

  # initial movement
  for(a in 2:Amax) Nage[1,a,] <- Nage[1,a,] %*% pmove_init[,,a]
  VB[1,] <- colSums(S_fa*Nage[1,,]*W_a)

  # utility metrics
  CPUE_Ut[1,] <- colSums(Nage[1,,]*S_fa)
  HPUE_Ut[1,] <- colSums(Nage[1,,]*S_fa*ret)
  Size_Ut[1,] <- colSums(Nage[1,,]*S_fa*L_a)/colSums(Nage[1,,]*S_fa)
  U_cpue[1,] <- Ut_list$cpue[1]/(1+exp(Ut_list$cpue[2]*(Ut_list$cpue[3]-CPUE_Ut[1,])))
  U_hpue[1,] <- Ut_list$hpue[1]/(1+exp(Ut_list$hpue[2]*(Ut_list$hpue[3]-HPUE_Ut[1,])))
  U_size[1,] <- Ut_list$size[1]/(1+exp(Ut_list$size[2]*(Ut_list$size[3]-Size_Ut[1,])))
  U_dist[1,] <- Ut_list$dist[1]/(1+exp(Ut_list$dist[2]*(Ut_list$dist[3]-dist.wt)))
  U_crd[1,] <- Ut_list$crowd[1]/(1+exp(Ut_list$crowd[2]*(Ut_list$crowd[3]-Eff_space[1,])))
  Tot_Ut[1,] <-  U_cpue[1,]+U_hpue[1,]+U_size[1,]+U_dist[1,]+U_crd[1,]
  Tot_Ut_sum[1] <- sum(Tot_Ut[1,])
  Tot_init <- 1*Tot_Ut_sum[1]



  #####################
  ## Time Dynamics ####
  #####################		
  for(y in 2:nyears){
    ## Effort dynamics
    for(k in 1:nsites) grav_wt[y,k] <- (Tot_Ut[y-1,k]/cost[k])^eff_param$w_cost
    pmax_eff[y] <- 1/(1+exp(-(Tot_Ut_sum[y-1]-Tot_init)/(sig1e*Tot_init)))
    persis_pmax_eff[y] <- pmax_eff[y]*(1-persis)+pmax_eff[y-1]*persis # accounts for some "stickyness" in the effort time step to time step
    if(y > fish_strt[1]) for(k in 1:nsites) Eff_space[y,k] <- (tot_eff*persis_pmax_eff[y]*grav_wt[y,k])/sum(grav_wt[y,])

    ## Fishing mortality
    FM[y,,] <- q_t[y,]*t(t(Sel)*Eff_space[y,])

    ## Recruitment
    Larv[y,] <- Eggs[y-1,]*(pref_depth_std[round(depth),1]*hab_suit_per)/sum(pref_depth_std[round(depth),1]*hab_suit_per)
    Recruit[y,] <- Nage[y,1,] <- ahab[y,]*Larv[y,]/(1+bhab[y,]*Larv[y,])

    ## Numbers at age
    Nage[y,2:Amax,] <- Nage[y-1,1:(Amax-1),]*exp(-(M_site[y-1,1:(Amax-1),]+FM[y-1,1:(Amax-1),]))
    # plus group
    Nage[y,Amax,] <- Nage[y,Amax,] + Nage[y-1,Amax,]*exp(-(M_site[y-1,Amax,]+FM[y-1,Amax,]))

    ## metric calcs
    VB[y,] <- colSums(Sel*Nage[y,,]*W_a)
    Catch_bio[y,] <- colSums((FM[y,,]/(FM[y,,]+M_site[y,,]))*Nage[y,,]*(1-exp(-(M_site[y,,]+FM[y,,])))*W_a)
    Catch_num[y,] <- colSums((FM[y,,]/(FM[y,,]+M_site[y,,]))*Nage[y,,]*(1-exp(-(M_site[y,,]+FM[y,,]))))
    Catch_numage[y,,] <- (FM[y,,]/(FM[y,,]+M_site[y,,]))*Nage[y,,]*(1-exp(-(M_site[y,,]+FM[y,,])))


    ## Utilty calcs
    CPUE_Ut[y,] <- Catch_num[y,]/Eff_space[y,]
      CPUE_Ut[is.na(CPUE_Ut)] <- 0
    HPUE_Ut[y,] <- colSums(Catch_numage[y,,]*ret)/Eff_space[y,]
      HPUE_Ut[is.na(HPUE_Ut)] <- 0
    Size_Ut[y,] <- colSums(Nage[y,,]*L_a*S_fa)/colSums(Nage[y,,]*S_fa)
      Size_Ut[is.na(Size_Ut)] <- 0
    U_cpue[y,] <- Ut_list$cpue[1]/(1+exp(Ut_list$cpue[2]*(Ut_list$cpue[3]-CPUE_Ut[y,])))
    U_hpue[y,] <- Ut_list$hpue[1]/(1+exp(Ut_list$hpue[2]*(Ut_list$hpue[3]-HPUE_Ut[y,])))
    U_size[y,] <- Ut_list$size[1]/(1+exp(Ut_list$size[2]*(Ut_list$size[3]-Size_Ut[y,])))
    U_dist[y,] <- Ut_list$dist[1]/(1+exp(Ut_list$dist[2]*(Ut_list$dist[3]-dist.wt))) 
    U_crd[y,] <- Ut_list$crowd[1]/(1+exp(Ut_list$crowd[2]*(Ut_list$crowd[3]-Eff_space[y-1,])))
    Tot_Ut[y,] <- U_cpue[y,]+U_hpue[y,]+U_size[y,]+U_dist[y,]+U_crd[y,]+U_AR[y,]+U_NR[y,]
    Tot_Ut_sum[y] <- sum(Tot_Ut[y,])

    ## Eggs calculation
    Eggs[y,] <- colSums(agefec*Nage[y,,])

    ## Movement
    ind_vec <- ifelse(colSums(t(t(Nage[y,,])*L_a^2)) > quantile(colSums(t(t(Nage[1,,])*L_a^2)), 0.75), 1, 0)
    for(a in 1:Amax){
      pmove[,,a] <- t(t(exp(-lam*dist_mat)^mvt_param$w_dist) * (pref_depth_std[round(depth),a]^mvt_param$w_dep) * (hab_suit_per^mvt_param$w_hab) * ((1/colSums(t(t(Nage[y,,]*L_a^2))/quantile(colSums(t(t(Nage[y,,])*L_a^2)), 0.75))^0.5)^ind_vec))
      pmove[,,a] <- pmove[,,a]/rowSums(pmove[,,a])
    }
    for(a in 2:Amax) Nage[y,a,] <- Nage[y,a,] %*% pmove[,,a]
  } # end of y loop



  #################################
  ## Out of loop calculations ####
  #################################
  model_years <- (nyears_init+1):(nyears_forward+nyears_init)
  FM <- FM[model_years,,]
  Depletion <- rowSums(Eggs[model_years,])/sum(E0)
  Eggs <- Eggs[model_years,]
  VB <- VB[model_years,]
  Nage <- Nage[model_years,,]
  Recruit <- Recruit[model_years,]
  Catch_bio <- Catch_bio[model_years,]
  Catch_num <- Catch_num[model_years,]
  Catch_numage <- Catch_numage[model_years,,]
  Eff_space <- Eff_space[model_years,]
  Tot_Ut <- Tot_Ut[model_years,]
  CPUE_Ut <- CPUE_Ut[model_years,]
  HPUE_Ut <- HPUE_Ut[model_years,]
  Size_Ut <- Size_Ut[model_years,]



  ###################
  ## output list ####
  ###################
  output <- NULL
    output<- list(dist = dist.wt,
                  AR_sites = AR_sites,
                  Eff_space = Eff_space,
                  FM = FM,
                  Tot_Ut = Tot_Ut,
                  CPUE_Ut = CPUE_Ut,
                  HPUE_Ut = HPUE_Ut,
                  Size_Ut = Size_Ut,
                  Tot_init = Tot_init,
                  Catch_num = Catch_num,
                  Catch_bio = Catch_bio,
                  Catch_numage = Catch_numage,
                  Depletion = Depletion,
                  Eggs = Eggs,
                  VB = VB,
                  Nage = Nage,
                  Recruit = Recruit
    ) # end of output list

  return(output)
  }) # end with function
} # end of run_OM function


