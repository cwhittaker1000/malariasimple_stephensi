#' @title Set equilibrium
#' @description Updates parameter list to include equilibrium values for a given EIR (entomological inoculation rate) for the malariasimple parameter list.
#' This function also includes some 'finishing touches' for the smooth running of the simulation and so is essential that this function is performed last.
#' @param params List of parameters
#' @param init_EIR Value of EIR at equilbrium

#' @examples
#' params <- get_parameters() |>
#'             set_equilibrium(init_EIR = 6)
#' @importFrom stats rlnorm
#' @importFrom stats rbinom
#' @export

set_equilibrium <- function(params, init_EIR)
{
  ##---------------------------------------------------------------------------------------
  ##                                     SOME HOUSEKEEPING
  ##---------------------------------------------------------------------------------------
  #Interventions
  if(params$itn_set == 0){
    params$max_itn_cov <- 0
    params$r_itn_daily <- rep(0,(params$n_days+1))
    params$s_itn_daily <- rep(1,(params$n_days+1))
    params$itn_eff_cov_daily <- rep(0,(params$n_days+1))
  }

  if(params$irs_set == 0){
    params$max_irs_cov <- 0
    params$rs <- rep(0,(params$n_days+1))
    params$ss <- rep(1,(params$n_days+1))
    params$irs_eff_cov_daily <- rep(1,(params$n_days+1))
  }

  #Set intervention group coverages
  if(params$itn_set == 0 & params$irs_set == 0){
    num_int <- 1
    cov <- 1
  } else if(params$itn_set == 1 & params$irs_set == 0){
    num_int <- 2
    cov <- c((1 - params$max_itn_cov), params$max_itn_cov)
  } else if(params$itn_set == 0 & params$irs_set == 1){
    num_int <- 3
    cov <- c((1 - params$max_irs_cov),0,params$max_irs_cov)
  } else if(params$itn_set == 1 & params$irs_set == 1){
    num_int <- 4
    cov <- c((1 - params$max_itn_cov) * (1 - params$max_irs_cov),
             params$max_itn_cov * (1 - params$max_irs_cov),
             (1 - params$max_itn_cov) * params$max_irs_cov,
             params$max_itn_cov * params$max_irs_cov)
  }
  params$cov <- cov
  params$num_int <- num_int
  age_vector <- params$age_vector
  het_brackets <- params$het_brackets

  ## Check Parameters
  if(!is.numeric(init_EIR) | init_EIR < 0) stop("init_EIR must be a numeric value greater than zero")

  ##---------------------------------------------------------------------------------------
  ##                                  POPULATION DEMOGRAPHICS
  ##---------------------------------------------------------------------------------------
  na <- as.integer(length(age_vector))  # number of age groups
  nh <- as.integer(het_brackets)  # number of heterogeneity groups
  h <- statmod::gauss.quad.prob(nh, dist = "normal")

  age_rate <- age_width <- age_mid_point <- den <- c()
  age_vector_plus <- c(age_vector, 100*365) #Assume 100 years is maximum human age
  for (i in 1:(na))
  {
    age_width[i] <- age_vector_plus[i+1] - age_vector_plus[i]
    age_rate[i] <- 1/(age_vector_plus[i + 1] - age_vector_plus[i])  # vector of rates at which people leave each age group (1/age group width)
    age_mid_point[i] <- 0.5 * (age_vector_plus[i] + age_vector_plus[i + 1])  # set age group vector to the midpoint of the group

  }
  age_rate[na] = 0


  den <- 1/(1 + age_rate[1]/params$eta)
  for (i in 1:(na-1))
  {
    den[i+1] <- age_rate[i] * den[i]/(age_rate[i+1] + params$eta)  # proportion in each age_vector group
  }

  ##---------------------------------------------------------------------------------------
  ##                        HETEROGENEITIES IN TRANSMISSION PARAMETERS
  ##---------------------------------------------------------------------------------------
  ## force of infection
  psi <- c()
  for (i in 1:na)
  {
    psi[i] <- 1 - (params$rho * exp(-age_vector[i]/params$a0))  #force of infection for each age group
  }
  fden <- psi * den
  omega <- sum(fden)  #normalising constant

  ## heterogeneity
  het_x <- h$nodes
  het_wt <- h$weights
  den_het <- outer(den, het_wt)
  zeta <- exp(-params$sigma2/2 + sqrt(params$sigma2) * het_x)/sum(het_wt * exp(-params$sigma2/2 + sqrt(params$sigma2) * het_x))

  ## EIR
  EIRY_eq <- init_EIR  # initial annual EIR
  EIRd_eq <- EIRY_eq/params$DY
  EIR_eq <- outer(psi, zeta) * EIRd_eq

  ## Immunity and FOI
  x_I <- den[1]/params$eta
  for (i in 2:na)
  {
    x_I[i] <- den[i]/(den[i - 1] * age_rate[i - 1])  #temporary variables
  }
  fd <- 1 - (1 - params$fD0)/(1 + (age_vector/params$aD)^params$gammaD)

  # maternal immunity begins at a level proportional to the clinical
  # immunity of a 20 year old, this code finds that level
  age20i <- rep(0, na)
  for (i in 2:na)
  {
    age20i[i] <- ifelse(age_vector[i] >= (20 * params$DY) & age_vector[i - 1] < (20 * params$DY),
                        i, age20i[i - 1])
  }
  age20u <- as.integer(age20i[na])
  age20l <- as.integer(age20u - 1)
  age_20_factor <- (20 * params$DY - age_vector[age20l] - 0.5 * age_width[age20l]) *
    2/(age_width[age20l] + age_width[age20u])

  ##---------------------------------------------------------------------------------------
  ##                        INITIAL VALUES FOR IMMUNITY STATES
  ##---------------------------------------------------------------------------------------
  init_IB <- matrix(0, na, nh)
  FOI_eq <- matrix(0, na, nh)
  init_ID <- matrix(0, na, nh)
  init_ICA <- matrix(0, na, nh)
  ICM_init_eq <- vector(length = nh, mode = "numeric")
  init_ICM <- matrix(0, na, nh)
  cA_eq <- matrix(0, na, nh)
  FOIvij_eq <- matrix(0, na, nh)
  p_det_eq <- matrix(0, na, nh)
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      init_IB[i, j] <- (ifelse(i == 1, 0, init_IB[i - 1, j]) +
                          EIR_eq[i,j]/(EIR_eq[i, j] * params$uB + 1) * x_I[i])/(1 + x_I[i]/params$dB)
      FOI_eq[i, j] <- EIR_eq[i, j] * ifelse(init_IB[i, j] == 0, params$b0,
                                            params$b0 * ((1 - params$b1)/(1 + (init_IB[i, j]/params$IB0)^params$kB) + params$b1))
      init_ID[i, j] <- (ifelse(i == 1, 0, init_ID[i - 1, j]) +
                          FOI_eq[i, j]/(FOI_eq[i, j] * params$uD + 1) * x_I[i])/(1 + x_I[i]/params$dID)
      init_ICA[i, j] <- (ifelse(i == 1, 0, init_ICA[i - 1, j]) +
                           FOI_eq[i,j]/(FOI_eq[i, j] * params$uCA + 1) * x_I[i])/(1 + x_I[i]/params$dCA)
      p_det_eq[i, j] <- params$d1 + (1 - params$d1)/(1 + fd[i] * (init_ID[i, j]/params$ID0)^params$kD)
      cA_eq[i, j] <- params$cU + (params$cD - params$cU) * p_det_eq[i, j]^params$gamma1
    }
  }
  # needs to be calculated after because it references ICA
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      ICM_init_eq[j] <- params$PM * (init_ICA[age20l, j] + age_20_factor *
                                       (init_ICA[age20u, j] - init_ICA[age20l, j]))
      init_ICM[i, j] <- ifelse(i == 1,
                               ICM_init_eq[j], init_ICM[i - 1,j])/(1 + x_I[i]/params$dCM)
    }
  }

  IC_eq <- init_ICM + init_ICA
  phi_eq <- params$phi0 * ((1 - params$phi1)/(1 + (IC_eq/params$IC0)^params$kC) + params$phi1)

  ##---------------------------------------------------------------------------------------
  ##                        INITIAL VALUES FOR HUMAN STATES
  ##---------------------------------------------------------------------------------------
  gamma <- params$eta + c(age_rate[1:(na - 1)], 0)
  delta <- c(params$eta, age_rate[1:(na - 1)])

  betaT <- matrix(rep(params$rT + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaD <- matrix(rep(params$rD + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaP <- matrix(rep(params$rP + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)

  aT <- FOI_eq * phi_eq * params$ft/betaT
  aD <- FOI_eq * phi_eq * (1 - params$ft)/betaD
  aP <- params$rT * aT/betaP

  Z_eq <- array(dim = c(na, nh, 4))
  Z_eq[1, , 1] <- den_het[1, ]/(1 + aT[1, ] + aD[1, ] + aP[1, ])
  Z_eq[1, , 2] <- aT[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 3] <- aD[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 4] <- aP[1, ] * Z_eq[1, , 1]

  for (j in 1:nh)
  {
    for (i in 2:na)
    {
      Z_eq[i, j, 1] <- (den_het[i, j] - delta[i] * (Z_eq[i - 1, j, 2]/betaT[i, j] +
                                                      Z_eq[i - 1, j, 3]/betaD[i, j] +
                                                      (params$rT *  Z_eq[i - 1, j, 2]/betaT[i, j]
                                                       + Z_eq[i - 1, j, 4])/betaP[i, j]))/(1 + aT[i, j] + aD[i, j] + aP[i, j])
      Z_eq[i, j, 2] <- aT[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 2]/betaT[i, j]
      Z_eq[i, j, 3] <- aD[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 3]/betaD[i, j]
      Z_eq[i, j, 4] <- aP[i, j] * Z_eq[i, j, 1] + delta[i] * (params$rT *
                                                                Z_eq[i - 1, j, 2]/betaT[i, j] + Z_eq[i - 1, j, 4])/betaP[i,j]

    }
  }

  init_Y <- matrix(Z_eq[, , 1], nrow = na, ncol=nh)
  init_T <- matrix(Z_eq[, , 2], nrow = na, ncol=nh)
  init_D <- matrix(Z_eq[, , 3], nrow = na, ncol=nh)
  init_P <- matrix(Z_eq[, , 4], nrow = na, ncol=nh)

  betaS <- apply(FOI_eq, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaA <- apply(FOI_eq * phi_eq + params$rA, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaU <- apply(FOI_eq + params$rU, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)

  init_A <- matrix(ncol = nh, nrow = na)
  init_U <- matrix(ncol = nh, nrow = na)
  init_S <- matrix(ncol = nh, nrow = na)

  for (i in 1:na)
  {
    for (j in 1:nh)
    {
      init_A[i, j] <- (delta[i] * ifelse(i == 1, 0, init_A[i - 1, j]) +
                         FOI_eq[i, j] * (1 - phi_eq[i, j]) * init_Y[i, j] +
                         params$rD * init_D[i,j])/(betaA[i, j] + FOI_eq[i, j] * (1 - phi_eq[i, j]))
      init_U[i, j] <- (params$rA * init_A[i, j] + delta[i] * ifelse(i == 1,
                                                                    0, init_U[i - 1, j]))/betaU[i, j]
      init_S[i, j] <- init_Y[i, j] - init_A[i, j] - init_U[i, j]
      FOIvij_eq[i, j] <- psi[i] * params$av0 * (params$cT * init_T[i, j] + params$cD *
                                                      init_D[i, j] + cA_eq[i, j] * init_A[i, j] + params$cU * init_U[i, j]) * zeta[j]/omega
    }
  }

  ##---------------------------------------------------------------------------------------
  ##                        INITIAL VALUES FOR MOSQUITO STATES
  ##---------------------------------------------------------------------------------------
  FOIv_eq <- sum(FOIvij_eq)
  init_Iv <- FOIv_eq * params$Surv0/(FOIv_eq + params$mum)
  init_Sv <- params$mum * init_Iv/(FOIv_eq * params$Surv0)
  init_Pv <- 1 - init_Sv - init_Iv

  # mosquito density needed to give this EIR
  mv0 <- omega * EIRd_eq/(init_Iv * params$av0)

  # larval states
  K0 <- 2 * mv0 * params$dLL * params$mum * (1 + params$dPL * params$muPL) * params$gammaL * (params$lambda + 1)/(params$lambda/(params$muLL *
                                                                                                                                   params$dEL) - 1/(params$muLL * params$dLL) - 1)
  init_PL <- 2 * params$dPL * params$mum * mv0
  init_LL <- params$dLL * (params$muPL + 1/params$dPL) * init_PL
  init_EL <- (init_LL/params$dLL + params$muLL* init_LL * (1 + params$gammaL * init_LL/K0))/(1/params$dEL - params$muLL * params$gammaL * init_LL/K0)

  # add in final dimension - interventions
  cov <- params$cov

  mat <- matrix(0, na, nh)

  init_S <- vapply(cov, FUN = function(x)
  {
    x * init_S
  }, mat)
  init_T <- vapply(cov, FUN = function(x)
  {
    x * init_T
  }, mat)
  init_D <- vapply(cov, FUN = function(x)
  {
    x * init_D
  }, mat)
  init_A <- vapply(cov, FUN = function(x)
  {
    x * init_A
  }, mat)
  init_U <- vapply(cov, FUN = function(x)
  {
    x * init_U
  }, mat)
  init_P <- vapply(cov, FUN = function(x)
  {
    x * init_P
  }, mat)

  init_IB = array(init_IB, c(na, nh, num_int))
  init_ID = array(init_ID, c(na, nh, num_int))
  init_ICA = array(init_ICA, c(na, nh, num_int))
  init_ICM = array(init_ICM, c(na, nh, num_int))

  if(!is.null(params$ncc)){
    init_IB = array(init_IB, c(na, nh, num_int, params$ncc))
    init_ID = array(init_ID, c(na, nh, num_int, params$ncc))
    init_ICA = array(init_ICA, c(na, nh, num_int, params$ncc))
    init_ICM = array(init_ICM, c(na, nh, num_int, params$ncc))

    # add in final dimension - interventions
    all_rounds = params$MDA_grp_prop*params$MDA_cov
    ccov = c(all_rounds, 1-all_rounds)

    mat2 <- array(0, c(na,nh, num_int))
    init_S <- vapply(ccov,FUN = function(x){x * init_S},mat2)
    init_T <- vapply(ccov,FUN = function(x){x * init_T},mat2)
    init_D <- vapply(ccov,FUN = function(x){x * init_D},mat2)
    init_A <- vapply(ccov,FUN = function(x){x * init_A},mat2)
    init_U <- vapply(ccov,FUN = function(x){x * init_U},mat2)
    init_P <- vapply(ccov,FUN = function(x){x * init_P},mat2)
  }

  # better het bounds for equilbirum initialisation in individual model
  zetas <- rlnorm(n = 1e5,meanlog = -params$sigma2/2, sdlog = sqrt(params$sigma2))
  while(sum(zetas>100)>0){
    zetas[zetas>100] <- rlnorm(n = sum(zetas>100),meanlog = -params$sigma2/2, sdlog = sqrt(params$sigma2))
  }

  wt_cuts <- round(cumsum(het_wt)*1e5)
  zeros <- which(wt_cuts==0)
  wt_cuts[zeros] <- 1:length(zeros)
  larges <- which(wt_cuts==1e5)
  wt_cuts[larges] <- (1e5 - (length(larges)-1)):1e5
  wt_cuts <- c(0,wt_cuts)
  het_bounds <- sort(zetas)[wt_cuts]
  het_bounds[length(het_bounds)] <- (params$max_age/365)+1

  ##---------------------------------------------------------------------------------------
  ##                        UNUSED DEFAULT PARAMETERS
  ##---------------------------------------------------------------------------------------


  params$equilibrium_set <- 1

  #Climate
  if (is.null(params$daily_rain_input))
    params$daily_rain_input <- rep(1, (params$n_days+1))
  if (is.null(params$daily_temp))
    params$daily_temp <- rep(1, (params$n_days+1))

  human_pop <- params$human_pop
  if(params$stochastic){
    human_init_list <- list(init_S = init_S, init_T = init_T,
                            init_D = init_D, init_A = init_A,
                            init_U = init_U, init_P = init_P,
                            init_Y = init_Y)
    init_count_list <- lapply(human_init_list, function(arr) {
      array(rbinom(length(arr), size = human_pop, prob = arr), dim = dim(arr))
    })
    params$init_S <- init_count_list$init_S
    params$init_T <- init_count_list$init_T
    params$init_D <- init_count_list$init_D
    params$init_A <- init_count_list$init_A
    params$init_U <- init_count_list$init_U
    params$init_P <- init_count_list$init_P
    params$init_Y <- init_count_list$init_Y
  } else {
    params$init_S <- init_S*human_pop
    params$init_T <- init_T*human_pop
    params$init_D <- init_D*human_pop
    params$init_A <- init_A*human_pop
    params$init_U <- init_U*human_pop
    params$init_P <- init_P*human_pop
    params$init_Y <- init_Y*human_pop
  }

  params$init_IB <- init_IB
  params$init_ID <- init_ID
  params$init_ICA <- init_ICA
  params$init_ICM <- init_ICM

  params$ICM_init_eq = ICM_init_eq

  params$days <- 0:(params$n_days)
  params$iterations <- 0:(params$n_ts)

  ## collate init
  res <- list( init_Iv = init_Iv, init_Sv = init_Sv,
               init_Pv = init_Pv, init_PL = init_PL, init_LL = init_LL, init_EL = init_EL,
               age_rate = age_rate, het_wt = het_wt, het_x = het_x,
               omega = omega, psi = psi, zeta = zeta, K0 = K0, 
               mv0 = mv0, 
               x_I = x_I,
               FOI_eq = FOI_eq, EIR_eq = EIR_eq, cA_eq = cA_eq,
               den = den, FOIv_eq = FOIv_eq,
               betaS = betaS, betaA = betaA, betaU = betaU, FOIvij_eq=FOIvij_eq,
               age_mid_point = age_mid_point, het_bounds = het_bounds, pi = pi,
               age20l = age20l, age20u = age20u, age_20_factor = age_20_factor)

  res <- append(res,params)

  return(res)
}

