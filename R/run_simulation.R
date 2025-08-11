#' @title Run deterministic malariasimple
#' @description Generates and runs the deterministic malariasimple model
#' @param params List of model parameters generated using get_parameters() and other helper functions
#' @param full_output Boolean variable stating whether model should output model outputs, or a selected few.
#' @examples
#' params <- get_parameters() |>
#'           set_equilibrium(init_EIR = 10)
#' simulation_output <- run_simulation(params)
#' @export

run_simulation <- function(params, full_output = FALSE){
  gen <- malariasimple_deterministic_ITN_IRS_two_species

  sys <- dust2::dust_system_create(gen(), params, n_particles = 1, dt = 1/params$tsd)

  dust2::dust_system_set_state_initial(sys)
  time <- 0:params$n_days
  out <- dust2::dust_system_simulate(sys, time)

  out <- aperm(out,c(2,1))
  out <- out[-1,]

  #Add time to output
  time <- 1:params$n_days
  out <- cbind(time,out)

  #Add colnames to output
  colnames <- get_output_colnames(sys, params)
  dimnames(out)[[2]] <- c("time",colnames)

  #Only output select variables (unless requested otherwise)
  if(!full_output){
    selected_cols <- c("time","mv_species1", "mv_species2",
                       grep("n_any_malaria*", colnames(out), value = TRUE),
                       grep("n_clin_inc*", colnames(out), value = TRUE),
                       grep("^n_pop*", colnames(out), value = TRUE))
    out <- out[,selected_cols]

  }
  return(out)
}

get_output_colnames <- function(sys, params){
  ## Get column names
  colname_ls <- dust2::dust_unpack_index(sys)
  colname_length <- sum(sapply(colname_ls, length))
  index <- vector("character", colname_length)
  pos <- 1

  for (i in seq_along(colname_ls)) {
    len <- length(colname_ls[[i]])
    index[pos:(pos + len - 1)] <- names(colname_ls)[i]
    pos <- pos + len
  }

  #Improve naming of user-defined outputs
  index[index == "n_ud_pop"] <- paste("n_pop",params$age_vector[params$min_age_n], params$population_rendering_max_ages,sep="_")         ## this is population size
  index[index == "n_ud_any_prev"] <- paste("n_any_malaria",params$age_vector[params$min_age_prev], params$prevalence_rendering_max_ages,sep="_") ## this is any malaria
  index[index == "n_ud_detect_prev"] <- paste("n_detectable_malaria",params$age_vector[params$min_age_prev], params$prevalence_rendering_max_ages,sep="_") ## this is detectable malaria
  index[index == "n_ud_inc"] <- paste("n_clin_inc",params$age_vector[params$min_age_inc], params$clin_inc_rendering_max_ages,sep="_")          ## this is clinical incidence
  return(index)
}
