# A two-species multispecies fixture whose first species carries input
# numbers-at-age: fit once at estimateMode = 3, write the fitted N-at-age back as
# NByageFixed, then set the requested estDynamics per species.
fixed_natage_build <- function(d, recFun, msmMode = 1, HCR = NULL) {
  args <- list(data_list = d, inits = NULL, estimateMode = 3, msmMode = msmMode,
               suitMode = 0, initMode = "NonEquilibrium", random_rec = FALSE,
               fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))
  if (!missing(recFun)) args$recFun <- recFun
  if (!is.null(HCR)) args$HCR <- HCR
  suppressMessages(suppressWarnings(do.call(fit_mod, args)))
}

fixed_natage_data <- function(estDynamics) {
  set.seed(123)
  d  <- make_msm_test_data()$data_list
  m0 <- fixed_natage_build(d)
  N  <- m0$quantities$N_at_age
  yrs <- d$styr:(d$styr + dim(N)[4] - 1)
  nb <- do.call(rbind, lapply(seq_along(yrs), function(i)
    data.frame(Species_name = "Species1", Species = 1, Sex = 0, Year = yrs[i],
               t(N[1, 1, , i]))))
  colnames(nb) <- c("Species_name", "Species", "Sex", "Year",
                    paste("Age", seq_len(dim(N)[3])))
  d$NByageFixed <- nb
  d$estDynamics <- estDynamics
  d
}

fixed_dyn_build <- function(estDynamics) {
  fixed_natage_build(fixed_natage_data(estDynamics))
}
