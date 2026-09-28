# Comparison of Rceattle and WHAM fit to conditional age-at-length (CAAL) data
# Modified from https://giancarlomcorrea.netlify.app/labs/OFI_WK_2023/examples/case3.R
# Note that the plus group in year 1 is wrong on WHAM line 1720 of the cpp and needs adjustment

# OSA residuals: Rceattle vs WHAM --------------------------------------------
# Cross-check Rceattle's one-step-ahead (OSA) residuals against WHAM's (the
# reference implementation that Rceattle's composition decomposition was ported
# from). Both packages are fit to the same GOAcod + whamGrowthData; with the
# fits matched (above), the aggregate index/catch OSA residuals (lognormal) and
# the composition residuals should fall on the 1:1 line. Conditioning order can
# differ slightly between packages, so expect a tight correlation rather than an
# exact equality.
#
# NOTE: manual comparison using the WHAM growth fork installed at the top of this
# script (GiancarloMCorrea/wham, ref = 'growth'). Its OSA residuals come from
# make_osa_residuals(), and its composition decomposition (src/age_comp_osa.hpp)
# is the same conditional binomial / beta-binomial that Rceattle ported into
# comp_osa.hpp -- so the residuals should agree (up to conditioning order).

# remotes::install_github(repo = 'GiancarloMCorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))
library(readxl)
library(ggplot2)
library(dplyr)
source("tests/comparison/WHAM_growth_comparison/helper.R")
dir.create("tests/comparison/WHAM_growth_comparison")
runmodels = FALSE

# Data ------------------------------------------------------------------
catch_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 1)
index_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 2)
fsh_lcomp_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 3)
index_lcomp_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 4)
caal_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 5)
maturity_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 6)


# Rceattle -------------------------------------------------------------
library(Rceattle)
data("GOAcod")
simData <- GOAcod


# * Make input data ----
# * Data controls
simData$spnames <- "WHAM data"
simData$nspp <- 1
years = 1976:2020 # years
nyrs <- length(years)
simData$styr <- 1976
simData$endyr <- 2020
simData$projyr <- 2030
simData$nsex <- 1
simData$nages <- 10
simData$minage <- 1
lengths = seq(from = 2, to = 130, by = 2) # lengths
simData$nlengths <- length(lengths)
simData$pop_wt_index <- 1
simData$ssb_wt_index <- 1
simData$spawn_month <- 0
simData$alpha_wt_len <- 5.56e-06
simData$beta_wt_len <- 3.2
simData$pop_age_transition_index <- 1
simData$sigma_rec <- 0.6
nages <- simData$nages
nlengths <- simData$nlengths

# * Fleet control
simData$fleet_control <- simData$fleet_control[c(1,3),] # BT and Trawl Fishery are both simple logistic
simData$fleet_control$Fleet_name <- c("Survey", "Fishery")
simData$fleet_control$Fleet_code <- 1:2
simData$fleet_control$Selectivity <- "Logistic"                  # Length-based logistic
simData$fleet_control$Selectivity_dimension <- "Length"          # Length-based logistic
simData$fleet_control$Selectivity_index <- 1:2
simData$fleet_control$Weight_index <- 1
simData$fleet_control$Observation_units <- 1
simData$fleet_control$Month <- c(0, 6)
simData$fleet_control$Bin_first_selected <- 1
simData$fleet_control$Catchability <- c(0, NA)
simData$fleet_control$Catchability_init <- c(1, NA)

# * Index data
simData$index_data <- data.frame(Fleet_name = "Survey",
                                 Fleet_code = 1,
                                 Species = 1,
                                 Year = years,
                                 Month = index_df$fr_yr,
                                 Selectivity_block = 1,
                                 Q_block = 1,
                                 Observation = index_df$index,
                                 Log_sd = sqrt(log(index_df$cv^2 + 1))) # WHAM converts CV -> sigma

# * Catch data
simData$catch_data <- data.frame(Fleet_name = "Fishery",
                                 Fleet_code = 2,
                                 Species = 1,
                                 Year = years,
                                 Month = 0,
                                 Selectivity_block = 1,
                                 Catch = catch_df$catch,
                                 Log_sd = sqrt(log(catch_df$cv^2 + 1))) # WHAM converts CV -> sigma

# * Comp data
# - Index length comp
tmp <- as.matrix(index_lcomp_df[,3:67])
colnames(tmp) <- paste0("Comp_",1:simData$nlengths)
index_comp <- cbind(data.frame(Fleet_name = "Survey",
                               Fleet_code = 1,
                               Species = 1,
                               Sex = 0,
                               Age0_Length1 = 1,
                               Year = years, # Negative turns data off
                               Month = 0,
                               Sample_size = index_lcomp_df$Nsamp),
                    tmp
)

# - Fishery length comp
tmp <- as.matrix(fsh_lcomp_df[,3:67])
colnames(tmp) <- paste0("Comp_",1:simData$nlengths)
fishery_comp <- cbind(data.frame(Fleet_name = "Fishery",
                                 Fleet_code = 2,
                                 Species = 1,
                                 Sex = 0,
                                 Age0_Length1 = 1,
                                 Year = years, # Negative turns data off
                                 Month = 0,
                                 Sample_size = fsh_lcomp_df$Nsamp),
                      tmp
)

simData$comp_data <- rbind(index_comp, fishery_comp)


# * CAAL data
# - Index
index_caal <- caal_df
colnames(index_caal) <- c("Year", "Length", "Sample_size", paste0("CAAL_",1:simData$nages))
simData$caal_data <- index_caal |>
  dplyr::mutate(Fleet_name = "Survey",
                Fleet_code = 1,
                Species = 1,
                Sex = 0) |>
  select(Fleet_name, Fleet_code, Species, Sex, Year, Length, Sample_size, paste0("CAAL_",1:simData$nages)) |>
  as.data.frame()
# - No fishery CAAL

# * Empirical selectivity
simData$emp_sel[] <- NA

# * Fixed numbers
simData$NByageFixed[] <- NA


# * Age transition matrix
tmp <- as.data.frame(diag(1,simData$nlengths))[1:simData$nages,]
colnames(tmp) <- paste0("Length_",1:simData$nlengths)
simData$age_trans_matrix <- cbind(data.frame(Age_transition_name = "Base",
                                             Age_transition_index = 1,
                                             Species = 1,
                                             Sex = 0,
                                             Age = 1:simData$nages),
                                  tmp
)


# * Age error
tmp <- as.data.frame(diag(1,simData$nages))
colnames(tmp) <- paste0("Obs_age",1:simData$nages)
simData$age_error <- cbind(data.frame(Species = 1,
                                      True_age = 1:simData$nages),
                           tmp
)

# * Empirical weight-at-age
WAA <- as.data.frame(matrix(1:simData$nages, ncol = simData$nages))
colnames(WAA) <- paste0("Age",1:simData$nages)
simData$weight <- cbind(data.frame(Wt_name = "Base",
                                   Wt_index = 1,
                                   Species = 1,
                                   Sex = 0,
                                   Year = 0), # 0 fills all years
                        WAA
)

# * Maturity (time-invariant in Rceattle; the WHAM input repeats one ogive over
# all 45 years, so row 1 is the whole schedule). as.matrix() first: a tibble row
# is a list, and matrix() would build a list-matrix that reads back as NA.
MatAA <- matrix(as.matrix(maturity_df[1,2:11]), ncol = simData$nages)
colnames(MatAA) <- paste0("Age",1:simData$nages)
simData$maturity <- cbind(data.frame(Species = 1),
                          MatAA
) |>
  as.data.frame()

# * Sex ratio
sexratio <- as.data.frame(matrix(0.5, ncol = simData$nages))
colnames(sexratio) <- paste0("Age",1:simData$nages)
simData$sex_ratio <- cbind(data.frame(Species = 1),
                           sexratio
)

# * Input mortality-at-age
mort <- as.data.frame(matrix(0.35, ncol = simData$nages))
colnames(mort) <- paste0("Age",1:simData$nages)
simData$M1_base <- cbind(data.frame(Species = 1,
                                    Sex = 0),
                         mort
)

# * Environmental data
simData$env_data <- data.frame(Year = years,
                               EnvData = 1)


# * Ration data
pyrs <- as.data.frame(matrix(1, nrow = 1, ncol = simData$nages))
colnames(pyrs) <- paste0("Age",1:simData$nages)
simData$ration_data <- cbind(data.frame(Species = 1,
                                        Sex = 0,
                                        Year = 0), # Zero fills out all years
                             pyrs
)


# Growth fn matched to WHAM: fix the growth SD (WHAM fixes SDgrowth_par) via the
# sd_L1 / sd_Linf linkage endpoints with est_phase = 0 (fix at init).  Init values
# are log-SD endpoints, threaded onto growth_log_sd (= WHAM SDgrowth_par = log(c(1,9))).
growthFun_vb <- build_growth(fun = "vonBertalanffy",
                             linkages = list(
                               sd_L1   = linkage_spec(~1, init = list(`(Intercept)` = 1), est_phase = 0L),
                               sd_Linf = linkage_spec(~1, init = list(`(Intercept)` = 9), est_phase = 0L)))

# * Estimate Rceattle ----
# Rceattle OSA residuals: estimateMode = 1 (objective differentiable, hindcast
# only). osa_residuals() builds the composition OSA data on demand.
# Cache the Rceattle OSA (fit + one-step-ahead residuals) to a gitignored RDS;
# delete it (or set runmodels = TRUE) to recompute after the model changes.
rce_osa_cache <- "tests/comparison/WHAM_growth_comparison/rce_osa.RDS"
if (!runmodels && file.exists(rce_osa_cache)) {
  rce_osa <- readRDS(rce_osa_cache)
} else {
  ss_osa <- Rceattle::fit_mod(data_list = simData, inits = NULL,
                              estimateMode = 1,
                              growthFun = growthFun_vb, # Von Bert, growth SD fixed to match WHAM
                              random_rec = TRUE, msmMode = 0,
                              initMode = 1, # Equilibrium -- same matched config as ss_est (init_dev fixed at 0)
                              # comp_offset = 0 -> WHAM-style multinomial (matches the fitting
                              # AND the OSA obsvec to WHAM; default 1e-5 would not match exactly)
                              fit_control = fit_control(phase = TRUE, verbose = 0,
                                                        comp_offset = 0))
  rce_osa <- Rceattle::osa_residuals(ss_osa, source = c("index", "catch", "comp", "caal"))
  saveRDS(rce_osa, rce_osa_cache)
}


# WHAM ----
library(wham)
# * Load model ----
wham_model <- readRDS(file = "tests/comparison/WHAM_growth_comparison/wham_model.RDS")


# WHAM OSA residuals on the existing fit. Growth-fork `$osa` type labels:
# "logcatch"/"logindex" (aggregate), "catchpal"/"indexpal" (length comp),
# "catchcaal"/"indexcaal" (conditional age-at-length).
# Cache the (slow) composition OSAs to a gitignored RDS; delete it (or set
# runmodels = TRUE) to recompute after the model changes.
wham_osa_cache <- "tests/comparison/WHAM_growth_comparison/wham_osa.RDS"
if (!runmodels && file.exists(wham_osa_cache)) {
  wham_osa <- readRDS(wham_osa_cache)
} else {
  wham_osa_fit <- wham::make_osa_residuals(
    wham_model, osa.opts = list(method = "oneStepGaussianOffMode", parallel = TRUE))
  wham_osa <- wham_osa_fit$osa          # = input$data$obs + a `residual` column
  saveRDS(wham_osa, wham_osa_cache)
}
print(unique(wham_osa$type))            # confirm the labels for your WHAM version


# Compare ----
# * Aggregate index & catch (one lognormal residual per year) ----
wham_agg_type <- c(catch = "logcatch", index = "logindex")
par(mfrow = c(1, 2))
for (src in c("index", "catch")) {
  rce <- rce_osa[rce_osa$source == src, c("year", "residual")]
  wh  <- wham_osa[wham_osa$type == wham_agg_type[[src]], c("year", "residual")]
  wh$year <- wh$year + simData$styr - 1
  m   <- merge(rce, wh, by = "year", suffixes = c(".rce", ".wham"))
  if (nrow(m) == 0) next
  plot(m$residual.wham, m$residual.rce,
       xlab = paste("WHAM", src, "OSA residual"),
       ylab = paste("Rceattle", src, "OSA residual"),
       main = sprintf("%s OSA  (r = %.3f)", src,
                      stats::cor(m$residual.wham, m$residual.rce)))
  abline(0, 1)
}

# * Length composition (by fleet) ----
# This is a length/growth model, so Rceattle "comp" is length composition and
# "caal" is conditional age-at-length. IMPORTANT pairing details:
#   (a) Rceattle's single "comp" source contains BOTH fleets (fleet 1 = survey,
#       fleet 2 = fishery); WHAM splits them into "indexpal" / "catchpal", so
#       filter Rceattle by fleet before merging or the fleets get conflated.
#   (b) WHAM stores the length BIN as its VALUE (2, 4, ..., 130); Rceattle stores
#       the bin INDEX (1..65), so convert with bin_value = 2 * age_length_bin.
par(mfrow = c(1, 2))
comp_pairs <- list(survey = list(fleet = 1, type = "indexpal"),
                   fishery = list(fleet = 2, type = "catchpal"))
for (nm in names(comp_pairs)) {
  p   <- comp_pairs[[nm]]
  rce <- rce_osa[rce_osa$source == "comp" & rce_osa$fleet == p$fleet, ]
  rce$bin <- 2 * rce$age_length_bin
  wh  <- wham_osa[wham_osa$type == p$type, ]
  wh$year <- wh$year + simData$styr - 1
  wh$bin  <- as.numeric(wh$bin)
  m <- merge(rce[, c("year", "bin", "residual")],
             wh[,  c("year", "bin", "residual")],
             by = c("year", "bin"), suffixes = c(".rce", ".wham"))
  m <- m[is.finite(m$residual.rce) & is.finite(m$residual.wham), ]
  plot(m$residual.wham, m$residual.rce,
       xlab = paste("WHAM", nm, "lencomp OSA"), ylab = paste("Rceattle", nm, "lencomp OSA"),
       main = sprintf("%s lencomp  (r = %.4f)", nm, stats::cor(m$residual.wham, m$residual.rce)))
  abline(0, 1)
}

# * Conditional age-at-length (CAAL) ----
# WHAM labels each CAAL bin "age_lengthvalue" (e.g. "1_8" = age 1 at length 8);
# Rceattle stores age in `age_length_bin` and the conditioning length BIN INDEX in
# `length`, so the WHAM length value = 2 * Rceattle `length`.
par(mfrow = c(1, 1))
rce <- rce_osa[rce_osa$source == "caal", ]
rce$age <- rce$age_length_bin
rce$lval <- 2 * rce$length                       # WHAM length value = 2 * bin index
wh  <- wham_osa[wham_osa$type == "indexcaal", ]
wh$year <- wh$year + simData$styr - 1
sp  <- do.call(rbind, strsplit(as.character(wh$bin), "_"))
wh$age  <- as.integer(sp[, 1]); wh$lval <- as.integer(sp[, 2])
m <- merge(rce[, c("year", "age", "lval", "residual")],
           wh[,  c("year", "age", "lval", "residual")],
           by = c("year", "age", "lval"), suffixes = c(".rce", ".wham"))
m <- m[is.finite(m$residual.rce) & is.finite(m$residual.wham), ]
plot(m$residual.wham, m$residual.rce,
     xlab = "WHAM CAAL OSA residual", ylab = "Rceattle CAAL OSA residual",
     main = sprintf("survey CAAL  (r = %.4f)", stats::cor(m$residual.wham, m$residual.rce)))
abline(0, 1)

