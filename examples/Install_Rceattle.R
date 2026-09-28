
# Install from GitHub ----
# - Recommended

# 1. Install devtools and TMB
install.packages("devtools")
install.packages("TMB", type = "source")
# (See https://github.com/kaskr/adcomp/wiki/Download for TMB toolchain setup.)
# Try "TMB::runExample(all = TRUE)" to confirm TMB compiles on your machine.

# 2. Install Rceattle. CRAN dependencies are pulled automatically.
devtools::install_github("afsc-assessments/Rceattle", build_vignettes = TRUE)

# 3. (Optional) TMBhelper provides richer optimization diagnostics.
# Rceattle uses it transparently if present, otherwise falls back to plain
# nlminb + sdreport.
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
