#if defined(__clang__)
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <TMB.hpp>
#if defined(__clang__)
# pragma clang diagnostic pop
#endif

// GCC/Eigen emits false -Warray-bounds positives where its vectorized reductions
// inline TMB's multi-dimensional indexing (e.g. the 5-D suitability() access in
// predation.hpp). A source pragma rather than a -Wno-array-bounds compiler flag,
// so CRAN's compilation-flags check stays clean and other diagnostics still show.
#if defined(__GNUC__)
# pragma GCC diagnostic ignored "-Warray-bounds"
#endif

#include "helper_functions.hpp"
#include "comp_osa.hpp"
#include "comp_sim.hpp"
#include "growth.hpp"
#include "selectivity.hpp"
#include "recruitment.hpp"
#include "spr.hpp"
#include "bioenergetics.hpp"
#include "predation.hpp"
#include "diet_data.hpp"
#include "linkage.hpp"
#include "likelihood.hpp"

/**
 * @file ceattle.cpp
 * @brief The CEATTLE objective function: configuration, population dynamics,
 *   predation, likelihood and report. The numbered section index is in the
 *   banner below; the process equations live in the headers this file includes.
 */

// Reads an R list() of matrices into a vector<matrix<Type>>: the per-fleet
// survey-index covariance matrices Sigma supplied when Index_loglike is
// "MVN"/"MVNORM" (the AMAK/ebswp DoCovBTS covariance survey likelihood; sections
// 8.2b and 13.1). Other fleets pass a 1x1 inert dummy, so the list is always
// length n_flt and indexes by fleet code.
template<class Type>
struct LOM_t : vector<matrix<Type> > {
  LOM_t(SEXP x) {  // x = R list of matrices
    (*this).resize(LENGTH(x));
    for(int i = 0; i < LENGTH(x); i++){
      SEXP m = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(m);
    }
  }
};

/* ------------------------------------------------------------------------- //
 * CEATTLE: Bioenergetic-based Assessment for Understanding Biomass Linkages To
 * The Environment -- a multispecies statistical model in Template Model
 * Builder. Version: see DESCRIPTION.
 *
 * Follows the Google C++ style guide
 * (https://google.github.io/styleguide/cppguide.html#General_Naming_Rules),
 * except that free functions stay snake_case.
 *
 * CITATIONS:
 * 1. Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2015). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134, 360–378. https://doi.org/10.1016/j.dsr2.2015.08.001
 * 2. Adams, G. D., Holsman, K. K., Barbeaux, S. J., Dorn, M. W., Ianelli, J. N., Spies, I., ... & Punt, A. E. (2022). An ensemble approach to understand predation mortality for groundfish in the Gulf of Alaska. Fisheries Research, 251, 106303.
 * 3. Wassermann, S. N., Adams, G. D., Haltuch, M. A., Kaplan, I. C., Marshall, K. N., & Punt, A. E. (2025). Even low levels of cannibalism can bias population estimates for Pacific hake. ICES Journal of Marine Science, 82(1), fsae064.
 * ------------------------------------------------------------------------- //
 *  INDEX:
 *  0. Load dependencies
 *  1. Model configuration
 *  2. Model inputs
 *  3. Model parameters
 *  4. Derived quantities
 *  5. Initial calculations
 *  6. Population dynamics
 *  7. Predation mortality equations
 *  8. Index components equations
 *  9. Fishery components equations
 *  10. Composition equations
 *  11. Predicted stomach content
 *  12. Derived quantities
 *  13. Likelihood components
 *  14. Report section
 *  15. Model return/end
 *
 */


template<class Type>
Type objective_function<Type>::operator() () {
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE

  /** ------------------------------------------------------------------------ //
   * 1. MODEL CONFIGURATION                                                    //
   * ------------------------------------------------------------------------- */
  // 1.1. CONFIGURE MODEL (this section sets up the switches)
  DATA_INTEGER(estimateMode);             // 0 = hindcast + HCR projection, 1 = hindcast, 2 = projection from inits, 3 = build only, 4 = smoke test
  DATA_INTEGER(msmMode);
  //    0 = single species (no predation mortality)
  //    1 = Type II MSVPA based (sensu Holsman et al 2015)
  //    2 = Type III MSVPA based
  //    3-9 = not yet implemented (Holling forms, predator interference and
  //          preemption, Hassell-Varley, Ecosim)
  DATA_IVECTOR(suitMode);                 // Estimate suitability
  // DATA_INTEGER(est_diet);              // Include diet data in the likelihood
  // DATA_INTEGER(avgnMode);              // Not read; predation uses avgN_at_age
  DATA_INTEGER(initMode);                 // How the age-structure is initialized
  DATA_IVECTOR(forecast);                 // Run the model in forecast mode


  // 1.2. Temporal dimensions
  DATA_INTEGER( styr );                   // Start year
  DATA_INTEGER( endyr );                  // End of estimation years
  DATA_INTEGER( projyr );                 // End year of projection

  DATA_INTEGER( srr_mse_switchyr );             // Last year of the average-recruitment window; for MSE runs.
  DATA_IVECTOR( suit_styr );              // First year of the suitability average, per predator species.
  DATA_IVECTOR( suit_endyr );             // Last year of the suitability average, per predator species.
  DATA_INTEGER( srr_hat_styr );           // First year of the stock-recruitment-penalty / env-rec window.
  DATA_INTEGER( srr_hat_endyr );          // Last year of the stock-recruitment-penalty / env-rec window.

  int nyrs = projyr - styr + 1;
  int nyrs_hind = endyr - styr + 1;

  // suit_styr / suit_endyr are per-predator (offset to model-year indices below,
  // once `nspp` is available).
  int nyrs_srrmean = srr_mse_switchyr - styr + 1;

  srr_hat_styr = srr_hat_styr - styr;
  srr_hat_endyr = srr_hat_endyr - styr;
  if(nyrs_srrmean > nyrs_hind){nyrs_srrmean = nyrs_hind;}
  if(srr_hat_styr == 0){srr_hat_styr = 1;} // R_hat starts at year 1.

  // 1.3. Number of species
  DATA_INTEGER( nspp );                   // Number of species (prey)

  // Offset per-predator suitability reference years to model-year indices and
  // count the averaging window length for each predator species.
  vector<int> nyrs_suit(nspp);
  for(int sp = 0; sp < nspp; sp++){
    suit_endyr(sp) = suit_endyr(sp) - styr;
    suit_styr(sp)  = suit_styr(sp) - styr;
    nyrs_suit(sp)  = suit_endyr(sp) - suit_styr(sp) + 1;
  }

  DATA_IVECTOR( pop_wt_index );           // Dim 1 of weight_obs to use for population dynamics
  DATA_IVECTOR( ssb_wt_index );           // Dim 1 of weight_obs to use for spawning stock biomass
  DATA_IVECTOR( pop_age_transition_index );// Dim 1 of age_trans_matrix to use for population dynamics
  DATA_IVECTOR( estDynamics );            // 0 = population parameters estimated, 1 = numbers-at-age provided, 2 = provided numbers-at-age times an estimated scalar
  pop_wt_index -= 1;                      // Indexing starts at 0
  ssb_wt_index -= 1;                      // Indexing starts at 0
  pop_age_transition_index -= 1;          // Indexing starts at 0


  // 1.4. RECRUITMENT SETTINGS
  DATA_INTEGER(srr_fun);                  // Stock recruit relationship for hindcast estimation
  DATA_INTEGER(srr_pred_fun);             // Stock recruit relationship for projection/brps/penalty
  DATA_INTEGER(proj_mean_rec);
  //    0 = project recruitment using log_R0 and rec devs
  //    1 = project recruitment using mean rec (can also have adjusted rec devs)
  DATA_INTEGER(srr_est_mode);             // Logical of wether to add normal prior to stock recruit-relationship
  DATA_VECTOR(srr_prior);                 // Prior mean for stock recruit relationship parameter
  DATA_VECTOR(srr_prior_sd);              // Prior sd for stock recruit relationship parameter
  DATA_INTEGER( niter );                  // Number of loops for MSM mode
  DATA_VECTOR( Bmsy_lim );                // Upper limit for Bmsy in ricker function. Will add penalty if 1/beta > lim


  // 1.5. HARVEST CONTROL RULE (HCR) SETTINGS
  DATA_INTEGER(HCR);                      // Function to be used for harvest control rule
  DATA_INTEGER(DynamicHCR);               // TRUE/FALSE. Whether to use static or dynamic reference points (default = FALSE)
  DATA_VECTOR(Ptarget);                   // Target SSB as a PROPORTION of static or dynamic SSB at F = 0 (B40% is 0.40)
  DATA_VECTOR(Plimit);                    // Limit SSB as a PROPORTION of static or dynamic SSB at F = 0 (B10% is 0.10)
  DATA_VECTOR(Ftarget_percent);           // Proportion of SSB-per-recruit or SB0 at F = 0 that sets the target F (F40% is 0.40)
  DATA_VECTOR(Flimit_percent);            // Proportion of SSB-per-recruit or SB0 at F = 0 that sets the limit F (F35% is 0.35)
  DATA_VECTOR(Alpha);                     // Parameter used in NPFMC Tier 3 HCR
  DATA_VECTOR(Fmult);                     // Multiplier for target fishing mortality (Fmult * Ftarget). For example 0.75 * F40%.
  DATA_VECTOR(QnormHCR);                  // Add to Flimit to set Ftarget based on the Pstar approach: Flimit + qnorm(Pstar, 0, Sigma)


  // 1.6. MODEL OBJECTS
  // 1.6.1. LOOPING INDICES -- k = observation, sp = species, sex = sex (0 = combined; 1 = females; 2 = males)
  // age = age, ln = length, yr = year
  int sp, sex, age, ln, yr;
  int index, flt;                                                         // Survey and fishery indices
  int flt_yr, flt_sex, comp_type;
  int fsh_ind, index_ind, comp_ind, yr_ind;                       // Indices for survey sets
  int wt_idx_pop, wt_idx_ssb, wt_idx_flt; // Indices for weight indices
  Type mo = 0;                                                           // Month float
  if(msmMode == 0) { niter = 1; }                                        // Number of iterations for SS mode


  /** ------------------------------------------------------------------------ //
   * 2. MODEL INPUTS                                                           //
   * ------------------------------------------------------------------------- */

  // -- 2.1. Species attributes
  DATA_IVECTOR( nsex );                   // Sexes modelled per species: 1 = combined, 2 = both.
  int max_sex = imax(nsex);               // Max sexes, for array dimensions
  DATA_VECTOR( spawn_month );             // Month of spawning to adjust mortality
  DATA_IVECTOR( nages );                  // Age bins per species.
  DATA_IVECTOR( minage );                 // Minimum age of each species
  DATA_IVECTOR( nlengths );               // Observation length bins per species.
  DATA_MATRIX(lengths);                   // Length bins for each species [sp, nlengths]
  DATA_IVECTOR( nlengths_pop );           // Population length bins per species. Equals nlengths when no population grid is given
  DATA_MATRIX( lengths_pop );             // Population length bins, lower edges (cm) [sp, nlengths_pop]: the grid the age-length key, weight and maturity are integrated on
  DATA_IMATRIX( pop_to_data_bin );        // Data length bin (0-based) that each population bin sums into [sp, nlengths_pop]
  DATA_ARRAY( NByageFixed );              // Input numbers- or index-at-age, multiplied (or not) by pop_scalar to get N_at_age
  DATA_VECTOR( MSSB0 );                   // SB0 from projecting forward in multi-species mode under no fishing
  DATA_VECTOR( MSB0 );                    // B0 from projecting forward in multi-species mode under no fishing

  int max_nlengths = imax(nlengths);      // Max nlengths, for array dimensions
  // The age-length key and selectivity-at-length live on the POPULATION length bins,
  // which are the data bins whenever no pop_lengths is supplied.
  int max_nlengths_pop = imax(nlengths_pop);

  // -- 2.3c. First and last population bin in each data bin, per species. SS3 forms
  // sel(L) * P(L|age) at population resolution and bins the result, so pred_CAAL and the
  // length compositions sum each data bin's run rather than applying one bin-average
  // selectivity across it. The map is monotone, so the runs are contiguous; with one grid
  // every run is one bin.
  matrix<int> pop_bin_lo(nspp, max_nlengths); pop_bin_lo.setZero();
  matrix<int> pop_bin_hi(nspp, max_nlengths); pop_bin_hi.setZero();
  for(int sp_i = 0; sp_i < nspp; sp_i++){
    for(int d = 0; d < nlengths(sp_i); d++){ pop_bin_lo(sp_i, d) = -1; pop_bin_hi(sp_i, d) = -2; }
    for(int lp = 0; lp < nlengths_pop(sp_i); lp++){
      int d = pop_to_data_bin(sp_i, lp);
      if(pop_bin_lo(sp_i, d) < 0) pop_bin_lo(sp_i, d) = lp;
      pop_bin_hi(sp_i, d) = lp;
    }
  }
  int max_age = imax(nages);              // Max nages, for array dimensions
  // Sized on the POPULATION bins: the non-parametric curves are built over nlengths_pop,
  // so a finer grid would otherwise write past non_par_sel.
  int max_bin = (max_age > max_nlengths_pop) ? max_age : max_nlengths_pop;

  // -- 2.2. M1_at_age specifications
  DATA_IVECTOR(M1_model);
  DATA_IVECTOR(M1_re);
  DATA_IVECTOR(M1_use_prior);
  DATA_IVECTOR(M2_use_prior);
  DATA_VECTOR(M_prior);
  DATA_VECTOR(M_prior_sd);

  // -- 2.3. Growth model specifications
  DATA_IVECTOR(growth_model); // 0: "input", 1: "vB-classic", 2: "Richards", 3: "nonparametric LAA" [sp]
  DATA_IVECTOR(growth_sd_style); // Plus-group SD-at-age treatment [sp]: 1 = WHAM (pin to exp(sd_Linf)), 2 = interpolate by length
  DATA_IVECTOR(growth_sd_form);  // Growth variability endpoints [sp]: 1 = SDs in cm (SS3 CV_Growth_Pattern 2), 2 = CVs, SD = CV * L (SS3 pattern 0)
  DATA_IVECTOR(growth_plus_length); // Plus-group mean length [sp]: 1 = M1-weighted, 2 = none (SS3 Linf_decay -998), 3 = SS3.24 (-999), 4 = SS3 decay rate
  DATA_VECTOR(plus_group_decay); // Decay rate (per year) for growth_plus_length == 4 [sp]
  DATA_IMATRIX(sel_dn6_ends);    // DoubleNormalSS3, per fleet [flt, 2]: 1 = scale the initial / final end by P5 / P6, 0 = SS3's -999 (unscaled)
  DATA_IVECTOR(mat_len_use);     // 1 = maturity-at-length replaces the age-based maturity in spawning output [sp]
  DATA_MATRIX(mat_len_pars);     // Maturity-at-length [sp, 2]: L50 (cm), logistic slope (per cm)
  DATA_VECTOR(growth_age_L1); // VB anchor age (= SS3 Growth_Age_for_L1) per sp; defaults to max(0.5, minage[sp]) in R-side fit_mod()

  // -- 2.3b. Linkage table (R/0-linkage_encode.R): covariate and year effects on a
  //          process parameter p in one stratum (species, sex, age bin, fleet; 0 = all).
  //          Row i adds beta_i * x_i(y), with x_i(y) = linkage_X(y, linkage_X_col(i)) its
  //          design column in model year y (row 0 = styr; offset 0 past the last row).
  //          With S_log, S_id and S_exp the sums of beta_i * x_i(y) over p's log-,
  //          identity- and exponential-link rows, and p's own deviates (log_M1_dev,
  //          index_q_dev, sel_inf_dev, log_sel_slp_dev) still applied:
  //            most parameters:        p_y = exp(log p + S_log) + S_id
  //            sel_inf, DoubleNormalSS3 P1-P6:  p_y = (p + S_id) * exp(S_log), stored scale
  //            catchability (5.3):     q_y = exp((log q + S_log) * exp(S_exp) + index_q_dev) + S_id
  //          (Intercept) rows are held at 0, so p carries the level, except R_init's,
  //          which is estimated: a log multiplier on the initial age structure, not on
  //          year-1 recruitment (6.4). Composition and growth-sd rows add no offset;
  //          an (Intercept)'s init, bounds and prior go to the base parameter.
  //          Terms used below:
  //            row    one beta_i; fixed-effect rows read beta_linkage.
  //            slot   one random-effect deviation, e.g. one year of (1 | Year) in one
  //                   stratum; its row's beta is beta_linkage_re_all(slot) (5.2b), times
  //                   beta_linkage_obs if observed.
  //            position  where a slot sits in beta_linkage_re (integrated) or
  //                   beta_linkage_re_pen (penalized): linkage_re_slot.
  //            group  the slots of one random term (re_group, e.g. Year, and structure)
  //                   on one parameter and stratum: one sd, log_sigma_linkage, and one
  //                   structure, IID, random walk or AR1. An observed AR1 is QAR1 (Rogers
  //                   et al. 2024), its latent also fitted to an env_data series.
  //          Empty when no build_*() supplied `linkages`.
  DATA_IVECTOR(linkage_process);       // process code, RCEATTLE_PROC_* (linkage.hpp): recruitment, M, growth, q, sel, comp
  DATA_IVECTOR(linkage_param);         // per-process parameter code
  DATA_IVECTOR(linkage_species);       // 1-based sp id; 0 = all
  DATA_IVECTOR(linkage_sex);           // 1-based sex id; 0 = all
  DATA_IVECTOR(linkage_age_bin);       // 1-based age bin, age - minage + 1; 0 = all
  DATA_IVECTOR(linkage_fleet);         // 1-based Fleet_code; 0 = all
  DATA_IVECTOR(linkage_X_col);         // 0-based column of linkage_X
  DATA_IVECTOR(linkage_link);          // link code, RCEATTLE_LINK_* (linkage.hpp)
  // Only an identity-link recruitment offset can make the curve non-positive, so
  // the recruitment floors (sections 6.3, 6.5, 6.6 and 6.9) run only then; every
  // other model's AD tape is left exactly as it was.
  int rec_floor_on = 0;
  for (int i = 0; i < linkage_link.size(); ++i) {
    if (linkage_process(i) == RCEATTLE_PROC_RECRUIT && linkage_link(i) == RCEATTLE_LINK_IDENTITY) rec_floor_on = 1;
  }
  DATA_IVECTOR(linkage_re_index);      // per row: -1 = fixed; else its 0-based slot in beta_linkage_re_all
  // linkage_re_sigma, _rho and _obs give positions in log_sigma_linkage, trans_rho_linkage
  // and beta_linkage_obs; one sd per group, so linkage_re_sigma is the slot's group.
  DATA_IVECTOR(linkage_re_sigma);      // per RE slot: its 0-based group (an index, not an sd)
  DATA_IVECTOR(linkage_re_integrate);  // per RE slot: 1 = integrated out (Laplace), 0 = penalized
  DATA_IVECTOR(linkage_re_slot);       // per RE slot: its position (2.3b)
  DATA_IVECTOR(linkage_re_struct);     // per RE group: RCEATTLE_RE_IID (a plain (1 | Year)), _RW or _AR1 (linkage.hpp)
  DATA_IVECTOR(linkage_re_rho);        // per RE group: 0-based position in trans_rho_linkage (ar1 only; -1 otherwise)
  DATA_IVECTOR(linkage_re_sigma_prior_family); // per RE group: prior on the SD (0=none,1=normal,2=lognormal,3=gamma,4=beta)
  DATA_IVECTOR(linkage_re_rho_prior_family);   // per RE group: prior on rho (natural (-1,1) scale)
  DATA_VECTOR(linkage_re_rho_prior_p1);        // per RE group: rho prior param 1
  DATA_VECTOR(linkage_re_rho_prior_p2);        // per RE group: rho prior param 2
  DATA_IVECTOR(linkage_re_obs);        // per RE group: 0-based position in beta_linkage_obs (Rogers QAR1 observed ar1; -1 otherwise)
  // linkage_re_obs_sd is unread here; log_obs_sd_linkage starts at its log.
  DATA_VECTOR(linkage_re_obs_sd);      // per RE group: the obs_sd supplied (0 if unobserved)
  DATA_VECTOR(linkage_re_obs_value);   // per RE slot: observed covariate value (0 if unobserved)
  DATA_IVECTOR(linkage_re_obs_mask);   // per RE slot: 1 if the covariate is observed that year, 0 otherwise
  DATA_IVECTOR(linkage_re_obsvec_idx); // per RE slot: 0-based obsvec position for OSA (-1 if unobserved / not built)
  DATA_VECTOR(linkage_re_sigma_prior_p1);      // per RE group: prior param 1
  DATA_VECTOR(linkage_re_sigma_prior_p2);      // per RE group: prior param 2
  DATA_IVECTOR(linkage_is_intercept);  // 1 if design_col == "(Intercept)", 0 otherwise
  DATA_IVECTOR(linkage_prior_family);  // none=0, normal=1, lognormal=2, gamma=3, beta=4
  DATA_VECTOR(linkage_prior_p1);       // family-specific prior param 1
  DATA_VECTOR(linkage_prior_p2);       // family-specific prior param 2
  DATA_MATRIX(linkage_X);              // design matrix [model year from styr, design column]

  // -- 2.4. Fleet controls (i.e. how to assign data to objects)
  DATA_IVECTOR(flt_type);                 // Fleet_type: 0 = Off, 1 = Fishery, 2 = Survey (NOT a boolean)
  DATA_VECTOR(flt_month);
  DATA_IVECTOR(flt_sel_type);             // Fleet selectivity parameterization
  DATA_IVECTOR(flt_sel_dim);              // Fleet selectivity dimension (0 = age, 1 = length)
  DATA_IVECTOR(flt_n_sel_bins);           // Age/length bins for non-parametric selectivity
  DATA_IVECTOR(flt_sel_cap_bin);          // NonParametricPM (type 9): first bin (0-based, age or length per flt_sel_dim) at/after which realized selectivity is held flat (RTMB cap_old_age). < 0 = no cap.
  DATA_IVECTOR(flt_varying_sel);          // Whether time-varying selectivity is estimated
  DATA_IVECTOR(flt_spp);                  // Fleet species
  DATA_IVECTOR(bin_first_selected);       // First selected age/length BIN, 0-based (selectivity below it = 0)
  DATA_IVECTOR(sel_norm_bin1);            // Normalization BIN, 0-based, or a control flag when negative (-999 = do not normalize, -99 = normalize by the maximum); rearrange_data() converts, so NA never arrives. LogisticPM (type 11): lower bin of the selectivity-penalty range.
  DATA_IVECTOR(sel_norm_bin2);            // Upper normalization BIN, same encoding as sel_norm_bin1. LogisticPM (type 11): upper bin of the selectivity-penalty range.
  DATA_IVECTOR(sel_norm_scope);           // Whether normalization pools its reference across sexes: 0 = WithinSex (each sex divided by its own reference; both reach 1), 1 = AcrossSexes (one pooled reference; the lower sex stays below 1). Orthogonal to sel_norm_bin1/2, which say WHERE the reference is taken.
  DATA_IVECTOR(flt_sel_start_yr);         // Per-fleet selectivity start year (0-based from styr); selectivity penalties start the year after this. Default 0 (= styr).
  DATA_IVECTOR(flt_sel_pen_first_bin);    // Per-fleet first bin (0-based, age or length per flt_sel_dim) for the non-parametric shape/monotonicity penalty. < 0 -> defaults to bin_first_selected. Lets the shape constraint span a narrower range than the (possibly non-zero) first selected bin (e.g. ATS mina_ats > first selected age).
  DATA_IVECTOR(flt_sel_lead);             // 1 if this fleet's selectivity penalty should be accumulated; 0 if it mirrors an earlier fleet's selectivity (same Selectivity_index + type) so the shared penalty is counted once.
  DATA_IVECTOR(flt_q_lead);               // As flt_sel_lead, for catchability: 1 if this fleet carries the q prior / deviate penalties, 0 if it shares an earlier fleet's Q_index so the shared block is counted once.
  DATA_IVECTOR(flt_sel_pen_last_bin);     // Per-fleet last bin (0-based, = left bin of the last adjacent pair) for the non-parametric shape penalty. < 0 -> defaults to nbins-2 (whole range).
  DATA_IVECTOR(flt_sel_shape_mode);       // Non-parametric shape-penalty mode: 0 = directional (sign of Sel_curve_pen1 -> penalize decreasing/increasing, one-sided, ADMB/AMAK); 1 = smooth (two-sided d^2 over adjacent ages, RTMB "rpm").
  DATA_VECTOR(flt_sel_avgsel_pen);        // Per-fleet weight on the AMAK "avgsel" base-level penalty (type 9 only); the form is penalty (5) in section 13. 0 = off (default).
  DATA_IVECTOR(comp_ll_type);             // Composition log-likelihood type
  DATA_IVECTOR(comp_accum_young);         // Per-fleet age/length-composition young-tail accumulation bin (1-based ordinal on the comp dimension); bins below it fold into it. 1 (or <1) = no young accumulation (AFSC ac_yng).
  DATA_IVECTOR(comp_accum_old);           // Per-fleet age/length-composition old-tail accumulation bin (1-based); bins above it fold into it. 0/NA or >= nbins = no old accumulation (AFSC ac_old).
  DATA_IVECTOR(caal_ll_type);             // CAAL composition log-likelihood type
  DATA_IVECTOR(flt_units);                // Fleet units (1 = weight, 2 = numbers)
  DATA_IVECTOR(flt_wt_index);             // Vector to save 1st dim of weight to use for weight-at-age
  DATA_IVECTOR(flt_age_transition_index); // Dim 1 of age_trans_matrix each fleet reads for its ALK
  DATA_IVECTOR(flt_ageing_error_index);   // Ageing error matrix each fleet reads (0-based); the fleet's own species unless fleet_control gives one
  DATA_IVECTOR(est_index_q);              // Catchability form per fleet; see q_map and sections 8.2/8.2b
  DATA_IVECTOR(index_varying_q);          // Whether time-varying q is estimated
  DATA_IVECTOR(est_sigma_index);          // Whether survey sigma is estimated
  DATA_IVECTOR(est_sigma_fsh);            // Whether fishery sigma is estimated

  // -- 2.4.1 Fishery Components
  DATA_IMATRIX( catch_ctl );              // Fishery catch metadata; columns = Fleet_code, Species, Year
  DATA_IMATRIX( catch_n );                // Info for fishery biomass; columns = Month
  DATA_MATRIX( catch_obs );               // Observed fishery catch (mt or thousands of fish per the fleet's Observation_units) and log_sd; columns = Observation, Error

  // -- 2.4.2 Survey components
  DATA_IMATRIX( index_ctl );              // Index metadata; columns = Fleet_code, Species, Year
  DATA_MATRIX( index_n );                 // Info for index; columns = Month
  DATA_MATRIX( index_obs );               // Observed index and log_sd; columns = Observation, Error
  DATA_VECTOR( index_log_q_prior );        // Prior mean for catchability
  DATA_IVECTOR( index_ll_type );          // Survey index likelihood family per fleet (0 = lognormal IID, 1 = MVN bare quadratic form, 2 = MVNORM full density, 3 = natural-scale normal with absolute sd)
  DATA_STRUCT( index_cov_mat, LOM_t );    // Per-fleet covariance matrices Sigma for the MVN/MVNORM survey likelihood (1x1 dummy if unused)
  DATA_VECTOR( index_cov_const );         // Per-fleet 0.5*(logdet(Sigma) + n*log(2pi)); subtracted for "MVN" so it reports the bare quadratic form

  // -- 2.4.2b One-step-ahead (OSA) residual support
  // `obsvec` is every observation entering the likelihood, flat: log catch, log index,
  // each comp/caal row's bin counts, each stomach's diet composition. `keep` is the
  // indicator TMB::oneStepPredict() toggles element by element; it defaults to all ones,
  // so ordinary fitting is unchanged. Each `*_obsvec_idx` gives its `*_obs` row's 0-based
  // position in `obsvec` (for a composition, where the row's bins start), or -1 when the
  // row is not fitted. osa_mode 1 switches the comp/caal/diet branches to a proper,
  // unweighted, keep-gated density; mode 2 adds the conditional CDF terms
  // oneStepPredict(method = "cdf") reads -- log F(x) gated by keep.cdf_lower,
  // log(1 - F(x)) by keep.cdf_upper -- computed only in mode 2, so no other path pays
  // for the pbeta/pnorm. The Dirichlet-multinomial families have no CDF term
  // (comp_osa.hpp) and osa_residuals() keeps those fleets off this method.
  DATA_VECTOR( obsvec );                    // Flat observations for OSA residuals
  DATA_VECTOR_INDICATOR( keep, obsvec );    // oneStepPredict keep indicator (defaults to 1 when fitting)
  DATA_IVECTOR( catch_obsvec_idx );         // obsvec position for each catch_obs row (-1 = excluded)
  DATA_IVECTOR( index_obsvec_idx );         // obsvec position for each index_obs row (-1 = excluded)
  DATA_IVECTOR( comp_obsvec_idx );          // obsvec start position for each comp_obs row's bins (-1 = excluded)
  DATA_IVECTOR( caal_obsvec_idx );          // obsvec start position for each caal_obs row's bins (-1 = excluded)
  DATA_IVECTOR( diet_obsvec_idx );          // obsvec start position for each stomach's prey bins (incl. "other prey"); length n_stomach_obs (-1 = excluded)
  DATA_INTEGER( osa_mode );                 // 0 = normal fitting (default); 1 = OSA build (unweighted keep-gated comp/caal/diet densities); 2 = OSA build + conditional CDF terms for method = "cdf"
  DATA_INTEGER( adreport_sel );             // 0 = no selectivity standard errors (default); 1 = ADREPORT log_sel_at_age; set via fit_control(selectivity_se =)
  DATA_SCALAR( comp_offset );               // proportion offset added to comp/caal/diet obs & pred before the multinomial; set via rearrange_data()/fit_control()

  // -- 2.4.2c. Simulation switches (sim_mod(simulate = TRUE))
  // Observations are always drawn; process error is a choice, so it is switched.
  //
  // simulate_state: which process. Slots 0 recruitment, 1 M, 2 growth, 3 catchability,
  // 4 selectivity -- the codes linkage_process uses, so a random linkage gates on
  // simulate_state(linkage_process(i)) untranslated. Default off: redrawing a process
  // changes what self_test() measures, from recovering parameters to recovering a
  // process. Recruitment covers init_dev too, the initial age structure being recruitment
  // pre-styr.
  //
  // simulate_period: slot 0 the fitted window (0 < Year <= endyr), slot 1 outside it.
  // Only slot 0 is read, and only by the PROCESS draws, whose densities cover the
  // hindcast alone. Observations are drawn in every period deliberately, because
  // run_mse() splices the negative-Year rows back in as the next assessment's data. Slot
  // 1 is reserved and inert; sim_mod() exposes neither slot.
  DATA_IVECTOR( simulate_state );
  DATA_IVECTOR( simulate_period );

  // -- 2.4.3. Composition data
  DATA_IMATRIX( comp_ctl );               // Age/length comp metadata; columns = Fleet_code, Species, Sex, Age0_Length1, Year
  DATA_MATRIX( comp_n );                  // Month and sample size on observed age/length comp; columns = Month, Sample size
  DATA_MATRIX( comp_obs );                // Observed age/length comp; cols = Comp_1, Comp_2, etc. can be proportion
  DATA_IMATRIX( caal_ctl );               // CAAL metadata; columns = Fleet_code, Species, Sex, Year, Length_bin
  DATA_MATRIX( caal_n );                  // Sample size on CAAL; ONE column = Sample size. Unlike comp_n there is no month: a CAAL observation is placed at its fleet's Month, and the age-length key is annual.

  // -- 2.4.4. Initial equilibrium catch: the catch yielded under the initial F,
  // in the year before the hindcast. 0 rows unless the data supply one.
  DATA_IMATRIX( equil_catch_ctl );        // columns = Fleet_code, Species
  DATA_MATRIX( equil_catch_obs );         // columns = Catch, Log_sd
  DATA_MATRIX( caal_obs );                // Observed CAAL; cols = Comp_1, Comp_2, etc. can be proportion

  // -- 2.4.5 Age and selectivity
  DATA_IMATRIX( emp_sel_ctl );            // Empirical selectivity metadata; columns = Fleet_code, Species, Sex, Year
  DATA_MATRIX( emp_sel_obs );             // Observed emprical fishery selectivity; columns = Compe_1, Comp_2, etc.
  DATA_ARRAY( age_trans_matrix);          // observed sp_age/size compositions; n = [nspp, nages, index_age_bins]
  DATA_ARRAY( age_error );                // Aging error matrices; n = [n_matrices, nages, nages]. One per species unless fleet_control$Ageing_error_index selects otherwise.

  // -- 2.3.5. Growth
  DATA_ARRAY( weight_obs );               // Weight-at-age by year; n = [nweight, sex, nages, nyrs]

  // 2.3.6. Diet data
  DATA_IVECTOR(diet_ll_type);             // Diet composition log-likelihood type
  DATA_VECTOR( fday );                    // number of foraging days for each predator
  DATA_ARRAY( ration_data );              // Relative-foraging rate
  DATA_MATRIX( diet_obs );                // Pred, prey, pred-age, prey-age for diet matrix (weight of prey in pred stomach)
  DATA_IMATRIX( diet_ctl );               // Per diet_obs row: predator, prey, both sexes, both ages, year
  DATA_INTEGER(n_stomach_obs);            // The total number of unique stomach samples (groups)
  DATA_IVECTOR(stomach_id);               // A vector mapping each diet data row to a stomach ID

  // 2.3.7. Environmental data
  DATA_MATRIX( env_index );               // Matrix of environmental predictors such as bottom temperature

  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR( other_food );              // Biomass of other prey (mt)
  DATA_VECTOR( Pvalue );                  // Proportion of Cmax; the P in Cmax*fT*Pvalue*PAge
  DATA_IVECTOR( Ceq );                    // Which consumption equation to use; currently all sp = 1
  DATA_IVECTOR( Cindex );                 // Which env_index column drives bioenergetics.
  DATA_VECTOR( CA );                      // Wt specific intercept of Cmax=CA*W^CB
  DATA_VECTOR( CB );                      // Wt specific slope of Cmax=CA*W^CB
  DATA_VECTOR( Qc );                      // used in fT, QC value
  DATA_VECTOR( Tco );                     // used in fT, thermal optimum
  DATA_VECTOR( Tcm );                     // used in fT, thermal max
  DATA_VECTOR( Tcl );                     // used in fT eq 3, limit
  DATA_VECTOR( CK1 );                     // used in fT eq 3, limit where C is .98 max (ascending)
  DATA_VECTOR( CK4 );                     // used in fT eq 3, temp where C is .98 max (descending)

  // -- 2.4.3. Others
  DATA_MATRIX( sex_ratio );               // Proportion-at-age of females of population
  DATA_MATRIX( maturity );                // Proportion of mature females at age; [nspp, nages]
  DATA_SCALAR(bias_adjust_obs);           // Whether to apply bias adjustment (-sigma^2/2) in index likelihoods
  DATA_SCALAR(bias_adjust_proc);           // Whether to apply bias adjustment (-sigma^2/2) in process likelihoods


  /** ------------------------------------------------------------------------ //
   * 3. PARAMETER SECTION                                                      //
   * ------------------------------------------------------------------------- */

  PARAMETER( dummy );                             // Tests derived quantities given input parameters; n = [1]
  PARAMETER_VECTOR( log_pop_scalar );              // Log multiplier on input numbers-at-age, per species (estDynamics = 2)

  // -- 3.1. Recruitment parameters
  PARAMETER_MATRIX( rec_pars );                   // Stock-recruit parameters: col1 = mean rec, col2 = SRR alpha, col3 = SRR beta
  PARAMETER_VECTOR( R_log_sd );                    // Standard deviation of recruitment deviations
  PARAMETER_MATRIX( rec_dev );                    // Annual recruitment deviation; n = [nspp, nyrs]
  PARAMETER_MATRIX( init_dev );                   // Log initial numbers-at-age under initMode 0, log deviations off equilibrium under 2/3/4/6; n = [nspp, max_age]

  // -- 3.2. Natural mortality (M1)
  PARAMETER_ARRAY( log_M1 );                       // Natural mortality: residual under multispecies mode, total under single species; n = [nspp, nsex, nages]
  PARAMETER_ARRAY( log_M1_dev );                   // Natural mortality annual deviate; n = [nspp, nsex, nyrs]
  PARAMETER_ARRAY( M1_beta );                     // Regression coefficients for environmnetally linked M1; n = [nspp, nsex, n indices]
  PARAMETER_ARRAY( M1_rho );                      // Correlation for AR1 random effects on age and year; n = [nspp, nsex, 2]
  PARAMETER_ARRAY( M1_dev_log_sd );                // Innovation sd of the M1 random effects, one per species/sex; n = [nspp, nsex]

  // -- 3.3. Growth
  PARAMETER_ARRAY(log_growth_pars);                // Mean growth curve parameters [sp, sex, par]
  // Time-varying growth comes from growth_linkage_offset below.
  PARAMETER_ARRAY(growth_log_sd);                  // Log standard deviation of length-at- min and max age [sp, sex, 2]
  PARAMETER_MATRIX(weight_length_pars);           // Length-weight parameters [sp, (alpha, beta)]

  // -- 3.3b. Linkage-table coefficients (see section 2.3b). Aligned row-for-row with the
  //          linkage_* DATA vectors above. May be length 0 (no linkages supplied).
  PARAMETER_VECTOR(beta_linkage);

  // Random-effect linkage coefficients and their hyperparameters. All length 0 unless a
  // `~ (1|group)` / rw() / ar1() linkage was supplied; when present, beta_linkage_re
  // enters the Laplace approximation (added to random=) and its density goes to
  // JNLL_LINKAGE_RE (13.14).
  PARAMETER_VECTOR(beta_linkage_re);     // RE deviation coefficients (Laplace-integrated)
  // Penalized RE deviations: same density, outside `random` (TMB picks random effects by
  // name). Fixed SD (and AR1 rho) only; see the vignette, "Integrated or penalized".
  PARAMETER_VECTOR(beta_linkage_re_pen);
  PARAMETER_VECTOR(log_sigma_linkage);   // one log-SD per RE group
  PARAMETER_VECTOR(trans_rho_linkage);   // one per AR1 group; rho = tanh(trans_rho_linkage), via rho_trans()
  PARAMETER_VECTOR(beta_linkage_obs);    // Rogers QAR1 effect size: one per observed ar1 group
  PARAMETER_VECTOR(log_obs_sd_linkage);  // Rogers QAR1 observation log-SD: one per observed ar1 group (estimated if obs_sd_est = TRUE)

  // -- 3.4. Fishing mortality parameters
  PARAMETER_VECTOR( log_Flimit );                  // Target fishing mortality for projections on log scale; n = [nspp, nyrs]
  PARAMETER_VECTOR( log_Ftarget );                 // Target fishing mortality for projections on log scale; n = [nspp, nyrs]
  PARAMETER_VECTOR( log_Finit );                   // Fishing mortality for initial population to induce non-equilibrium; n = [nspp]
  PARAMETER_VECTOR( proj_F_prop );                // Proportion of fishing mortality from each fleet for projections; n = [n_fsh]
  PARAMETER_MATRIX( log_F );                       // Annual fishing mortality; n = [n_fsh, nyrs]

  // -- 3.4. Survey catchability parameters
  PARAMETER_VECTOR( index_log_q );                 // Survey catchability; n = [n_index]
  PARAMETER_MATRIX( index_q_beta );               // Survey catchability regression coefficients on env_data columns
  PARAMETER_MATRIX( index_q_dev );                // Annual survey catchability deviates; n = [n_index, nyrs_hind]
  PARAMETER_VECTOR( index_q_log_sd );              // Log standard deviation of prior on survey catchability; n = [1, n_index]
  PARAMETER_VECTOR( index_q_dev_log_sd );// Log standard deviation of time varying survey catchability; n = [1, n_index]

  // -- 3.5. Selectivity parameters
  PARAMETER_ARRAY( sel_coff );                    // selectivity parameters for non-parametric; n = [n_selectivities, nsex, n_sel_bins]
  PARAMETER_ARRAY( sel_coff_dev );                // Annual deviates for non-parametric selectivity parameters; n = [n_selectivities, nsex, n_sel_bins]
  PARAMETER_ARRAY( log_sel_slp );                  // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( sel_inf );                     // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( log_sel_slp_dev );              // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_ARRAY( sel_inf_dev );                 // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_ARRAY( log_sel_apical );              // per-sex log multiplier on the whole curve, after the form and before normalization; n = [n_selectivities, nsex]
  PARAMETER_ARRAY( sel_dn6 );                     // DoubleNormalSS3 (SS3 pattern 24) parameters on SS3's scales: peak, logit top, log asc, log desc, logit init, logit final; n = [6, n_selectivities, nsex]
  PARAMETER_VECTOR( sel_dev_log_sd );              // Log standard deviation of selectivity; n = [1, n_selectivities]
  PARAMETER_MATRIX( sel_curve_pen );              // Non-parametric selectivity penalty columns (Sel_curve_pen1-3), read differently per form -- see the penalty blocks in section 13.7

  // -- 3.6. Data variance
  //FIXME: remove log_sd terms below
  PARAMETER_VECTOR( index_log_sd );                // Log standard deviation of survey index time-series; n = [1, n_index]
  PARAMETER_VECTOR( catch_log_sd );                // Log standard deviation of fishery catch time-series; n = [1, n_fsh]
  PARAMETER_VECTOR( comp_weights );               // Weights for composition data
  PARAMETER_VECTOR( caal_weights );               // Weights for CAAL data
  vector<Type>  DM_pars_comp = exp(comp_weights); // Dirichlet-multinomial scalars
  vector<Type>  DM_pars_caal = exp(caal_weights); // Dirichlet-multinomial scalars

  PARAMETER_VECTOR( diet_comp_weights );          // Weights for diet composition data
  vector<Type>  DM_diet_pars = exp(diet_comp_weights);// Dirichlet-multinomial scalars

  // 3.8. Gamma selectivity parameters
  PARAMETER_VECTOR( log_gam_a );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = 1.0e-10 and 19.9
  PARAMETER_VECTOR( log_gam_b );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = -5.2 and 10

  // 3.9. Preference
  PARAMETER_MATRIX( log_phi );                    // Species preference coefficient; n = [nspp, nspp]


  /** ------------------------------------------------------------------------ //
   * 4. DERIVED QUANTITIES SECTION ----                                        //
   * ------------------------------------------------------------------------- */

  // 4.1. Derived indices
  int n_flt = flt_type.size();
  vector<int> joint_adjust(comp_obs.rows()); joint_adjust.setZero();
  Type penalty = 0.0;
  Type ricker_intercept = 0.0;

  // -- 4.2. Growth
  array<Type> growth_matrix(nspp * 2 + n_flt, max_sex, max_age, max_nlengths_pop, nyrs); growth_matrix.setZero(); // Growth transition matrix on the POPULATION length bins, per fleet and per species derived quantity (biomass, ssb)
  array<Type> weight_hat(nspp * 2 + n_flt, max_sex, max_age, nyrs); weight_hat.setZero(); // Estimated weight-at-age, per fleet and per species derived quantity (biomass, ssb)
  array<Type> mat_weight_hat(nspp * 2 + n_flt, max_sex, max_age, nyrs); mat_weight_hat.setZero(); // Mature weight-at-age (kg), maturity-at-length species only
  array<Type> length_hat(nspp * 2 + n_flt, max_sex, max_age, nyrs); length_hat.setZero(); // Estimated length-at-age, per fleet and per species derived quantity (biomass, ssb)

  // -- 4.3. Estimated population quantities
  vector<Type>  pop_scalar = exp(log_pop_scalar);                                   // Multiplier on input numbers-at-age, per species
  vector<Type>  avg_R(nspp); avg_R.setZero();                                       // Mean recruitment of hindcast
  matrix<Type>  R_hat(nspp, nyrs); R_hat.setZero();                                 // Expected recruitment given SR curve
  matrix<Type>  mort_sum(nspp, max_age); mort_sum.setZero();
  // Fishery selectivity at age used by initMode 6 to spread Finit over ages;
  // filled in section 5.12.1, and 1 everywhere for every other mode.
  array<Type>   sel_init(nspp, max_sex, max_age); sel_init.setZero();
  // Initial age structure WITHOUT the initial deviates: the equilibrium they are applied
  // to. Only the equilibrium catch reads it, and that is the quantity it is defined on --
  // SS3 computes the equilibrium catch from the deviation-free equilibrium, not from the
  // year-1 numbers. Filled alongside N_at_age in section 6.5 so the two cannot drift
  // apart.
  array<Type>   N_eq(nspp, max_sex, max_age); N_eq.setZero();
  matrix<Type>  R0(nspp, nyrs); R0.setZero();                                       // Equilibrium recruitment at F = 0.
  matrix<Type>  alpha(nspp, nyrs); alpha.setZero();                                 // Stock recruit alpha
  matrix<Type>  Beta(nspp, nyrs); Beta.setZero();                                   // Stock recruit beta
  matrix<Type>  steepness(nspp, nyrs); steepness.setZero();                         // Expected % of R0 at 20% SSB0.
  vector<Type>  R_init(nspp); R_init.setZero();                                     // Equilibrium recruitment at F = Finit (non-equilibrium).
  matrix<Type>  R(nspp, nyrs); R.setZero();                                         // Estimated recruitment (thousands of fish; n-at-age is in thousands)
  array<Type>   biomass_at_age(nspp, max_sex, max_age, nyrs); biomass_at_age.setZero();// Estimated biomass-at-age (mt; thousands of fish x kg)
  matrix<Type>  biomass(nspp, nyrs); biomass.setZero();                             // Estimated biomass (mt)
  matrix<Type>  exploitable_biomass(nspp, nyrs); exploitable_biomass.setZero();     // Estimated exploitable biomass (mt)
  matrix<Type>  ssb(nspp, nyrs); ssb.setZero();                                     // Estimated spawning stock biomass (mt)
  matrix<Type>  biomass_depletion(nspp, nyrs); biomass_depletion.setZero();         // Estimated biomass biomass_depletion
  matrix<Type>  ssb_depletion(nspp, nyrs); ssb_depletion.setZero();                 // Estimated biomass_depletion of spawning stock biomass
  array<Type>   M_at_age(nspp, max_sex, max_age, nyrs); M_at_age.setZero();         // Total natural mortality at age
  array<Type>   M_at_age_dB0(nspp, max_sex, max_age, nyrs); M_at_age_dB0.setZero(); // Total natural mortality at age (dynamic B0)
  array<Type>   M_at_age_dBF(nspp, max_sex, max_age, nyrs); M_at_age_dBF.setZero(); // Total natural mortality at age (dynamic BF)
  array<Type>   M1_at_age(nspp, max_sex, max_age, nyrs); M1_at_age.setZero();       // Residual or total natural mortality at age
  array<Type>   N_at_age(nspp, max_sex, max_age, nyrs); N_at_age.setZero();         // Numbers at age
  array<Type>   avgN_at_age(nspp, max_sex, max_age, nyrs); avgN_at_age.setZero();   // Average numbers-at-age
  array<Type>   avgN_at_age_dB0(nspp, max_sex, max_age, nyrs); avgN_at_age_dB0.setZero();   // Average numbers-at-age (dynamic B0)
  array<Type>   avgN_at_age_dBF(nspp, max_sex, max_age, nyrs); avgN_at_age_dBF.setZero();   // Average numbers-at-age (dynamic BF)
  array<Type>   Z_at_age(nspp, max_sex, max_age, nyrs); Z_at_age.setZero();         // Total mortality at age
  vector<Type>  R_sd(nspp); R_sd.setZero();                                         // Standard deviation of recruitment variation
  vector<Type>  zero_N_pen(nspp); zero_N_pen.setZero();                             // Additional penalty to add to likelihood if n-at-age goes < 0

  // -- 4.4. Selectivity parameters
  array<Type>   sel_at_age(n_flt, max_sex, max_age, nyrs); sel_at_age.setZero();    // Estimated selectivity at age
  array<Type>   sel_at_length(n_flt, max_sex, max_nlengths_pop, nyrs); sel_at_length.setZero();// Estimated selectivity at length, on the POPULATION length bins
  array<Type>   avg_sel(n_flt, max_sex, nyrs_hind); avg_sel.setZero();              // Average selectivity for non-parametric up to n_sel_bins
  array<Type>   non_par_sel(n_flt, max_sex, max_bin, nyrs); non_par_sel.setZero();  // Estimated selectivity for AMAK non-parametric (pre-normalization)
  array<Type>   log_non_par_sel(n_flt, max_sex, max_bin, nyrs); log_non_par_sel.setZero(); // Same curve on the log scale, where the penalties are written (exp() underflows below about -745, so it is carried, not recovered)
  vector<Type>  sel_dev_sd(n_flt); sel_dev_sd.setZero();                            // Standard deviation of selectivity deviates

  // -- 4.5. Fishery components
  matrix<Type>  F_spp(nspp, nyrs); F_spp.setZero();                                 // Fully selected fishing mortality by species
  matrix<Type>  F_flt(n_flt, nyrs); F_flt.setZero();                                // Fully selected fishing mortality by fleet
  array<Type>   F_flt_age(n_flt, max_sex, max_age, nyrs); F_flt_age.setZero();            // Estimated fishing mortality-at-age/sex for each fishery
  array<Type>   Flimit_at_age(nspp, max_sex, max_age, nyrs); Flimit_at_age.setZero();   // Estimated limit fishing mortality-at-age/sex for each species (FOFL proxy)
  array<Type>   Ftarget_at_age(nspp, max_sex, max_age, nyrs); Ftarget_at_age.setZero(); // Estimated target fishing mortality-at-age/sex for each species (max FABC proxy)
  array<Type>   F_at_age(nspp, max_sex, max_age, nyrs); F_at_age.setZero();       // Sum of annual estimated fishing mortalities for each species-at-age
  vector<Type>  catch_hat(catch_obs.rows()); catch_hat.setZero();                   // Estimated fishery yield/numbers (mt, or thousands of fish)
  vector<Type>  max_catch_hat(catch_obs.rows()); max_catch_hat.setZero();           // Estimated exploitable biomass/numbers by fleet (mt, or thousands of fish)
  // As index_sd below, for the fishery. Catch is lognormal throughout, so this is
  // always a log-scale sd (unitless), not a log of anything. rename_output()
  // also writes it as log_catch_sd, the name fits saved before 5.9.0 carry.
  vector<Type>  catch_sd(catch_obs.rows()); catch_sd.setZero();
  // Concentrated (analytical) catch sd used when Estimate_catch_sd = 2; a
  // log-scale sd (unitless), not a log. Mirrors index_analytical_sd below.
  vector<Type>  catch_analytical_sd(n_flt); catch_analytical_sd.setZero();
  vector<Type>  catch_n_obs(n_flt); catch_n_obs.setZero();  // fitted catch rows per fleet

  // -- 4.6. Biological reference points
  array<Type>   NByage0(nspp, max_sex, max_age, nyrs); NByage0.setZero();                 // Numbers at age at mean recruitment and F = 0
  array<Type>   NByageF(nspp, max_sex, max_age, nyrs); NByageF.setZero();                 // Numbers at age at mean recruitment and F = Ftarget; the age structure behind SBF
  array<Type>   N_at_age_dB0(nspp, max_sex, max_age, nyrs); N_at_age_dB0.setZero();   // Numbers at age at F = 0 (accounts for annual recruitment)
  array<Type>   N_at_age_dBF(nspp, max_sex, max_age, nyrs); N_at_age_dBF.setZero();   // Female numbers at age at F = Ftarget (accounts for annual recruitment)
  matrix<Type>  DynamicSB0(nspp, nyrs); DynamicSB0.setZero();                       // Estimated dynamic spawning biomass at F = 0 (accounts for S-R curve)
  matrix<Type>  DynamicB0(nspp, nyrs); DynamicB0.setZero();                         // Estimated dynamic  biomass at F = 0 (accounts for S-R curve)
  matrix<Type>  DynamicSBF(nspp, nyrs); DynamicSBF.setZero();                       // Estimated dynamic spawning biomass at F = Ftarget (accounts for S-R curve)
  // Zeroed explicitly, like every other array here: the recursion writes only
  // ages below a species' own nages, and nothing at all under predation, so the
  // rest is read from the REPORT as whatever it was initialised to.
  array<Type>   NbyageSPR(4, nspp, max_age); NbyageSPR.setZero();                   // Numbers at age per recruit for the SPR reference points; slots 0-3 are F = 0, Flimit, Ftarget, Finit
  vector<Type>  SPRlimit(nspp); SPRlimit.setZero();                                 // Estimated Plimit SPR
  vector<Type>  SPRtarget(nspp); SPRtarget.setZero();                               // Estimated Ptarget SPR
  vector<Type>  SPR0(nspp); SPR0.setZero();                                         // Estimated spawning biomass per recruit at F = 0
  vector<Type>  SPRFinit(nspp); SPRFinit.setZero();                                 // Estimated spawning biomass per recruit at Finit
  matrix<Type>  SB0(nspp, nyrs); SB0.setZero();                                     // Estimated spawning stock biomass at F = 0 (Accounts for S-R)
  matrix<Type>  SBF(nspp, nyrs); SBF.setZero();                                     // Estimated spawning stock biomass at F = target (Accounts for S-R)
  matrix<Type>  B0(nspp, nyrs); B0.setZero();                                       // Estimated biomass at F = 0 (Accounts for S-R)
  vector<Type>  Flimit = exp(log_Flimit);                                            // Limit F (FOFL proxy, e.g. F35%) on natural scale
  vector<Type>  Ftarget = exp(log_Ftarget);                                          // Target F (max FABC proxy, e.g. F40%) on natural scale
  vector<Type>  Finit = exp(log_Finit);                                              // Initial F for non-equilibrium age-structure
  vector<Type>  index_q_mult(index_q_beta.cols()); index_q_mult.setZero();          // Environmental design matrix for q
  vector<Type>  M1_mult(M1_beta.dim(2)); M1_mult.setZero();                         // Environmental design matrix for M1
  vector<Type>  beta_q_tmp(index_q_beta.cols()); beta_q_tmp.setZero();              // Q beta parameters by species, for the matrix multiply
  vector<Type>  env_q_tmp(index_q_beta.cols()); env_q_tmp.setZero();                // Q env data by year, for the matrix multiply
  vector<Type>  beta_M1_tmp(M1_beta.dim(2)); beta_M1_tmp.setZero();                 // Beta parameters by species, for the matrix multiply
  vector<Type>  env_M1_tmp(M1_beta.dim(2)); env_M1_tmp.setZero();                   // Env data by year, for the matrix multiply
  matrix<Type>  proj_F(nspp, nyrs); proj_F.setZero();                               // Projected F (Fabc/Ftac/etc) using harvest control rule

  // -- 4.7. Survey components
  vector<Type>  index_q_sd(n_flt); index_q_sd.setZero();                            // Vector of standard deviation of survey catchability prior
  vector<Type>  index_q_dev_sd(n_flt); index_q_dev_sd.setZero();                    // Vector of standard deviation of time-varying survey catchability deviation
  vector<Type>  index_hat(index_obs.rows()); index_hat.setZero();                   // Estimated survey index; mt or thousands of fish per the fleet's Observation_units
  // Rejection bookkeeping for the natural-scale survey draws (sim_mod only), per
  // index_obs row: a non-zero tries count marks the rows a rejection-capable branch drew
  // -- the untruncated "Normal" family and MVN/MVNORM, which must redraw a non-positive
  // index. "TruncatedNormal" draws by inverse CDF and never enters here. sim_mod() sizes
  // the draw/density gap from index_hat and index_sd, not from these counts.
  vector<Type>  index_trunc_tries_sim(index_obs.rows());   index_trunc_tries_sim.setZero();
  vector<Type>  index_trunc_rejects_sim(index_obs.rows()); index_trunc_rejects_sim.setZero();
  // The marginal sd the draw itself used, so sim_mod() can size the truncated mass as
  // Phi(-index_hat/sd). Separate from index_sd because a covariance fleet draws from
  // index_cov_mat, leaving its index_sd whatever the (unused) Index_sd column carried.
  vector<Type>  index_trunc_sd_sim(index_obs.rows());      index_trunc_sd_sim.setZero();
  // The redraw budget ACTUALLY applied to each row, which is what sim_mod() quotes: the
  // constant below on the univariate branch, and that scaled by the fleet's row count on
  // the correlated one, which rejects the whole vector at once, so its joint rejection
  // probability climbs with the row count.
  vector<Type>  index_trunc_budget_sim(index_obs.rows()); index_trunc_budget_sim.setZero();
  const int index_trunc_max_tries = 100;
  // The observation sd the index likelihood actually used for each row, whichever of the
  // three est_sigma_index routes supplied it. NOT a log: the sd itself, on the scale that
  // fleet's Index_distribution works on -- a log-scale sd (unitless) for Lognormal, an
  // ABSOLUTE sd in the units of the index for the natural-scale families. rename_output()
  // also writes it as log_index_sd, the name fits saved before 5.9.0 carry.
  vector<Type>  index_sd(index_obs.rows()); index_sd.setZero();
  vector<Type>  index_analytical_sd(n_flt); index_analytical_sd.setZero();          // Concentrated (analytical) index sd used when Estimate_index_sd = 2; a log-scale sd, not a log
  vector<Type>  index_q_analytical(n_flt); index_q_analytical.setZero();            // Analytical catchability (not an sd) used when Catchability = 3; "AnalyticalArith" (7) has its own block at 8.2b
  matrix<Type>  index_q(n_flt, nyrs_hind); index_q.setZero();                       // Estimated survey catchability //FIXME: extend out to full time-series
  vector<Type>  index_n_obs(n_flt); index_n_obs.setZero();                          // Observations per survey time series

  // -- 4.8. Composition data
  // age_hat / age_obs_hat are indexed by AGE whatever dimension the fleet's composition
  // data are on -- a length composition is built at age, then converted through the
  // age-length transition -- so the width is the widest age index any row writes:
  // nages(sp), or nages(sp) * 2 for a joint-sex row (comp_ctl column 2 == 3), never
  // narrower than comp_obs. Sizing from comp_obs alone gives only the workbook's Comp_
  // columns, i.e. nlengths on a length-only model, and Eigen does not bounds-check in a
  // release build, so nlengths < nages would write silently into adjacent memory. Species
  // is range-checked because nothing between rearrange_data() and MakeADFun validates
  // comp_data$Species, and a non-numeric one arrives as NA_integer_.
  int max_age_cols = static_cast<int>(comp_obs.cols());
  for(int comp_row = 0; comp_row < comp_ctl.rows(); comp_row++){
    int comp_row_sp = comp_ctl(comp_row, 1) - 1;
    if((comp_row_sp < 0) | (comp_row_sp >= nspp)){
      error("comp_data Species on row %d is not in 1..nspp", comp_row + 1);
    }
    int comp_row_cols = (comp_ctl(comp_row, 2) == 3) ? nages(comp_row_sp) * 2 : nages(comp_row_sp);
    if(comp_row_cols > max_age_cols){ max_age_cols = comp_row_cols; }
  }
  vector<Type>  n_hat(comp_obs.rows()) ; n_hat.setZero() ;                          // Estimated catch (numbers)
  matrix<Type>  age_hat(comp_obs.rows(), max_age_cols); age_hat.setZero();          // Estimated catch at true age
  matrix<Type>  age_obs_hat(comp_obs.rows(), max_age_cols); age_obs_hat.setZero();  // Estimated catch at observed age (accounts for ageing error)
  matrix<Type>  comp_hat = comp_obs; comp_hat.setZero();                            // Estimated comp
  matrix<Type>  caal_hat = caal_obs; caal_hat.setZero();                            // Estimated CAAL
  array<Type>   pred_CAAL(n_flt, max_sex, max_age, max_nlengths, nyrs); pred_CAAL.setZero(); // Predicted CAAL for each fleet
  array<Type>   comp_pred(max_sex, max_age, max_nlengths); comp_pred.setZero();      // Predicted fleet numbers by age and length for one comp row

  // -- 4.9. Ration components
  array<Type>   consumption_at_age( nspp, max_sex, max_age, nyrs ); consumption_at_age.setZero(); // Annual ration at age (kg/yr)
  matrix<Type>  fT( nspp, nyrs ); fT.setZero();                                     // Pre-allocation of temperature function of consumption

  // -- 4.10. Diet components
  array<Type>   diet_prop(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); diet_prop.setZero();             // Stomach proportion by weight U
  array<Type>   diet_prop_hat(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); diet_prop_hat.setZero();     // Predicted stomach proportion by weight U
  array<Type>   other_food_diet_prop(nspp, max_sex, max_age, nyrs); other_food_diet_prop.setZero();                 // Other food diet proportion by weight
  matrix<Type>  diet_hat = diet_obs; diet_hat.setZero();                                                            // Estimated stomach proportion by weight U (formated following data input)

  // -- 4.11. Suitability components
  array<Type>   avail_food(nspp, max_sex, max_age, nyrs); avail_food.setZero();                               // Available food to predator
  array<Type>   avail_food_dB0(nspp, max_sex, max_age, nyrs); avail_food_dB0.setZero();                       // Available food to predator (dynamic B0)
  array<Type>   avail_food_dBF(nspp, max_sex, max_age, nyrs); avail_food_dBF.setZero();                       // Available food to predator (dynamic BF)
  array<Type>   stom_div_bio(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); stom_div_bio.setZero(); // Stomach proportion over biomass; U/ (W * N)
  array<Type>   suitability(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); suitability.setZero();   // Suitability/gamma selectivity of predator age u on prey age a
  array<Type>   suit_other(nspp, max_sex, max_age, nyrs); suit_other.setZero();                               // Suitability not accounted for by the included prey
  array<Type>   other_diet_prop_hat(nspp, max_sex, max_age, nyrs); other_diet_prop_hat.setZero();             // Diet of prey not included in the model

  // -- 4.12. Suitability parameters
  vector<Type> gam_a = exp(log_gam_a);                                    // Predator size-selectivity: shape parameter for gamma suitability, mean for normal of logs
  vector<Type> gam_b = exp(log_gam_b);                                    // Predator size-selectivity: scale parameter for gamma suitability, sd for normal of logs
  matrix<Type> vulnerability(nspp, nspp); vulnerability.setZero();        // Predator-prey preference coefficients
  vector<Type> vulnerability_other(nspp); vulnerability_other.setZero();  // Preference for other food

  // -- 4.13. Predation components
  array<Type>   M2_at_age(nspp, max_sex, max_age, nyrs); M2_at_age.setZero();                        // Predation mortality at age
  array<Type>   M2_at_age_dB0(nspp, max_sex, max_age, nyrs); M2_at_age_dB0.setZero();                // Predation mortality at age (dynamic B0)
  array<Type>   M2_at_age_dBF(nspp, max_sex, max_age, nyrs); M2_at_age_dBF.setZero();                // Predation mortality at age (dynamic BF)
  array<Type>   M2_prop(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); M2_prop.setZero();  // Each predator-at-age's contribution to prey M2 (yr^-1), not a share of it
  array<Type>   B_eaten(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); B_eaten.setZero();  // Biomass of prey eaten via predation by a predator at age
  array<Type>   B_eaten_as_prey(nspp, max_sex, max_age, nyrs); B_eaten_as_prey.setZero();            // Biomass eaten as prey via predation

  /** ------------------------------------------------------------------------ //
   * 5. INITIAL CALCULATIONS                                                   //
   * ------------------------------------------------------------------------- */

  // 5.1. DATA VARIANCE TERMS
  R_sd = exp(R_log_sd); // Convert log sd to natural scale
  sel_dev_sd = exp(sel_dev_log_sd) ;
  index_q_sd = exp(index_q_log_sd) ;
  index_q_dev_sd = exp(index_q_dev_log_sd) ;
  Cindex -=1; // Subtract 1 from Cindex to deal with indexing start at 0


  // -- 5.2. SIMULATE LINKAGE RANDOM EFFECTS (sim_mod(simulate = TRUE))
  // A random linkage on any of the five processes draws here, before they read it, from
  // its density in rceattle_re_group_nll() (linkage.hpp, called in 13.14). Slots are in
  // time order within a group, so a random walk starts from the fitted first slot, and an
  // AR1 draw is standardized then scaled by the MARGINAL sd sigma.
  //
  // beta_linkage_re_drawn_sim marks which slots this call wrote. The reported deviation
  // vector spans EVERY random linkage, drawn or not, so without the mask a caller would
  // score recovery on slots left at their fitted values.
  SIMULATE {
    vector<Type> beta_linkage_re_drawn_sim(linkage_re_sigma.size());
    beta_linkage_re_drawn_sim.setZero();
    if (log_sigma_linkage.size() > 0 && simulate_period(0) == 1) {
      int n_grp  = log_sigma_linkage.size();

      // Drawn only if EVERY row using the group names a process the caller asked for.
      // A group belongs to one process (.assign_re_registry() keys it on process,
      // parameter, stratum and random term), so this reads that process's simulate_state.
      vector<int> grp_used(n_grp), grp_draw(n_grp);
      grp_used.setZero();
      grp_draw.fill(1);
      for (int i = 0; i < linkage_re_index.size(); ++i) {
        if (linkage_re_index(i) < 0) continue;
        int grp = linkage_re_sigma(linkage_re_index(i));
        grp_used(grp) = 1;
        // linkage_process also codes composition (5), which simulate_state does not cover
        // because composition linkages carry priors, not random effects. Treat anything
        // outside the five as "not asked for" rather than reading past simulate_state.
        int proc = linkage_process(i);
        if (proc < 0 || proc >= simulate_state.size() || simulate_state(proc) != 1) {
          grp_draw(grp) = 0;
        }
      }

      for (int grp = 0; grp < n_grp; ++grp) {
        if (grp_used(grp) == 0 || grp_draw(grp) == 0) continue;

        // Rogers QAR1: the latent is measured by an observed covariate series, so
        // redrawing it alone would leave the two describing different histories.
        // Left as fitted; sim_mod() warns that it was.
        if (linkage_re_obs(grp) >= 0) continue;

        vector<int> slots = rceattle_group_slots(linkage_re_sigma, grp);
        int len = slots.size();
        if (len == 0) continue;

        Type sigma = exp(log_sigma_linkage(grp));
        vector<Type> re(len);
        int st = linkage_re_struct(grp);

        // Each slot sits in beta_linkage_re if integrated, else beta_linkage_re_pen, at
        // position linkage_re_slot(slot).
        if (st == RCEATTLE_RE_RW) {
          re(0) = linkage_re_integrate(slots(0))
              ? beta_linkage_re(linkage_re_slot(slots(0)))
              : beta_linkage_re_pen(linkage_re_slot(slots(0)));
          for (int t = 1; t < len; ++t) re(t) = re(t - 1) + rnorm(Type(0), sigma);
        } else if (st == RCEATTLE_RE_AR1) {
          Type rho = rho_trans(trans_rho_linkage(linkage_re_rho(grp)));
          AR1(rho).simulate(re);
          re *= sigma;
        } else {                                         // RCEATTLE_RE_IID
          for (int t = 0; t < len; ++t) re(t) = rnorm(Type(0), sigma);
        }

        for (int t = 0; t < len; ++t) {
          if (linkage_re_integrate(slots(t))) {
            beta_linkage_re(linkage_re_slot(slots(t))) = re(t);
          } else {
            beta_linkage_re_pen(linkage_re_slot(slots(t))) = re(t);
          }
          // A random walk keeps its first slot at the fitted level (it has no density
          // term of its own), so that slot is written but not drawn.
          if (!(st == RCEATTLE_RE_RW && t == 0)) beta_linkage_re_drawn_sim(slots(t)) = 1;
        }
      }
    }
    REPORT(beta_linkage_re_drawn_sim);
  }

  // -- 5.2b. EFFECTIVE LINKAGE COEFFICIENTS (beta_linkage_re_all, beta_linkage_eff)
  // Every deviation in one vector, indexed by RE slot: integrated deviations in
  // beta_linkage_re, penalized ones in beta_linkage_re_pen, with linkage_re_slot giving
  // the position within whichever holds it. With no penalized group it copies
  // beta_linkage_re.
  vector<Type> beta_linkage_re_all(linkage_re_sigma.size());
  for (int s = 0; s < beta_linkage_re_all.size(); ++s) {
    beta_linkage_re_all(s) = linkage_re_integrate(s)
        ? beta_linkage_re(linkage_re_slot(s))
        : beta_linkage_re_pen(linkage_re_slot(s));
  }

  // THE _sim NAMING RULE, stated once here and referred to below. TMB never clears the
  // report environment, so anything REPORTed inside a SIMULATE block stays visible in
  // obj$report() for this object's life -- under the observed object's own name, a random
  // replicate would read as the data. Every draw is reported under a name ending _sim,
  // and the composition, CAAL and diet draws go into separate copies (comp_sim, caal_sim,
  // diet_sim) for the same reason.
  SIMULATE {
    vector<Type> beta_linkage_re_sim = beta_linkage_re_all;
    REPORT(beta_linkage_re_sim);
  }

  // Effective linkage coefficient per table row. A fixed row uses its beta_linkage(i); a
  // random-effect row (linkage_re_index >= 0) takes its deviation from the slot-space
  // vector, whose density is JNLL_LINKAGE_RE (13.14). The accumulators below see one beta vector
  // either way; with no RE rows this copies beta_linkage.
  vector<Type> beta_linkage_eff = beta_linkage;
  for (int i = 0; i < beta_linkage_eff.size(); ++i) {
    if (linkage_re_index(i) >= 0) {
      Type z = beta_linkage_re_all(linkage_re_index(i));
      // Rogers QAR1: an observed AR1 latent enters the target scaled by an estimated
      // effect size, beta_linkage_obs. Unobserved groups keep the deviate as is.
      int grp = linkage_re_sigma(linkage_re_index(i));
      if (linkage_re_obs(grp) >= 0) z *= beta_linkage_obs(linkage_re_obs(grp));
      beta_linkage_eff(i) = z;
    }
  }

  // The table every rceattle_apply_*_linkages() call below reads (linkage.hpp): each row
  // adds offset = beta * X(yr, X_col) to one parameter in one stratum (2.3b).
  LinkageTable<Type> linkage_table = {
    linkage_process, linkage_param,                                // which parameter
    linkage_species, linkage_sex, linkage_age_bin, linkage_fleet,  // which stratum (0 = all)
    linkage_X_col, linkage_link, linkage_X, beta_linkage_eff};     // the offset


  // 5.3. CATCHABILITY
  // Environmental linkage offsets, one pass per link scale. Zero unless a
  // q linkage was supplied, so models without one are unaffected.
  matrix<Type> q_linkage_offset(n_flt, nyrs_hind);     q_linkage_offset.setZero();
  matrix<Type> q_linkage_offset_nat(n_flt, nyrs_hind); q_linkage_offset_nat.setZero();
  matrix<Type> q_linkage_log_mult(n_flt, nyrs_hind);   q_linkage_log_mult.setZero();

  // Log-link rows into the log-scale offset, identity-link rows into the natural-scale
  // one, exponential rows into the log of the multiplier on log q (applied as exp()).
  rceattle_apply_q_linkages(q_linkage_offset,     RCEATTLE_LINK_LOG,         linkage_table, n_flt, nyrs_hind);
  rceattle_apply_q_linkages(q_linkage_offset_nat, RCEATTLE_LINK_IDENTITY,    linkage_table, n_flt, nyrs_hind);
  rceattle_apply_q_linkages(q_linkage_log_mult,   RCEATTLE_LINK_EXPONENTIAL, linkage_table, n_flt, nyrs_hind);

  REPORT(q_linkage_offset);
  REPORT(q_linkage_offset_nat);
  REPORT(q_linkage_log_mult);

  for(flt = 0; flt < n_flt; flt++){
    for(yr = 0; yr < nyrs_hind; yr++){
      // q_y = exp((log q + S_log) * exp(S_exp) + index_q_dev) + S_id. The deviate is added
      // after the product, as SS3 applies its devs after the env link; inside, its sd
      // would vary by year.
      index_q(flt, yr) = exp((index_log_q(flt) + q_linkage_offset(flt, yr))
                               * exp(q_linkage_log_mult(flt, yr))
                             + index_q_dev(flt, yr))
                           + q_linkage_offset_nat(flt, yr);

      // Catchability 5 ("Environmental", legacy): log q_y = log q + beta * env_index,
      // replacing the line above; R refuses a q linkage on such a fleet.
      if(est_index_q(flt) == 5){
        beta_q_tmp = index_q_beta.row(flt);
        env_q_tmp = env_index.row(yr) ;
        index_q_mult =  env_q_tmp * beta_q_tmp;
        index_q(flt, yr) = exp(index_log_q(flt) + (index_q_mult).sum());
      }
    }
  }

  matrix<int> r_sexes(diet_obs.rows(), 2); r_sexes.setZero();
  matrix<int> k_sexes(diet_obs.rows(), 2); k_sexes.setZero();


  // 5.4. MATURITY, Finit, AND SEX RATIO
  // -- for SSB derivation, not SPR
  matrix<Type> mature_females = maturity;
  for( sp = 0; sp < nspp ; sp++) {

    if(initMode < 3 || initMode == 5){
      Finit(sp) = 0; // An equilibrium start takes Finit = 0, so R_init and R0 agree. initMode 5 (OffsetEquilibrium) is an F = 0 equilibrium seeded by first-year recruitment.
    }

    // Sex ratio for SSB derivation
    for( age = 0 ; age < nages(sp); age++ ) {
      if(nsex(sp) == 1){
        mature_females( sp, age ) = maturity( sp, age ) * sex_ratio(sp, age); // Multiply sex_ratio and maturity for 1 sex models
      }
    }
  }

  // 5.5. OFFSETS
  // One array per link scale per process: `*_linkage_offset` from log-link rows (S_log),
  // `*_linkage_offset_nat` from identity-link rows (S_id), combined as in section 2.3b.
  // Catchability has a third, `q_linkage_log_mult` (S_exp, 5.3). All zero with no linkages.
  //
  // - RECRUITMENT OFFSETS
  array<Type> recruitment_linkage_offset(nspp,
                                         int(RCEATTLE_N_REC_PARAMS),
                                         nyrs);
  array<Type> recruitment_linkage_offset_nat(nspp,
                                             int(RCEATTLE_N_REC_PARAMS),
                                             nyrs);
  recruitment_linkage_offset.setZero();
  recruitment_linkage_offset_nat.setZero();
  rceattle_apply_recruitment_linkages(recruitment_linkage_offset,     RCEATTLE_LINK_LOG,      linkage_table, nspp, nyrs);
  rceattle_apply_recruitment_linkages(recruitment_linkage_offset_nat, RCEATTLE_LINK_IDENTITY, linkage_table, nspp, nyrs);
  REPORT(recruitment_linkage_offset);
  REPORT(recruitment_linkage_offset_nat);

  // - M OFFSETS
  array<Type> M_linkage_offset(nspp, max_sex, max_age, nyrs);
  array<Type> M_linkage_offset_nat(nspp, max_sex, max_age, nyrs);
  M_linkage_offset.setZero();
  M_linkage_offset_nat.setZero();
  rceattle_apply_M_linkages(M_linkage_offset,     RCEATTLE_LINK_LOG,      linkage_table, nspp, nsex, nages, nyrs);
  rceattle_apply_M_linkages(M_linkage_offset_nat, RCEATTLE_LINK_IDENTITY, linkage_table, nspp, nsex, nages, nyrs);
  REPORT(M_linkage_offset);
  REPORT(M_linkage_offset_nat);


  // - GROWTH OFFSETS
  array<Type> growth_linkage_offset(nspp, max_sex, nyrs, RCEATTLE_N_GROWTH_PARAMS);
  array<Type> growth_linkage_offset_nat(nspp, max_sex, nyrs, RCEATTLE_N_GROWTH_PARAMS);
  growth_linkage_offset.setZero();
  growth_linkage_offset_nat.setZero();
  rceattle_apply_growth_linkages(growth_linkage_offset,     RCEATTLE_LINK_LOG,      linkage_table, nspp, nsex, nyrs);
  rceattle_apply_growth_linkages(growth_linkage_offset_nat, RCEATTLE_LINK_IDENTITY, linkage_table, nspp, nsex, nyrs);
  REPORT(growth_linkage_offset);
  REPORT(growth_linkage_offset_nat);


  // 5.6. RECRUITMENT PARAMETERS
  // R0, alpha, Beta: p_y = exp(log p + S_log) + S_id (2.3b).
  // TODO: a negative nat_offset can make R0, alpha or Beta non-positive; only hindcast R, R_hat and the penalty curve are floored (issue #234).
  for(sp = 0; sp < nspp; sp++){
    for(yr = 0; yr < nyrs; yr++){
      R0(sp, yr)    = exp(rec_pars(sp, RCEATTLE_REC_R0)    + recruitment_linkage_offset(sp, RCEATTLE_REC_R0,    yr))
                    + recruitment_linkage_offset_nat(sp, RCEATTLE_REC_R0,    yr);
      alpha(sp, yr) = exp(rec_pars(sp, RCEATTLE_REC_ALPHA) + recruitment_linkage_offset(sp, RCEATTLE_REC_ALPHA, yr))
                    + recruitment_linkage_offset_nat(sp, RCEATTLE_REC_ALPHA, yr);
      Beta(sp, yr)  = exp(rec_pars(sp, RCEATTLE_REC_BETA)  + recruitment_linkage_offset(sp, RCEATTLE_REC_BETA,  yr))
                    + recruitment_linkage_offset_nat(sp, RCEATTLE_REC_BETA,  yr);
    }
  }


  // 5.7. GROWTH
  // -- Rearrange growth parameters. Log-link offset enters additively
  //    on the log scale (multiplicative on natural K/L1/Linf/m);
  //    identity-link offset enters additively on the natural scale.
  array<Type> growth_parameters(nspp, max_sex, nyrs, RCEATTLE_N_GROWTH_PARAMS); growth_parameters.setZero(); // K, L1, Linf, m
  for(sp = 0; sp < nspp; sp++){
    for(sex = 0; sex < nsex(sp); sex ++){
      for(yr = 0; yr < nyrs; yr++){
        for(int par = 0; par < RCEATTLE_N_GROWTH_PARAMS; par++){
          growth_parameters(sp, sex, yr, par) = exp(
            log_growth_pars(sp, sex, par)
          + growth_linkage_offset(sp, sex, yr, par)
          ) + growth_linkage_offset_nat(sp, sex, yr, par);
        }
      }
    }
  }
  REPORT(growth_parameters);

  // -- Calculate weight
  calculate_weight(
    weight_hat,
    length_hat,
    growth_matrix,
    mat_weight_hat,
    weight_obs,
    growth_model,
    growth_sd_style,
    growth_sd_form,
    growth_plus_length,
    plus_group_decay,
    nspp,
    nyrs,
    nyrs_hind,
    n_flt,
    flt_spp,
    flt_month,
    nsex,
    minage,
    growth_age_L1,
    nages,
    nlengths,
    pop_wt_index,
    ssb_wt_index,
    flt_wt_index,
    spawn_month,
    lengths,
    nlengths_pop,
    lengths_pop,
    growth_parameters,
    growth_log_sd,
    weight_length_pars,
    mat_len_use,
    mat_len_pars,
    log_M1
  );

  // -- Spawning output per fish (kg) at spawning time, female-only for a
  //    two-sex species. Weight x maturity-at-age x sex ratio by default; with
  //    maturity-at-length it is the mature weight integrated over the length
  //    distribution at spawning, so maturity and weight vary together with length.
  array<Type> spawn_output(nspp, max_age, nyrs); spawn_output.setZero();
  for(sp = 0; sp < nspp; sp++){
    wt_idx_ssb = 2 * sp + 1;
    for(age = 0; age < nages(sp); age++){
      for(yr = 0; yr < nyrs; yr++){
        if(mat_len_use(sp) == 1){
          spawn_output(sp, age, yr) = mat_weight_hat(wt_idx_ssb, 0, age, yr) * ((nsex(sp) == 1) ? sex_ratio(sp, age) : Type(1.0));
        } else {
          spawn_output(sp, age, yr) = weight_hat(wt_idx_ssb, 0, age, yr) * mature_females(sp, age);
        }
      }
    }
  }
  REPORT(spawn_output);


  // 5.8. SELECTIVITY
  // Environmental-linkage offsets. Two arrays per parameter family (log and
  // natural scale), all zero unless a selectivity linkage was supplied.
  array<Type> sel_slp_off (2, n_flt, max_sex, nyrs);            sel_slp_off.setZero();
  array<Type> sel_slp_off_nat (2, n_flt, max_sex, nyrs);        sel_slp_off_nat.setZero();
  array<Type> sel_inf_off (2, n_flt, max_sex, nyrs);            sel_inf_off.setZero();
  array<Type> sel_inf_off_nat (2, n_flt, max_sex, nyrs);        sel_inf_off_nat.setZero();
  array<Type> sel_coff_off (n_flt, max_sex, max_bin, nyrs);     sel_coff_off.setZero();
  array<Type> sel_coff_off_nat (n_flt, max_sex, max_bin, nyrs); sel_coff_off_nat.setZero();
  array<Type> sel_apical_off (n_flt, max_sex, nyrs);            sel_apical_off.setZero();
  array<Type> sel_apical_off_nat (n_flt, max_sex, nyrs);        sel_apical_off_nat.setZero();
  array<Type> sel_dn6_off        (6, n_flt, max_sex, nyrs);     sel_dn6_off.setZero();
  array<Type> sel_dn6_off_nat    (6, n_flt, max_sex, nyrs);     sel_dn6_off_nat.setZero();

  // The per-sex apical multiplier is compiled in only when a selectivity
  // linkage names it (RCEATTLE_SEL_APICAL); otherwise the curve is untouched and the AD
  // tape of every existing model is unchanged.
  int sel_apical_on = 0;
  for (int i = 0; i < linkage_process.size(); i++) {
    if (linkage_process(i) == RCEATTLE_PROC_SEL && linkage_param(i) == RCEATTLE_SEL_APICAL) sel_apical_on = 1;
  }

  rceattle_apply_sel_linkages(   // log link, into the log-scale arrays
    sel_slp_off, sel_inf_off, sel_coff_off, sel_apical_off, sel_dn6_off,
    RCEATTLE_LINK_LOG, linkage_table, n_flt, max_sex, max_bin, nyrs);
  rceattle_apply_sel_linkages(   // identity link, into the natural-scale arrays
    sel_slp_off_nat, sel_inf_off_nat, sel_coff_off_nat, sel_apical_off_nat, sel_dn6_off_nat,
    RCEATTLE_LINK_IDENTITY, linkage_table, n_flt, max_sex, max_bin, nyrs);

  // "sel_at_age" and "sel_at_length" modified via pass-by-reference
  // when "sel_at_length" is used, it is converted to "sel_at_age" using the growth matrix
  calculate_selectivity(
    nspp,                 // Number of species.
    n_flt,                // Number of fleets.
    nyrs,                 // Years, hindcast plus projection.
    nyrs_hind,            // Years in the hindcast.
    styr,                 // Start year
    nsex,                 // Sexes modelled per species: 1 = combined, 2 = both.
    nages,                // Age bins per species.
    nlengths,             // Observation length bins per species.
    lengths,              // Length bin boundaries matrix
    nlengths_pop,         // Population length bins per species (= nlengths when no pop grid)
    lengths_pop,          // Population length bin lower edges
    flt_spp,              // Fleet to species mapping
    flt_sel_type,         // Selectivity model type per fleet
    flt_varying_sel,      // Time_varying_sel per fleet (picks NonParametricIntegrable's construction)
    flt_sel_dim,          // Age or length based
    bin_first_selected,   // Min bin selected per fleet
    flt_n_sel_bins,       // Max estimated bins per fleet
    flt_sel_cap_bin,      // Bin (0-based) at/after which realized non-par sel is capped flat (NonParametricPM)
    sel_norm_bin1,        // Normalization control/bin 1
    sel_norm_bin2,        // Normalization control/bin 2
    sel_norm_scope,       // Normalization reference pooled across sexes?
    flt_sel_start_yr,     // Per-fleet selectivity start year (0-based)
    emp_sel_obs,          // Empirical observations matrix
    emp_sel_ctl,          // Empirical control matrix
    log_sel_slp,           // Logistic slope parameters
    log_sel_slp_dev,       // Slope deviations
    sel_inf,              // Inflection parameters
    sel_inf_dev,          // Inflection deviations
    sel_coff,             // Non-parametric coefficients
    sel_coff_dev,         // Coefficient deviations
    avg_sel,              // [Modified] Average selectivity
    non_par_sel,          // [Modified] Unnormalized non-parametric selectivity
    log_non_par_sel,      // [Modified] The same curve on the log scale, for the penalties
    sel_at_length,        // [Modified] Final length-based selectivity
    sel_at_age,           // [Modified] Final age-based selectivity
    growth_matrix,        // Length to age transition matrix
    sel_slp_off, sel_slp_off_nat,   // selectivity linkage offsets (log / natural)
    sel_inf_off, sel_inf_off_nat,
    sel_coff_off, sel_coff_off_nat,
    sel_apical_on, log_sel_apical,  // per-sex apical height, only when a linkage names it
    sel_apical_off, sel_apical_off_nat,
    sel_dn6, sel_dn6_off, sel_dn6_off_nat, sel_dn6_ends   // DoubleNormalSS3
  );

  // -- Selected body weight (kg): sum_l P(l|a) s(l) w(l) / sum_l P(l|a) s(l),
  //    SS3's "bodywt". A length-selective fleet takes the larger or smaller fish
  //    of an age class, so its catch and survey biomass weigh that class by the
  //    fish it selects, not by the class mean. Estimated growth only.
  for(flt = 0; flt < n_flt; flt++){
    sp = flt_spp(flt);
    if(growth_model(sp) == 0 || flt_sel_dim(flt) != 1 || flt_sel_type(flt) == 0) continue;
    int wt_idx_sel = nspp * 2 + flt;
    for(sex = 0; sex < nsex(sp); sex++){
      for(age = 0; age < nages(sp); age++){
        for(yr = 0; yr < nyrs; yr++){
          Type num = 0, den = 0;
          // Integrated over POPULATION bins: both the key and the curve live
          // there, so the selectivity weighting is applied at the resolution
          // the curve is defined on rather than to a bin average.
          for(int ln = 0; ln < nlengths_pop(sp); ln++){
            // The last bin takes the previous bin's width; a single-bin grid has
            // no neighbour to take one from, so it keeps the edge itself.
            Type lenmid = (ln < nlengths_pop(sp) - 1)
              ? (lengths_pop(sp, ln) + lengths_pop(sp, ln + 1)) / Type(2.0)
              : (nlengths_pop(sp) < 2
                   ? lengths_pop(sp, ln)
                   : lengths_pop(sp, ln) + (lengths_pop(sp, ln) - lengths_pop(sp, ln - 1)) / Type(2.0));
            Type p_sel = growth_matrix(wt_idx_sel, sex, age, ln, yr) * sel_at_length(flt, sex, ln, yr);
            num += p_sel * weight_length_pars(sp, 0) * pow(lenmid, weight_length_pars(sp, 1));
            den += p_sel;
          }
          // An age class the fleet does not select contributes nothing either way.
          weight_hat(wt_idx_sel, sex, age, yr) = num / (den + Type(1e-30));
        }
      }
    }
  }


  // 5.9. BIOENERGETICS AND CONSUMPTION
  // - Calculate temperature function of consumption
  calculate_temperature_function(fT, nspp, nyrs, Ceq, Cindex, Qc, Tcm, Tco, Tcl, CK1, CK4, env_index);

  // - Calculate historic ration
  calculate_ration(consumption_at_age, nspp, nyrs, nyrs_hind, nsex, nages, Ceq, CA, CB, fday, Pvalue, fT, ration_data, weight_hat);


  // 5.10. PARAMETRIC SUITABILITY
  // - Estimate vulnerability (Parametric Suitability)
  calculate_vulnerability(vulnerability, vulnerability_other, nspp, suitMode, log_phi);

  // - Parametric Suitability (Gamma, Lognormal, Normal)
  calculate_parametric_suitability(suitability, suit_other, nspp, nyrs, nsex, nages, suitMode,
                                   length_hat, weight_hat, vulnerability, vulnerability_other, gam_a, gam_b);


  // 5.11. REORGANIZE MSVPA DIET DATA
  // - Reorganize diet_obs content
  organize_diet_obs(nspp, nyrs, nyrs_hind, styr, minage, nsex,
                    diet_obs, diet_ctl, diet_prop);

  // - Calculate other food stomach content
  calculate_other_food_diet_prop(nspp, nyrs, nsex, nages, other_food,
                                 diet_prop, other_food_diet_prop);

  // 5.12. FISHING MORTALITY and FSPRs
  // Outside the multi-species iteration loop: hindcast F needs no M. The
  // forecast arm below is provisional -- section 6.7 recomputes it once the
  // HCR has its reference points.
  F_spp.setZero();
  F_flt_age.setZero();
  F_at_age.setZero();
  Ftarget_at_age.setZero();
  Flimit_at_age.setZero();
  for(flt = 0; flt < n_flt; flt++) {

    sp = flt_spp(flt);  // Species this fleet fishes or surveys

    if(flt_type(flt) == 1){
      for(age = 0; age < nages(sp); age++) {
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = 0; yr < nyrs; yr++) {

            // Hindcast
            if( yr < nyrs_hind){
              F_flt_age(flt, sex, age, yr) = sel_at_age(flt, sex, age, yr) * exp(log_F(flt, yr));
            }


            // Forecast
            if( yr >= nyrs_hind){
              // -- Apply HCRs
              switch(HCR){
              case 0: // No fishing
                proj_F(sp, yr) = 0;
                break;

              case 1: // CMSY
                proj_F(sp, yr) = Ftarget(sp);
                break;

              case 2: // Constant F
                proj_F(sp, yr) = Ftarget(sp);
                break;

              case 3: // Constant F that acheives X% of SSB0
                proj_F(sp, yr) = Ftarget(sp);
                break;

              case 4: // Constant Fspr
                proj_F(sp, yr) = Ftarget(sp) * Fmult(sp);
                break;

              case 5: // NPFMC Tier 3 HCR
                proj_F(sp, yr) = Ftarget(sp); // Used Fabc of Ftarget_age%
                break;

              case 6: // PFMC Category 1 HCR
                proj_F(sp, yr) = Ftarget(sp) = Flimit(sp) + QnormHCR(sp);
                break;

              case 7: // SESSF Tier 1 HCR
                proj_F(sp, yr) = Ftarget(sp); // Used Fabc of Ftarget_age%
                break;
              }

              // No projected F when the forecast is off or numbers-at-age are input (estDynamics > 0).
              if((forecast(sp) == 0) | (estDynamics(sp) > 0)){
                proj_F(sp, yr) = 0;
              }
              F_flt_age(flt, sex, age, yr) = sel_at_age(flt, sex, age, yr) * proj_F_prop(flt) * proj_F(sp, yr);
            }

            // -- Sum F across fleets
            F_at_age(sp, sex, age, yr) += F_flt_age(flt, sex, age, yr);

            // -- Calculate F target by age and sex for reference points
            Flimit_at_age(sp, sex, age, yr) += sel_at_age(flt, sex, age, yr) * proj_F_prop(flt) * Flimit(sp); // account for time-varying sel
            Ftarget_at_age(sp, sex, age, yr) += sel_at_age(flt, sex, age, yr) * proj_F_prop(flt) * Ftarget(sp); // account for time-varying sel
          }
        }
      }

      // F across fleets or species
      for(yr = 0; yr < nyrs; yr++) {
        // Hindcast
        if( yr < nyrs_hind){
          F_flt(flt, yr) = exp(log_F(flt, yr));
          F_spp(sp, yr) += exp(log_F(flt, yr)); // Fully selected fishing mortality
        }

        // Forecast
        if( yr >= nyrs_hind){
          F_flt(flt, yr) = proj_F_prop(flt) * proj_F(sp, yr);
          F_spp(sp, yr) +=  proj_F_prop(flt) * proj_F(sp, yr);
        }
      }
    }
  }


  // 5.12.1. FISHERY SELECTIVITY AT AGE FOR THE INITIAL STATE (initMode 6 only)
  // SS3's InitF convention: the initial equilibrium decays at Finit weighted by fishery
  // selectivity, so a size-selective fishery barely touches the young ages. Mean over the
  // species' fishery fleets in year 1, and 1 everywhere for every other mode. Filled
  // before the SPR reference points build Z_init, and so SPRFinit, from it.
  for(sp = 0; sp < nspp; sp++){
    for(sex = 0; sex < nsex(sp); sex++){
      for(age = 0; age < nages(sp); age++){
        sel_init(sp, sex, age) = 1.0;
      }
    }
  }
  if(initMode == 6){
    for(sp = 0; sp < nspp; sp++){
      int n_fsh = 0;
      for(flt = 0; flt < n_flt; flt++){
        if((flt_type(flt) == 1) && (flt_spp(flt) == sp)) n_fsh++;
      }
      if(n_fsh > 0){
        for(sex = 0; sex < nsex(sp); sex++){
          for(age = 0; age < nages(sp); age++){
            Type acc = 0.0;
            for(flt = 0; flt < n_flt; flt++){
              if((flt_type(flt) == 1) && (flt_spp(flt) == sp)){
                acc += sel_at_age(flt, sex, age, 0);
              }
            }
            sel_init(sp, sex, age) = acc / Type(n_fsh);
          }
        }
      }
    }
  }


  // 5.13. SIMULATE PROCESS ERROR (sim_mod(simulate = TRUE), simulate_state)
  // Drawn before the dynamics consume the deviations -- the initial-numbers block below
  // is the first reader of rec_dev and init_dev -- so the dynamics and every observation
  // draw downstream are consistent with the simulated process. No second copy of the
  // population dynamics is needed (WHAM's sim_pop), these deviations and their SDs all
  // being parameters. Each draw uses the density that scores it in section 13, bias
  // correction included, which is what makes simulated recruitment mean-unbiased.
  //
  // True if exactly one density scores rec_dev. False for the AMAK/Ianelli configuration
  // (srr_fun == 0 with srr_pred_fun > 0), where the stock-recruit curve is fitted as a
  // second penalty on the same deviation. REPORTed so sim_mod() can warn without
  // re-deriving it.
  int rec_srr_single_density = !((srr_fun == 0) && (srr_pred_fun > 0));
  REPORT(rec_srr_single_density);

  SIMULATE {
    // Masks marking which cells this call actually drew. Every deviation array
    // is REPORTed whole -- hindcast and projection alike -- so without these a
    // caller cannot tell a simulated cell from one left at its fitted value, and
    // a recovery statistic taken over the whole array scores both.
    matrix<Type> rec_dev_drawn_sim(rec_dev.rows(), rec_dev.cols());
    matrix<Type> init_dev_drawn_sim(init_dev.rows(), init_dev.cols());
    array<Type>  log_M1_dev_drawn_sim(log_M1_dev.dim);
    rec_dev_drawn_sim.setZero();
    init_dev_drawn_sim.setZero();
    log_M1_dev_drawn_sim.setZero();

    for(sp = 0; sp < nspp; sp++){

      // The density covers the hindcast only -- projection recruitment comes from the
      // harvest control rule or sample_rec() -- so the fitted window is the only period
      // this can honestly redraw. srr_fun == 0 with srr_pred_fun > 0
      // (rec_srr_single_density, above) is excluded: in that AMAK/Ianelli configuration
      // section 13 scores rec_dev through JNLL_REC_DEV and again through
      // JNLL_SRR_PENALTY, so there is no single distribution to draw from. Nothing is
      // drawn, sim_mod() says why, and the full argument is in
      // vignette("model-diagnostics").
      if(simulate_state(0) == 1 && simulate_period(0) == 1 && rec_srr_single_density){
        for(yr = 0; yr < nyrs_hind; yr++){
          rec_dev(sp, yr) = rnorm(-bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp));
          rec_dev_drawn_sim(sp, yr) = 1;
        }
      }

      // Initial abundance-at-age: the same process, for years before styr. The
      // gate is the density's own -- the equilibrium modes and OffsetEquilibrium
      // (5) fix init_dev and carry no penalty, so there is nothing to draw -- and
      // the rec_srr_single_density gate as well, a fresh initial age structure on
      // top of FITTED hindcast deviations being no history the model generated.
      if(simulate_state(0) == 1 && simulate_period(0) == 1 && rec_srr_single_density &&
         (initMode > 1) && (initMode != 5)){
        for(age = 1; age < nages(sp); age++){
          init_dev(sp, age - 1) = rnorm(-bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp));
          init_dev_drawn_sim(sp, age - 1) = 1;
        }
      }
    }

    // Natural mortality random effects. The IID modes are the AR1 code at rho = 0, so one
    // construction covers all six: draw a standardized AR1 and scale by the marginal SD
    // afterwards. M1_dev_log_sd is the INNOVATION sd; Sigma_M is the marginal.
    //
    // The draw is BROADCAST along whichever dimension the mode holds constant (modes 1/4
    // over years, modes 2/5 over ages) and across the sexes sharing one parameter --
    // build_map() separates the sexes only at M1_model = 2, so num_re_sexes is the number
    // of INDEPENDENT draws. This block writes the expanded array, so writing only the
    // cell the density reads would leave the rest at fitted values and the deviation
    // would not be constant.
    if(simulate_state(1) == 1 && simulate_period(0) == 1){
      for(sp = 0; sp < nspp; sp++){
        int num_re_sexes = (M1_model(sp) == 2 && nsex(sp) > 1) ? 2 : 1;

        if((M1_re(sp) == 1) || (M1_re(sp) == 4)){        // by age, constant over years
          Type sigma_M = exp(M1_dev_log_sd(sp, 0));
          Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
          Type Sigma_M = pow(pow(sigma_M, 2) / (1.0 - pow(rho_M_a, 2)), 0.5);
          for(int s = 0; s < num_re_sexes; s++){
            vector<Type> M_re_age(nages(sp));
            AR1(rho_M_a).simulate(M_re_age);
            M_re_age *= Sigma_M;
            int sex_end = (num_re_sexes == 1) ? nsex(sp) : s + 1;
            for(int sex = s; sex < sex_end; sex++){
              for(age = 0; age < nages(sp); age++){
                for(yr = 0; yr < nyrs; yr++){
                  log_M1_dev(sp, sex, age, yr) = M_re_age(age);
                  log_M1_dev_drawn_sim(sp, sex, age, yr) = 1;
                }
              }
            }
          }
        }

        if((M1_re(sp) == 2) || (M1_re(sp) == 5)){        // by year, constant over ages
          Type sigma_M = exp(M1_dev_log_sd(sp, 0));
          Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
          Type Sigma_M = pow(pow(sigma_M, 2) / (1.0 - pow(rho_M_y, 2)), 0.5);
          for(int s = 0; s < num_re_sexes; s++){
            vector<Type> M_re_yr(nyrs_hind);
            AR1(rho_M_y).simulate(M_re_yr);
            M_re_yr *= Sigma_M;
            int sex_end = (num_re_sexes == 1) ? nsex(sp) : s + 1;
            for(int sex = s; sex < sex_end; sex++){
              for(yr = 0; yr < nyrs_hind; yr++){
                for(age = 0; age < nages(sp); age++){
                  log_M1_dev(sp, sex, age, yr) = M_re_yr(yr);
                  log_M1_dev_drawn_sim(sp, sex, age, yr) = 1;
                }
              }
            }
          }
        }

        if((M1_re(sp) == 3) || (M1_re(sp) == 6)){        // by age and year
          Type sigma_M = exp(M1_dev_log_sd(sp, 0));
          Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
          Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
          Type Sigma_M = pow(pow(sigma_M, 2) / ((1.0 - pow(rho_M_y, 2)) * (1.0 - pow(rho_M_a, 2))), 0.5);
          for(int s = 0; s < num_re_sexes; s++){
            array<Type> M_re_a_yr(nages(sp), nyrs_hind);
            // Same argument order as the density in section 14 -- year on the
            // outermost dimension, age on the fastest-running one.
            SEPARABLE(AR1(rho_M_y), AR1(rho_M_a)).simulate(M_re_a_yr);
            int sex_end = (num_re_sexes == 1) ? nsex(sp) : s + 1;
            for(int sex = s; sex < sex_end; sex++){
              for(age = 0; age < nages(sp); age++){
                for(yr = 0; yr < nyrs_hind; yr++){
                  log_M1_dev(sp, sex, age, yr) = Sigma_M * M_re_a_yr(age, yr);
                  log_M1_dev_drawn_sim(sp, sex, age, yr) = 1;
                }
              }
            }
          }
        }
      }
    }

    // Reported so the draw can be checked directly. Inferring it from the data
    // it produces cannot tell "drew correctly" from "drew nothing, and the
    // observation error moved". Each goes out with its drawn-cell mask.
    matrix<Type> rec_dev_sim = rec_dev;
    matrix<Type> init_dev_sim = init_dev;
    array<Type>  log_M1_dev_sim = log_M1_dev;
    REPORT(rec_dev_sim);
    REPORT(init_dev_sim);
    REPORT(log_M1_dev_sim);
    REPORT(rec_dev_drawn_sim);
    REPORT(init_dev_drawn_sim);
    REPORT(log_M1_dev_drawn_sim);
  }

  /** ------------------------------------------------------------------------ //
   * 6. POPULATION DYNAMICS EQUATIONS                                          //
   * ------------------------------------------------------------------------- */
  // Start iterations for multi-species convergence (anything dependent on M)
  for(int iter = 0; iter < niter; iter++) {

    // 6.1. TOTAL MORTALITY-AT-AGE
    for(sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        // Element-wise: assigning a scalar to a vector would broadcast.
        for(int i = 0; i < M1_beta.dim(2); i++){
          beta_M1_tmp(i) = M1_beta(sp, sex, i);
        }
        for(age = 0; age < nages(sp); age++) {
          for(yr = 0; yr < nyrs; yr++) {
            // Matrix multiplication from sliced arrays doesn't work
            env_M1_tmp = env_index.row(yr);
            M1_mult = env_M1_tmp * beta_M1_tmp;

            // M1_mult.sum() is LEGACY (scheduled removal: v4.5.0). It carries
            // the environmental effect of the soft-deprecated M1_indices path,
            // so fits built that way still reproduce; `linkages` is the current
            // route. See "Scheduled removal" in NEWS.md (4.1.0 section).
            M1_at_age(sp, sex, age, yr) = exp(
              log_M1(sp, sex, age)
              + log_M1_dev(sp, sex, age, yr)
              + M1_mult.sum()
              + M_linkage_offset(sp, sex, age, yr)
            ) + M_linkage_offset_nat(sp, sex, age, yr);
            M_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
            M_at_age_dB0(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age_dB0(sp, sex, age, yr);
            M_at_age_dBF(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age_dBF(sp, sex, age, yr);
            Z_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + F_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
          }
        }
      }
    }


    // 6.2. SPR BASED REFERENCE POINTS
    //FIXME - make time-varying?
    SPR0.setZero();
    SPRFinit.setZero();
    SPRlimit.setZero();
    SPRtarget.setZero();
    // Not defined under predation: M carries M2, which scales with predator abundance, so
    // spawning output per recruit would not be a property of the prey stock alone.
    // data_check() refuses the switch combinations that read these downstream.
    if(msmMode == 0){
      for(sp = 0; sp < nspp; sp++) {

        int na = nages(sp);
        int term_yr = nyrs_hind - 1;      // rates for a reference point are the terminal hindcast year's
        wt_idx_ssb = 2 * sp + 1;

        // The four schedules. Flimit and Ftarget are selectivity-scaled and read the
        // terminal hindcast year; Finit is the flat F that generated the initial age
        // structure, so it is constant across ages and reads year 0, weight included.
        // FIXME: time-vary sel in the forecast
        vector<Type> Z_unfished(na), Z_limit(na), Z_target(na), Z_init(na);
        vector<Type> wt_term(na), wt_first(na), mature_at_age(na);
        Type female_split = (nsex(sp) == 1) ? Type(1.0) : sex_ratio(sp, 0);

        for(age = 0; age < na; age++){
          Z_unfished(age) = M_at_age(sp, 0, age, term_yr);
          Z_limit(age)    = M_at_age(sp, 0, age, term_yr) + Flimit_at_age(sp, 0, age, term_yr);
          Z_target(age)   = M_at_age(sp, 0, age, term_yr) + Ftarget_at_age(sp, 0, age, term_yr);
          Z_init(age)     = M_at_age(sp, 0, age, 0) + Finit(sp) * sel_init(sp, 0, age);

          // Spawning output per TOTAL recruit, so the female fraction enters once:
          // spawn_output (5.7) carries maturity and the age-varying ratio for a
          // one-sex species, female_split the recruitment split (6.6) for a two-sex one.
          wt_term(age)      = spawn_output(sp, age, term_yr);
          wt_first(age)     = spawn_output(sp, age, 0);
          mature_at_age(age) = female_split;
        }

        vector<Type> n_unfished = per_recruit_survivors(Z_unfished);
        vector<Type> n_limit    = per_recruit_survivors(Z_limit);
        vector<Type> n_target   = per_recruit_survivors(Z_target);
        vector<Type> n_init     = per_recruit_survivors(Z_init);

        SPR0(sp)      = spawning_biomass_per_recruit(n_unfished, Z_unfished, wt_term,  mature_at_age, spawn_month(sp));
        SPRlimit(sp)  = spawning_biomass_per_recruit(n_limit,    Z_limit,    wt_term,  mature_at_age, spawn_month(sp));
        SPRtarget(sp) = spawning_biomass_per_recruit(n_target,   Z_target,   wt_term,  mature_at_age, spawn_month(sp));
        SPRFinit(sp)  = spawning_biomass_per_recruit(n_init,     Z_init,     wt_first, mature_at_age, spawn_month(sp));

        // Reported for inspection; slot order matches the four calls above.
        for(age = 0; age < na; age++){
          NbyageSPR(0, sp, age) = n_unfished(age);
          NbyageSPR(1, sp, age) = n_limit(age);
          NbyageSPR(2, sp, age) = n_target(age);
          NbyageSPR(3, sp, age) = n_init(age);
        }
      }
    }


    // 6.3. INITIAL RECRUITMENT
    // -- For beverton-holt, steepness and R0 are derived from SPR0
    penalty = 0.0;
    zero_N_pen.setZero();
    for( sp = 0; sp < nspp ; sp++) {
      // Under predation SPR is undefined: a hindcast curve takes R_init from the
      // free level R0(sp, 0) and reports steepness 0, skipping the cases below.
      if((msmMode > 0) & (srr_fun > 1)){
        for(yr = 0; yr < nyrs; yr++){ steepness(sp, yr) = 0; }
        R_init(sp) = R0(sp, 0);
      } else switch(srr_fun){
      case 0: // Random about mean (e.g. Alaska)
        // No compensation, so steepness is constant across years.
        for(yr = 0; yr < nyrs; yr++){ steepness(sp, yr) = 0.99; }
        {
          R_init(sp) = R0(sp, 0);
        }
        break;

      case 1: // Random about mean with environmental linkage
        for(yr = 0; yr < nyrs; yr++){ steepness(sp, yr) = 0.99; }
        {
          R_init(sp) = R0(sp, 0);
        }
        break;

      case 2: // Beverton-Holt
        {
          // Steepness for every year: alpha may be time-varying through a linkage.
          for(yr = 0; yr < nyrs; yr++){
            steepness(sp, yr) = alpha(sp, yr) * SPR0(sp)/(4.0 + alpha(sp, yr) * SPR0(sp));
          }
          // Unfished recruitment implied by the curve, first year only. Filling every
          // year would make R0 an AD function of alpha and Beta across the whole series
          // -- ~6x the fit time for a value that is reported rather than used,
          // Beverton-Holt recruitment being a function of alpha, Beta and SSB, not of R0.
          R0(sp, 0) = (alpha(sp, 0) - 1.0/SPR0(sp)) / Beta(sp, 0); // (Alpha-1/SPR0)/beta
          R_init(sp) = (alpha(sp, 0) - 1.0/SPRFinit(sp)) / Beta(sp, 0); // (Alpha-1/SPRFinit)/beta
          // Same floor as R in 6.6.1: an identity-link alpha offset can make both negative.
          if (rec_floor_on) { Type pen_R = 0; R0(sp, 0) = posfun(R0(sp, 0), Type(1e-3), pen_R); R_init(sp) = posfun(R_init(sp), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }
        }
        break;

      case 3: // Retired alias of 2; build_srr() refuses it (alpha moves via a linkage)
        {
          // Steepness for every year: alpha may be time-varying through a linkage.
          for(yr = 0; yr < nyrs; yr++){
            steepness(sp, yr) = alpha(sp, yr) * SPR0(sp)/(4.0 + alpha(sp, yr) * SPR0(sp));
          }
          // Unfished recruitment implied by the curve, first year only -- see
          // case 2 for why it is not filled for every year.
          R0(sp, 0) = (alpha(sp, 0) - 1.0/SPR0(sp)) / Beta(sp, 0); // (Alpha-1/SPR0)/beta
          R_init(sp) = (alpha(sp, 0) - 1.0/SPRFinit(sp)) / Beta(sp, 0); // (Alpha-1/SPRFinit)/beta
          // Same floor as R in 6.6.1: an identity-link alpha offset can make both negative.
          if (rec_floor_on) { Type pen_R = 0; R0(sp, 0) = posfun(R0(sp, 0), Type(1e-3), pen_R); R_init(sp) = posfun(R_init(sp), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }
        }
        break;

      case 4: // Ricker
        {
          // Steepness for every year: alpha may be time-varying through a linkage.
          for(yr = 0; yr < nyrs; yr++){
            { Type aS = alpha(sp, yr) * SPR0(sp);   // kept positive under an identity-link alpha offset
              if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
              steepness(sp, yr) = 0.2 * exp(0.8*log(aS)); }
          }

          // - R at F0
          Type pen_rk = 0;   // this species' floor penalties alone
          ricker_intercept = alpha(sp, 0) * SPR0(sp) - 1.0;
          ricker_intercept =  posfun(ricker_intercept, Type(0.001), pen_rk) + 1.0;
          // Year 0 only: posfun() keeps the Ricker intercept positive, charging the
          // excursion to the objective. The Ianelli configuration (srr_fun < 2)
          // never reaches this case, so its R0 keeps its estimated value.
          R0(sp, 0) = log(ricker_intercept)/(Beta(sp, 0) * SPR0(sp)/1000000.0);

          // R at equilibrium F
          ricker_intercept = alpha(sp, 0) * SPRFinit(sp) - 1.0;
          ricker_intercept =  posfun(ricker_intercept, Type(0.001), pen_rk) + 1.0;
          R_init(sp) = log(ricker_intercept)/(Beta(sp, 0) * SPRFinit(sp)/1000000.0);
          zero_N_pen(sp) += pen_rk;
        }
        break;

      case 5: // Ricker with environmental impacts on alpha
        {
          // Steepness for every year: alpha may be time-varying through a linkage.
          for(yr = 0; yr < nyrs; yr++){
            { Type aS = alpha(sp, yr) * SPR0(sp);   // kept positive under an identity-link alpha offset
              if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
              steepness(sp, yr) = 0.2 * exp(0.8*log(aS)); }
          }

          // - R at F0
          Type pen_rk = 0;   // this species' floor penalties alone
          ricker_intercept = alpha(sp, 0) * SPR0(sp) - 1.0;
          ricker_intercept =  posfun(ricker_intercept, Type(0.001), pen_rk) + 1.0;
          // Year 0 only, as case 4.
          R0(sp, 0) = log(ricker_intercept)/(Beta(sp, 0) * SPR0(sp)/1000000.0);

          // R at equilibrium F
          ricker_intercept = alpha(sp, 0) * SPRFinit(sp) - 1.0;
          ricker_intercept =  posfun(ricker_intercept, Type(0.001), pen_rk) + 1.0;
          R_init(sp) = log(ricker_intercept)/(Beta(sp, 0) * SPRFinit(sp)/1000000.0);
          zero_N_pen(sp) += pen_rk;
        }
        break;

      default:
        error("Invalid 'srr_fun'");
      }

      // Under the Ianelli configuration (srr_fun 0 or 1, srr_pred_fun > 1) the curve is
      // estimated as a recruitment penalty rather than driving the hindcast, and the
      // switch above -- keyed on srr_fun -- would leave steepness at the mean-recruitment
      // constant. Derive it from the penalty curve's alpha instead, since the
      // stock-recruit prior (JNLL_SRR_PRIOR) is evaluated on it. R0 is deliberately NOT
      // re-derived: recruitment here is R0 * exp(rec_dev), so R0 keeps its estimated
      // value.
      if((srr_fun < 2) & (srr_pred_fun > 1) & (msmMode == 0)){
        for(yr = 0; yr < nyrs; yr++){
          if((srr_pred_fun == 2) | (srr_pred_fun == 3)){
            steepness(sp, yr) = alpha(sp, yr) * SPR0(sp)/(4.0 + alpha(sp, yr) * SPR0(sp));
          }
          if((srr_pred_fun == 4) | (srr_pred_fun == 5)){
            { Type aS = alpha(sp, yr) * SPR0(sp);   // kept positive under an identity-link alpha offset
              if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
              steepness(sp, yr) = 0.2 * exp(0.8*log(aS)); }
          }
        }
      }
      // SPR0 is undefined under predation, so steepness is reported as 0;
      // data_check() refuses a steepness prior there.
      if((srr_fun < 2) & (srr_pred_fun > 1) & (msmMode > 0)){
        for(yr = 0; yr < nyrs; yr++){ steepness(sp, yr) = 0; }
      }
    }


    // 6.4. INITIAL ABUNDANCE AT AGE, BIOMASS, AND SSB (YEAR 1)
    biomass.setZero();
    ssb.setZero();
    for(sp = 0; sp < nspp; sp++) {

      // Sex ratio at recruitment (1-sex model gets all R)
      if(nsex(sp)  == 1){
        sex_ratio(sp, 0) = 1.0;
      }

      // A free initial recruitment level, log scale; 0 without an R_init linkage.
      // Read at year 0 because the initial state is one year, and per species
      // because that is the stratum it is estimated on. Only the estimated branch
      // reads it: input numbers-at-age have no initial state to scale.
      Type init_level = recruitment_linkage_offset(sp, RCEATTLE_REC_R_INIT, 0);

      for(age = 0; age < nages(sp); age++){
        for(sex = 0; sex < nsex(sp); sex ++){

          switch(estDynamics(sp)){
          case 0: // Estimated

            // - Amin (i.e. recruitment).
            if(age == 0){
              // Year-1 recruitment is a recruitment year with its own deviate,
              // not part of the pre-hindcast level, so it stays at R_init.
              R(sp, 0) = R_init(sp) * exp(rec_dev(sp, 0));
              N_at_age(sp, 0, 0, 0) = R(sp, 0) * sex_ratio(sp, 0);
              // Guarded because the male slot exists only where some species is
              // two-sex: the value would be 0 anyway (sex_ratio 1, set above),
              // but the index is out of range when max_sex is 1.
              if(nsex(sp) > 1){
                N_at_age(sp, 1, 0, 0) = R(sp, 0) * (1-sex_ratio(sp, 0));
              }
              // The equilibrium carries no recruitment deviation: it is the mean
              // level the deviates depart from. It DOES carry the initial level,
              // on every age including this one, so the equilibrium that
              // equil_catch_hat integrates sits at one recruitment level.
              N_eq(sp, 0, 0) = R_init(sp) * exp(init_level) * sex_ratio(sp, 0);
              if(nsex(sp) > 1){
                N_eq(sp, 1, 0) = R_init(sp) * exp(init_level) * (1-sex_ratio(sp, 0));
              }
            }

            // - Estimate  as free parameters
            if(initMode == 0){
              if(age > 0){
                N_at_age(sp, 0, age, 0) = exp(init_dev(sp, age-1)) * sex_ratio(sp, 0);
                if(nsex(sp) > 1){
                  N_at_age(sp, 1, age, 0) = exp(init_dev(sp, age-1)) * (1-sex_ratio(sp, 0));
                }
                // No equilibrium underlies free initial numbers; Finit is 0 here
                // and an equilibrium catch is refused, so this is never fitted.
                N_eq(sp, 0, age) = N_at_age(sp, 0, age, 0);
                if(nsex(sp) > 1){
                  N_eq(sp, 1, age) = N_at_age(sp, 1, age, 0);
                }
              }
            }

            // - Equilibrium or non-equilibrium from Rinit, Finit, mortality and init devs
            // Finit is 0 for initMode 0, 1, 2 and 5; it is estimated for the
            // fished non-equilibrium modes 3, 4 and 6 (see section 6.1).
            if(initMode > 0){

              // How far the initial age-structure sits off the mean-recruitment
              // equilibrium R0, on the log scale; 0 under neither source below. Two
              // sources, never both (an R_init linkage is refused under mode 5): the free
              // R_init level -- SS3's SR_regime block on the year before the hindcast,
              // which SS3 shrinks at sigmaR/ave_age and Rceattle leaves free, so put a
              // prior on the coefficient to match SS3 -- or, under OffsetEquilibrium (5),
              // the first-year recruitment deviation, seeding the structure off
              // exp(rec_pars + rec_dev(sp, 0)) with init devs off (Cole Monnahan / AFSC
              // GOA pollock convention).
              Type init_log_scalar = init_level;
              if(initMode == 5){
                init_log_scalar += rec_dev(sp, 0);
              }

              // Sum M1 until age - 1. OffsetEquilibrium (5) uses the same
              // standard departing-age cumulative-M decay as the equilibrium
              // modes (Finit is 0 for mode 5).
              if((initMode == 1) | (initMode == 2) | (initMode == 3) | (initMode == 5)){
                mort_sum(sp, age) = 0;
                for(int age_tmp = 0; age_tmp < age; age_tmp++){
                  mort_sum(sp, age) += M1_at_age(sp, sex, age_tmp, 0) + Finit(sp);
                }
              }

              // SS3's InitF convention: the initial F is spread over ages by
              // the fishery's selectivity before it accumulates, so an
              // unselected young age carries almost none of it.
              if(initMode == 6){
                mort_sum(sp, age) = 0;
                for(int age_tmp = 0; age_tmp < age; age_tmp++){
                  mort_sum(sp, age) += M1_at_age(sp, sex, age_tmp, 0)
                                     + Finit(sp) * sel_init(sp, sex, age_tmp);
                }
              }

              if(initMode == 4){
                mort_sum(sp, age) = 0;
                for(int age_tmp = 0; age_tmp < age; age_tmp++){
                  mort_sum(sp, age) += M1_at_age(sp, sex, age_tmp, 0);
                }
                mort_sum(sp, age) += Finit(sp);
              }

              // -- 6.5.2. Age Amin+1:Amax-1 (initial abundance)
              if((age > 0) & (age < nages(sp) - 1)) {

                if(sex == 0){
                  N_at_age(sp, 0, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1) + init_log_scalar) * sex_ratio(sp, 0);
                  N_eq(sp, 0, age)        = R_init(sp) * exp( - mort_sum(sp, age) + init_log_scalar) * sex_ratio(sp, 0);
                }
                if(sex == 1){
                  N_at_age(sp, 1, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1) + init_log_scalar) * (1-sex_ratio(sp, 0));
                  N_eq(sp, 1, age)        = R_init(sp) * exp( - mort_sum(sp, age) + init_log_scalar) * (1-sex_ratio(sp, 0));
                }
              }

              // -- 6.5.3. Amax
              if(age == (nages(sp) - 1)) {

                // The plus group accumulates at its own total mortality. Under initMode
                // 6 the initial F reaches an age through that age's selectivity, so the
                // oldest age decays at M1 + Finit * sel (SS3's equilibrium convention);
                // every other mode charges the full Finit, as in mort_sum above.
                Type Z_plus = M1_at_age(sp, sex, nages(sp) - 1, 0) + Finit(sp) *
                  ((initMode == 6) ? sel_init(sp, sex, nages(sp) - 1) : Type(1.0));
                Type plus_grp = 1 - exp(-Z_plus);

                if(sex == 0){// NOTE: This solves for the geometric series
                  N_at_age(sp, 0, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1) + init_log_scalar) / plus_grp * sex_ratio(sp, 0);
                  N_eq(sp, 0, age)        = R_init(sp) * exp( - mort_sum(sp, age) + init_log_scalar) / plus_grp * sex_ratio(sp, 0);
                }

                if(sex == 1){
                  N_at_age(sp, 1, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1) + init_log_scalar) / plus_grp * (1-sex_ratio(sp, 0));
                  N_eq(sp, 1, age)        = R_init(sp) * exp( - mort_sum(sp, age) + init_log_scalar) / plus_grp * (1-sex_ratio(sp, 0));
                }
              }
            }
            break;

          case 1: // Numbers-at-age fixed exactly to NByageFixed (pop_scalar mapped to
            // NA in build_map so log_pop_scalar = 0 -> pop_scalar = 1.0)
            N_at_age(sp, sex, age, 0) = pop_scalar(sp) * NByageFixed(sp, sex, age, 0);
            break;

          case 2: // Numbers-at-age scaled by a single estimated age-independent scalar
            N_at_age(sp, sex, age, 0) = pop_scalar(sp) * NByageFixed(sp, sex, age, 0);
            break;

          default:
            error("Invalid 'estDynamics'");
          }

          // -- 6.5.3. Estimate total biomass in year 1
          wt_idx_pop = 2 * sp;
          biomass_at_age(sp, sex, age, 0) = N_at_age(sp, sex, age, 0) * weight_hat( wt_idx_pop, sex, age, 0);
          biomass(sp, 0) += N_at_age(sp, sex, age, 0) * weight_hat( wt_idx_pop, sex, age, 0);
        }

        // -- 6.5.4. Estimated initial female SSB
        wt_idx_ssb = 2 * sp + 1;
        ssb(sp, 0) += N_at_age(sp, 0, age, 0) * exp(-Z_at_age(sp, 0, age, 0) * spawn_month(sp)/12) * spawn_output(sp, age, 0); // 6.6.
      }
    }


    // 6.5. HINDCAST NUMBERS AT AGE, BIOMASS-AT-AGE (mt), and ssb-AT-AGE (mt)
    penalty = 0.0;
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 1; yr < nyrs_hind; yr++) {

        // Switch for srr (for MSEs if using Ianelli SRR method)
        int srr_switch = srr_fun;
        if(yr >= nyrs_srrmean){
          srr_switch = srr_pred_fun;
        }

        // -- 6.6.1. Recruitment

        // Recruits arriving at age minage in year yr were spawned in yr - minage,
        // before styr while yr < minage. With no modelled SSB those years take
        // R_init * exp(rec_dev) -- equilibrium recruitment at F = Finit, what
        // year 0 gets, so ages seeded from the same equilibrium agree. Stock
        // Synthesis treats the pre-start period the same way; WHAM avoids it by
        // fixing the lag at one year. Cannot fire at minage = 1.
        int spawn_yr = yr - minage(sp);
        int srr_use = (spawn_yr < 0) ? 0 : srr_switch;
        Type ssb_tmp = (spawn_yr < 0) ? Type(0.0) : ssb(sp, spawn_yr);
        // Pre-styr years take R_init, not R0(sp, yr): under a hindcast curve R0 is
        // set only at year 0 (derived, or the free level under predation).
        Type rec_mean = (spawn_yr < 0) ? R_init(sp) : R0(sp, yr);

        R(sp, yr) = calculate_recruitment(srr_use, rec_mean, ssb_tmp, alpha(sp, yr), Beta(sp, yr), rec_dev(sp, yr), SPR0(sp));
        // Under an identity-link recruitment offset the curve can go non-positive: keep R
        // positive (posfun, 0.001 barrier), charging the excursion to the zero-N row.
        if (rec_floor_on) { Type pen_R = 0; R(sp, yr) = posfun(R(sp, yr), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }

        N_at_age(sp, 0, 0, yr) = R(sp, yr) * sex_ratio(sp, 0);
        if(nsex(sp) > 1){
          N_at_age(sp, 1, 0, yr) = R(sp, yr) * (1.0-sex_ratio(sp, 0));
        }


        // -- 6.6.2. Ages beyond recruitment
        for(age = 0; age < nages(sp); age++){
          for(sex = 0; sex < nsex(sp); sex ++){

            switch(estDynamics(sp)){
            case 0: // Estimated numbers-at-age

              // -- Where Amin < age < Amax
              if(age < (nages(sp) - 1)) {
                N_at_age(sp, sex, age + 1, yr) = N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1));// S_at_age(sp, sex, age, yr - 1);
              }

              // -- Plus group where age = Amax.
              if(age == (nages(sp) - 1)) {
                N_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age - 1, yr - 1) * exp(-Z_at_age(sp, sex, age-1, yr-1)) + N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1)); // S_at_age(sp, sex, age, yr - 1);
              }
              break;

            case 1: // Numbers-at-age fixed exactly to NByageFixed (pop_scalar = 1.0 via map)
              N_at_age(sp, sex, age, yr) = pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              break;

            case 2: // Numbers-at-age scaled by a single estimated age-independent scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              break;

            default:
              error("Invalid 'estDynamics'");
            }

            // A fresh accumulator per cell: the penalty is this cell's excursion alone.
            { Type pen_N = 0; N_at_age(sp, sex, age, yr) = posfun(N_at_age(sp, sex, age, yr), Type(0.001), pen_N); zero_N_pen(sp) += pen_N; }

            // -- 6.6.3. Estimate total biomass
            wt_idx_pop = 2 * sp ;
            biomass_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * weight_hat( wt_idx_pop, sex, age, yr );
            biomass(sp, yr) += biomass_at_age(sp, sex, age, yr);
          }

          // -- 6.6.4. Estimated female ssb
          wt_idx_ssb = 2 * sp + 1;
          ssb(sp, yr) += N_at_age(sp, 0, age, yr) * exp(-Z_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0) * spawn_output(sp, age, yr); // 6.6.
        }
      }
    }


    // 6.6. DEPLETION BASED BIOMASS REFERENCE POINTS (i.e. SB0 and dynamic SB0)
    // -- calculate mean recruitment
    // FIXME: may want specified year ranges
    avg_R.setZero();
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 0; yr < nyrs_srrmean; yr++) {
        avg_R(sp) += R(sp, yr)/Type(nyrs_srrmean);
      }
    }

    // -- Calc biomass_depletion based BRPs
    NByage0.setZero();
    NByageF.setZero();
    N_at_age_dBF.setZero();
    N_at_age_dB0.setZero();
    B0.setZero();
    SB0.setZero();
    SBF.setZero();
    DynamicB0.setZero();
    DynamicSB0.setZero();
    DynamicSBF.setZero();
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 0; yr < nyrs; yr++) {

        // - Year 1
        // -- Initial abundance-at-age is the same as the hindcast
        if(yr == 0){
          for(sex = 0; sex < nsex(sp); sex ++){
            for(age = 0; age < nages(sp); age++){
              NByage0(sp, sex, age , 0) = NByageF(sp, sex, age , 0) = N_at_age_dBF(sp, sex, age, 0) = N_at_age_dB0(sp, sex, age, 0) = N_at_age(sp,sex,age,0);
            }
          }
        }

        // - Recruitment Year > 1
        if(yr > 0){
          // - Option 1a: Use mean rec
          if((proj_mean_rec == 1) & (srr_pred_fun < 2)){
            // Equilibrium RPs used mean rec across all years
            NByage0(sp, 0, 0, yr) = avg_R(sp);
            NByageF(sp, 0, 0, yr) = avg_R(sp);

            // Dynamic RPs use mean rec for the forecast and observed rec for the
            // hindcast. LEGACY (scheduled removal: v4.5.0): prefer the
            // stock-recruit arm below. This split feeds SB0(sp, nyrs-1), so it
            // wants a deliberate decision rather than a sweep.
            if(yr < nyrs_hind){
              N_at_age_dBF(sp, 0, 0, yr) = N_at_age_dB0(sp, 0, 0, yr) = R(sp, yr); // Hindcast use observed R
            }
            if(yr >= nyrs_hind){
              N_at_age_dBF(sp, 0, 0, yr) = N_at_age_dB0(sp, 0, 0, yr) = exp(log(avg_R(sp)) + rec_dev(sp, yr)); // Projections use mean R given bias in R0
            }
          }


          // - Option 2: Use SRR
          // Reached whenever a curve exists, not only when the projection runs on it. The
          // breadth is deliberate: Option 1a excludes srr_pred_fun >= 2, so proj_mean_rec
          // = 1 (the build_srr() default) with an estimated curve would match NEITHER
          // arm, leaving reference-point recruitment 0 after the first year. SB0 would
          // decay towards zero, and SB0(sp, nyrs-1) is the depletion reference HCRs 5 and
          // 6 read. Projection recruitment is unaffected: 6.8 reads that switch
          // separately, so mean-recruitment projections stay mean.
          if((proj_mean_rec == 0) | (srr_pred_fun >= 2)){

            // - Equilibrium reference points (No recruitment deviation: pass Type(0.0))
            // Equilibrium spawning biomass exists every year; yr < minage reads the
            // first year's value so the lag stays in bounds.
            int rp_yr = yr - minage(sp);
            if(rp_yr < 0){ rp_yr = 0; }

            // TODO: not floored under an identity-link recruitment offset, so SB0, SBF and B0 can go negative (issue #234).
            NByage0(sp, 0, 0, yr) = calculate_recruitment(srr_pred_fun, R0(sp, yr), SB0(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), Type(0.0), SPR0(sp));
            NByageF(sp, 0, 0, yr) = calculate_recruitment(srr_pred_fun, R0(sp, yr), SBF(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), Type(0.0), SPR0(sp));

            // -  Dynamic reference points (Includes annual recruitment deviation: pass rdev)
            // Cohorts spawned before the first year take the hindcast's realized R
            // (R_init with its deviation), so both runs start from the same cohorts.
            Type rdev = rec_dev(sp, yr);
            if(yr < minage(sp)){
              Type R_early = (yr < nyrs_hind) ? R(sp, yr) : R_init(sp) * exp(rdev);
              N_at_age_dB0(sp, 0, 0, yr) = R_early;
              N_at_age_dBF(sp, 0, 0, yr) = R_early;
            } else {
              // Under the Ianelli form the hindcast is on mean recruitment until
              // srr_mse_switchyr, so rec_dev there is measured from R0; the curve takes
              // the realized deviation from it, log R - log R_hat.
              if((srr_fun != srr_pred_fun) & (yr < nyrs_srrmean)){
                Type R_curve = calculate_recruitment(srr_pred_fun, R0(sp, yr), ssb(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), Type(0.0), SPR0(sp));
                // Same floor as R in 6.6.1, for an identity-link offset that makes the curve non-positive.
                if (rec_floor_on) { Type pen_R = 0; R_curve = posfun(R_curve, Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }
                rdev = log(R(sp, yr)) - log(R_curve);
              }
              // TODO: dynamic-B0 recruitment is not floored under an identity-link offset (issue #234).
              N_at_age_dB0(sp, 0, 0, yr) = calculate_recruitment(srr_pred_fun, R0(sp, yr), DynamicSB0(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), rdev, SPR0(sp));
              N_at_age_dBF(sp, 0, 0, yr) = calculate_recruitment(srr_pred_fun, R0(sp, yr), DynamicSBF(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), rdev, SPR0(sp));
            }

          } // End recruitment switch

          // Account for sex ratio. The male block is guarded because that slot
          // exists only where some species is two-sex; the value would be 0
          // anyway (sex_ratio 1, set in 6.4).
          NByage0(sp, 0, 0, yr) = NByage0(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          NByageF(sp, 0, 0, yr) = NByageF(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          N_at_age_dB0(sp, 0, 0, yr) = N_at_age_dB0(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          N_at_age_dBF(sp, 0, 0, yr) = N_at_age_dBF(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females

          if(nsex(sp) > 1){
            // Males (divide by sex_ratio because we multiplied in the lines above)
            NByage0(sp, 1, 0, yr) = NByage0(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0));
            NByageF(sp, 1, 0, yr) = NByageF(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0));
            N_at_age_dB0(sp, 1, 0, yr) = N_at_age_dB0(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0));
            N_at_age_dBF(sp, 1, 0, yr) = N_at_age_dBF(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0));
          }


          // N-at-age year > 1
          for(sex = 0; sex < nsex(sp); sex ++){
            for(age = 1; age < nages(sp)-1; age++) {
              NByage0(sp, sex, age, yr) =  NByage0(sp, sex, age-1, yr-1) * exp(-M_at_age(sp, sex, age-1, yr-1)); // F = 0

              NByageF(sp, sex, age, yr) =  NByageF(sp, sex, age-1, yr-1) * exp(-M_at_age(sp, sex, age-1, yr-1) - Ftarget_at_age(sp, sex, age-1, yr-1)); // F = target

              // No 0.001 floor here: the hindcast's is a numerical guard, and dynamic B0
              // must not assume a stock that fell to it would not have crashed without
              // fishing.
              N_at_age_dB0(sp, sex, age, yr) =  N_at_age_dB0(sp, sex, age-1, yr-1) * exp(-M_at_age_dB0(sp, sex, age-1, yr - 1)); // F = 0

              N_at_age_dBF(sp, sex, age, yr) =  N_at_age_dBF(sp, sex, age-1, yr-1) * exp(-M_at_age_dBF(sp, sex, age-1, yr - 1) - Ftarget_at_age(sp, sex, age-1, yr-1)); // F = Ftarget
            }

            // Plus group  -- No fishing
            NByage0(sp, sex, nages(sp)-1, yr) = NByage0(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-2, yr - 1))  + NByage0(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-1, yr - 1));

            NByageF(sp, sex, nages(sp)-1, yr) = NByageF(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-2, yr - 1) - Ftarget_at_age(sp, sex, nages(sp)-2, yr-1))  + NByageF(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-1, yr - 1) - Ftarget_at_age(sp, sex, nages(sp)-1, yr-1));

            N_at_age_dB0(sp, sex, nages(sp)-1, yr) = N_at_age_dB0(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age_dB0(sp, sex, nages(sp)-2, yr - 1))  + N_at_age_dB0(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age_dB0(sp, sex, nages(sp)-1, yr - 1));

            N_at_age_dBF(sp, sex, nages(sp)-1, yr) = N_at_age_dBF(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age_dBF(sp, sex, nages(sp)-2, yr - 1) - Ftarget_at_age(sp, sex, nages(sp)-2, yr-1))  + N_at_age_dBF(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age_dBF(sp, sex, nages(sp)-1, yr - 1) - Ftarget_at_age(sp, sex, nages(sp)-1, yr-1));

          }
        }

        // A species with input numbers-at-age (estDynamics > 0) keeps them in the
        // dynamic runs; projection years are read directly, as 6.8 fills N later.
        if(estDynamics(sp) > 0){
          for(sex = 0; sex < nsex(sp); sex ++){
            for(age = 0; age < nages(sp); age++){
              Type n_fixed = (yr < nyrs_hind) ? N_at_age(sp, sex, age, yr) :
                pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              N_at_age_dB0(sp, sex, age, yr) = n_fixed;
              N_at_age_dBF(sp, sex, age, yr) = n_fixed;
            }
          }
        }


        // Calculate Dynamic SB0 and SB at F target
        for(age = 0; age < nages(sp); age++) {

          wt_idx_ssb = 2 * sp + 1;
          // Multispecies: M_at_age carries the projection's realized M2. The one rule
          // that reads SBF, NPFMC (HCR 5), is refused there; SB0 is replaced by MSSB0
          // below.
          SB0(sp, yr) +=  NByage0(sp, 0, age, yr) *  spawn_output(sp, age, nyrs_hind - 1) * exp(-M_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0);
          SBF(sp, yr) +=  NByageF(sp, 0, age, yr) *  spawn_output(sp, age, nyrs_hind - 1) * exp(-(M_at_age(sp, 0, age, yr) + Ftarget_at_age(sp, 0, age, yr)) * spawn_month(sp)/12.0);
          DynamicSB0(sp, yr) +=  N_at_age_dB0(sp, 0, age, yr) *  spawn_output(sp, age, yr) * exp(-M_at_age_dB0(sp, 0, age, yr) * spawn_month(sp)/12.0);
          DynamicSBF(sp, yr) +=  N_at_age_dBF(sp, 0, age, yr) *  spawn_output(sp, age, yr) * exp(-(M_at_age_dBF(sp, 0, age, yr) + Ftarget_at_age(sp, 0, age, yr)) * spawn_month(sp)/12.0);

          for(sex = 0; sex < nsex(sp); sex ++){

            wt_idx_pop = 2 * sp;
            B0(sp, yr) +=  NByage0(sp, sex, age, yr) *  weight_hat( wt_idx_pop, sex, age, nyrs_hind - 1 );
            DynamicB0(sp, yr) +=  N_at_age_dB0(sp, sex,age, yr) *  weight_hat( wt_idx_pop, sex, age, yr );
          }
        }


        // Input SB0 (if running in multi-species mode)
        if(msmMode > 0){
          SB0(sp, yr) = MSSB0(sp);
          B0(sp, yr) = MSB0(sp);
        }
      }
    }

    // 6.7-6.9. FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (mt), and ssb-AT-AGE (mt)
    // Includes Harvest Control Rules
    for(sp = 0; sp < nspp; sp++) {
      for(yr = nyrs_hind; yr < nyrs; yr++){

        // 6.7. HARVEST CONTROL RULES FOR PROJECTION
        // Equilibrium or dynamic reference points (i.e. SB0 or dynamic SB0)
        Type ref_SB0 = (DynamicHCR == 1) ? DynamicSB0(sp, yr-1) : SB0(sp, nyrs-1);
        Type ref_SBF = (DynamicHCR == 1) ? DynamicSBF(sp, yr-1) : SBF(sp, nyrs-1);

        switch(HCR){
        case 0: // No fishing
          proj_F(sp, yr) = 0.0;
          break;

        case 1: // CMSY
          // proj_F(sp, yr) retains its incoming value
          break;

        case 2: // Constant F
          // proj_F(sp, yr) retains its incoming value
          break;

        case 3: // Constant F to acheive X% of SB0
          // proj_F(sp, yr) retains its incoming value
          break;

        case 4: // Constant Fspr with multiplier
          // proj_F(sp, yr) retains its incoming value
          break;

        case 5: // NPFMC Tier 3 HCR
          if(ssb(sp, yr-1) < ref_SBF){
            proj_F(sp, yr) = Ftarget(sp) * (((ssb(sp, yr-1)/ref_SBF)-Alpha(sp))/(1-Alpha(sp))); // Used Fabc of FtargetSPR%
          }
          if((ssb(sp, nyrs_hind-1) < ref_SB0 * Plimit(sp)) || (ssb(sp, yr-1) / ref_SBF < Alpha(sp))){ // If overfished
            proj_F(sp, yr) = 0.0;
          }
          break;

        case 6: // PFMC Category 1 HCR
          if(ssb(sp, yr-1) < ref_SB0 * Ptarget(sp)){
            proj_F(sp, yr) = (Flimit(sp) + QnormHCR(sp)) * (ref_SB0 * Ptarget(sp) * (ssb(sp, yr-1) - ref_SB0 * Plimit(sp))) / (ssb(sp, yr-1) * (ref_SB0 * (Ptarget(sp) - Plimit(sp))));
          }
          if(ssb(sp, yr-1) < ref_SB0 * Plimit(sp)){ // If overfished
            proj_F(sp, yr) = 0.0;
          }
          break;

        case 7: // SESSF Tier 1 HCR
          if(ssb(sp, yr-1) < ref_SB0 * Ptarget(sp)){
            proj_F(sp, yr) = Ftarget(sp) * ((ssb(sp, yr-1)/(ref_SB0 * Plimit(sp)))-1); // Used Fabc of FtargetSPR%
          }
          if(ssb(sp, yr-1) < ref_SB0 * Plimit(sp)){ // If overfished
            proj_F(sp, yr) = 0.0;
          }
          break;
        }


        // No projected F when the forecast is off or numbers-at-age are input (estDynamics > 0).
        if((forecast(sp) == 0) | (estDynamics(sp) > 0)){
          proj_F(sp, yr) =  0.0;
        }

        // -- Update F for the projection (account for selectivity and fleets)
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            F_at_age(sp, sex, age, yr) = 0.0;
          }
        }

        // -- Multiply F from HCR by selectivity and fleet proportion
        F_spp(sp, yr) = proj_F(sp, yr);
        for(flt = 0; flt < n_flt; flt++) {
          if((sp == flt_spp(flt)) & (flt_type(flt) == 1)){ // Fisheries only, as in 5.12
            F_flt(flt, yr) = proj_F_prop(flt) * proj_F(sp, yr);
            for(age = 0; age < nages(sp); age++) {
              for(sex = 0; sex < nsex(sp); sex ++){
                F_flt_age(flt, sex, age, yr) = sel_at_age(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * proj_F(sp, yr); // Projection F uses the terminal hindcast year's selectivity
                if(flt_type(flt) == 1){
                  F_at_age(sp, sex, age, yr) += F_flt_age(flt, sex, age, yr);
                }
              }
            }
          }
        }


        // -- Update mortality for forecast with new F
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            M_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
            Z_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + F_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
          }
        }


        // 6.8. FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (mt), and ssb-AT-AGE (mt)
        // -- 6.8.1. Forecasted recruitment
        // - Option 1: Use mean rec
        if(proj_mean_rec == 1){
          R(sp, yr) = avg_R(sp); //  // Projections use mean R given bias in R0
        }

        // - Option 2: Use SRR and rec devs
        if(proj_mean_rec == 0){

          // - Calculate recruitment (with linkage offsets pre-added).
          // Pre-styr spawning year takes R_init, as in the hindcast block;
          // reachable only if the hindcast is shorter than minage.
          int proj_spawn_yr = yr - minage(sp);
          int proj_srr_use = (proj_spawn_yr < 0) ? 0 : srr_pred_fun;
          Type ssb_tmp = (proj_spawn_yr < 0) ? Type(0.0) : ssb(sp, proj_spawn_yr);
          Type proj_rec_mean = (proj_spawn_yr < 0) ? R_init(sp) : R0(sp, yr);
          // TODO: projected recruitment is not floored under an identity-link offset (issue #234).
          R(sp, yr) = calculate_recruitment(proj_srr_use, proj_rec_mean, ssb_tmp, alpha(sp, yr), Beta(sp, yr), rec_dev(sp, yr), SPR0(sp));
        }

        N_at_age(sp, 0, 0 , yr) = R(sp, yr) * sex_ratio(sp, 0);
        if(nsex(sp) > 1){
          N_at_age(sp, 1, 0 , yr) = R(sp, yr) * (1 - sex_ratio(sp, 0));
        }


        // -- 6.8.2. Ages > recruitment
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            switch(estDynamics(sp)){
            case 0: // Estimated numbers-at-age

              // -- 6.8.2.  Where minage + 1 <= age < Ai
              if(age < (nages(sp) - 1)) {
                N_at_age(sp, sex, age + 1, yr) = N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1));
              }

              // -- 6.8.3. Plus group where age > Ai.
              if(age == (nages(sp) - 1)) {
                N_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age - 1, yr - 1) * exp(-Z_at_age(sp, sex, age-1, yr-1)) + N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1));
              }
              break;

            case 1: // Fixed numbers-at-age - fixed scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              break;

            case 2: // Fixed numbers-at-age age-independent scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              break;

            default: // Wrong estDynamics
              error("Invalid 'estDynamics'");
            }

            // -- 6.8.4. FORECAST BIOMASS
            wt_idx_pop = 2 * sp;
            biomass_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * weight_hat( wt_idx_pop, sex, age, nyrs_hind-1 ); // 6.5.
            biomass(sp, yr) += biomass_at_age(sp, sex, age, yr);
          } // End sex loop

          // -- 6.8.5. FORECAST SSB (SUM ACROSS AGES)
          wt_idx_ssb = 2 * sp + 1;
          ssb(sp, yr) += N_at_age(sp, 0, age, yr) * exp(-Z_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0) * spawn_output(sp, age, nyrs_hind-1); // 6.6.
        }
      }
    }


    // 6.9. EXPECTED RECRUITMENT
    // The function-scope `yr` is left at `nyrs` by the forecast loop above, so
    // the year-0 block uses the explicit `first_yr` constant: reading
    // R0/alpha/Beta at `yr` would run one column past the end and segfault
    // inside MakeADFun. The assertion pins those matrices to the year horizon.
    const int first_yr = 0;
    if (R0.cols() != nyrs || alpha.cols() != nyrs || Beta.cols() != nyrs) {
      error("Section 6.9: R0/alpha/Beta column count does not equal nyrs; "
            "recruitment-parameter matrices were sized inconsistently.");
    }
    for(sp = 0; sp < nspp; sp++) {

      // Year 1 (arent fit in likelihood)
      // Under predation there is no SPRFinit (6.2) to place the first year on
      // the curve, so it takes R_init, the anchor its realised recruitment uses.
      if((msmMode > 0) & (srr_pred_fun > 1)){
        R_hat(sp, first_yr) = R_init(sp);
      } else {
        switch(srr_pred_fun){
        case 0: // Random about mean (e.g. Alaska)
          R_hat(sp, first_yr) = R0(sp, first_yr);
          break;

        case 1: // Random about mean with environmental effects
          R_hat(sp, first_yr) = R0(sp, first_yr);
          break;

        case 2: // Beverton-Holt
          R_hat(sp, first_yr) = (alpha(sp, first_yr) - 1/SPRFinit(sp)) / Beta(sp, first_yr); // (Alpha-1/SPR0)/beta
          break;

        case 3: // Retired alias of 2; build_srr() refuses it (alpha moves via a linkage)
          R_hat(sp, first_yr) = (alpha(sp, first_yr) - 1/SPRFinit(sp)) / Beta(sp, first_yr); // (Alpha-1/SPR0)/beta
          break;

        case 4: // Ricker
          { Type aS = alpha(sp, first_yr) * SPRFinit(sp);   // kept positive so the log is finite; R_hat itself is floored below
            if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
            R_hat(sp, first_yr) = log(aS) / (Beta(sp, first_yr) * SPRFinit(sp)/1000000.0); }
          break;

        case 5: // Ricker with environmental impacts on alpha
          { Type aS = alpha(sp, first_yr) * SPRFinit(sp);   // kept positive so the log is finite; R_hat itself is floored below
            if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
            R_hat(sp, first_yr) = log(aS) / (Beta(sp, first_yr) * SPRFinit(sp)/1000000.0); }
          break;
        default:
          error("Invalid 'srr_pred_fun'");
        }
      }
      // Same floor as R in 6.6.1: the first year's R_hat is reported and read by retrospective().
      if (rec_floor_on) { Type pen_R = 0; R_hat(sp, first_yr) = posfun(R_hat(sp, first_yr), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }

      // Year 1+
      for(yr = 1; yr < nyrs; yr++){
        // Pre-styr spawning year takes R_init, the SAME anchor the realised recruitment
        // above uses, so expected and realised differ only by rec_dev and the
        // stock-recruit penalty (13.9) scores the deviation alone -- these years have no
        // spawning biomass to carry information about the curve's level. Anchoring on the
        // curve instead (R_hat's own first year) leaves a mean-versus-curve level gap the
        // penalty reads as signal: under the Ianelli configuration (srr_fun mean,
        // srr_pred_fun BevertonHolt) a fixed 0.83 nats per guarded year, pulling rec_pars
        // against alpha and Beta from years with no data. Identical for
        // BevertonHolt/BevertonHolt and Ricker/Ricker, where R_init is already R_hat's
        // first-year value.
        int hat_spawn_yr = yr - minage(sp);
        int hat_srr_use = (hat_spawn_yr < 0) ? 0 : srr_pred_fun;
        Type ssb_tmp = (hat_spawn_yr < 0) ? Type(0.0) : ssb(sp, hat_spawn_yr);
        Type hat_rec_mean = (hat_spawn_yr < 0) ? R_init(sp) : R0(sp, yr);

        // Expected recruitment carries no deviation, hence Type(0.0).
        R_hat(sp, yr) = calculate_recruitment(hat_srr_use, hat_rec_mean, ssb_tmp, alpha(sp, yr), Beta(sp, yr), Type(0.0), SPR0(sp));
        // Same floor as R in 6.6.1: log(R_hat) enters the stock-recruit penalty.
        if (rec_floor_on) { Type pen_R = 0; R_hat(sp, yr) = posfun(R_hat(sp, yr), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }
      }
    }



    // 6.11. ESTIMATE AVERAGE NUMBERS AT AGE
    for(sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        for(age = 0; age < nages(sp); age++) {
          for(yr = 0; yr < nyrs; yr++) {
          // MSVPA average N over the year.
            avgN_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) / Z_at_age(sp, sex, age, yr);
            avgN_at_age_dB0(sp, sex, age, yr) = N_at_age_dB0(sp, sex, age, yr) * (1 - exp(-M_at_age_dB0(sp, sex, age, yr))) / M_at_age_dB0(sp, sex, age, yr);
            avgN_at_age_dBF(sp, sex, age, yr) = N_at_age_dBF(sp, sex, age, yr) * (1 - exp(-M_at_age_dBF(sp, sex, age, yr) - Ftarget_at_age(sp, sex, age, yr))) / (M_at_age_dBF(sp, sex, age, yr) + Ftarget_at_age(sp, sex, age, yr));
          }
        }
      }
    }


    /** ------------------------------------------------------------------------ //
     // 7. MSVPA PREDATION MORTALITY EQUATIONS                                   //
     * ------------------------------------------------------------------------- */

    // START PREDATION
    if(msmMode > 0) {

      // NOTE -- LOOPING INDICES -- sp = species, age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length

      // 7.1. Empirical MSVPA Suitability
      calculate_msvpa_suitability(stom_div_bio, suitability, suit_other,
                                  nspp, nyrs, suit_styr, suit_endyr, nyrs_suit, msmMode,
                                  nsex, nages, suitMode, diet_prop, avgN_at_age, weight_hat, other_food_diet_prop);


      // 7.2. MSVPA PREDATION MORTALITY
      if((msmMode == 1) | (msmMode == 2)) {
        calculate_msvpa_predation(avail_food, avail_food_dB0, avail_food_dBF,
                                  M2_at_age, M2_at_age_dB0, M2_at_age_dBF, M2_prop,
                                  B_eaten, B_eaten_as_prey, diet_prop_hat, other_diet_prop_hat,
                                  nspp, nyrs, nsex, nages, msmMode,
                                  avgN_at_age, avgN_at_age_dB0, avgN_at_age_dBF,
                                  suitability, suit_other, weight_hat, consumption_at_age, other_food);
      }
    } // End 7. Predation mortality
  } // End population dynamics iterations



  /** ------------------------------------------------------------------------ //
   // 8. INDEX COMPONENTS EQUATIONS                                             //
   * ------------------------------------------------------------------------- */

  // -- 8.1. Index of abundance/biomass
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){

    index = index_ctl(index_ind, 0) - 1;          // Survey index
    sp = index_ctl(index_ind, 1) - 1;             // Species
    flt_yr = index_ctl(index_ind, 2);             // Data year
    wt_idx_flt = nspp * 2 + index;                // Fleet weight index
    mo = index_n(index_ind, 0);                   // Month

    index_hat(index_ind) = 0.0;

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;
    }
    else if(flt_yr < 0){
      flt_yr = -flt_yr - styr;
    }

    for(age = 0; age < nages(sp); age++) {
      for(sex = 0; sex < nsex(sp); sex++){

        // Index by fleet type. Month is not used for a fishery: a year-average
        // has no instant to be taken at. SS3 splits the same way on its
        // per-fleet survey timing (SS_expval.tpl: timing >= 0 -> exp(-Z*timing);
        // timing < 0, required of every fishing fleet -> (1 - exp(-Z))/Z).
        //
        // A seasonal fishery takes the same form. Its exact predictor over [t1, t1+D] is
        //     N_a * exp(-Zpre_a * t1) * (1 - exp(-Zwin_a * D)) / (Zwin_a * D)
        //     Zpre_a = M_a + other fleets' F_a      (this fishery has yet to fish)
        //     Zwin_a = M_a + F_a/D + other fleets' F_a
        // Since Zwin*D = M*D + F, the window factor tracks F almost as the annual
        // one does, and q absorbs the F-free exp(-M*t1) offset. Trend error
        // against the exact form is ~1.5% over F 0.05-0.8, against -29% to +33%
        // for the snapshot. TODO: fit the window form directly if a fleet needs it.
        Type n_index;
        if(flt_type(index) == 1){
          n_index = N_at_age(sp, sex, age, flt_yr) * (1.0 - exp( - Z_at_age(sp, sex, age, flt_yr))) / Z_at_age(sp, sex, age, flt_yr);
        } else {
          n_index = N_at_age(sp, sex, age, flt_yr) * exp( - (mo/12.0) * Z_at_age(sp, sex, age, flt_yr));
        }

        // Weight
        if(flt_units(index) == 1){
          index_hat(index_ind) += n_index * sel_at_age(index, sex, age, flt_yr) * weight_hat(wt_idx_flt, sex, age, flt_yr );
        }
        // Numbers
        if(flt_units(index) == 2){
          index_hat(index_ind) += n_index * sel_at_age(index, sex, age, flt_yr);
        }
      }
    }
  }


  // -- 8.2. Analytical survey q following Walters and Ludwig 1994
  index_n_obs.setZero();
  index_q_analytical.setZero();
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){



    index = index_ctl(index_ind, 0) - 1;            // Survey index
    sp = index_ctl(index_ind, 1) - 1;             // Species
    flt_yr = index_ctl(index_ind, 2);             // Data year

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;


      mo = index_n(index_ind, 0);                    // Month
      if(flt_yr < nyrs_hind){

        // If etimated standard deviation or analytical sigma (non-time-varying)
        if(est_sigma_index(index) > 0) {
          index_n_obs(index) += 1;
          index_q_analytical(index) += log(index_obs(index_ind, 0) / index_hat(index_ind));
        }

        // If time-varying sigma
        if(est_sigma_index(index) == 0 ) {
          index_n_obs(index) += 1 / square(index_obs(index_ind, 1));
          index_q_analytical(index) += log(index_obs(index_ind, 0) / index_hat(index_ind)) / square(index_obs(index_ind, 1));
        }
      }
    }
  }

  // Take average
  for(index = 0 ; index < n_flt; index ++){
    index_q_analytical(index) = exp(index_q_analytical(index) / index_n_obs(index));

    // Set index_q to analytical if used
    if(est_index_q(index) == 3){
      for(yr = 0; yr < nyrs_hind; yr++){
        index_q(index, yr) = index_q_analytical(index);
      }
    }
  }


  // -- 8.2b. Arithmetic-mean analytical q (AMAK/ebswp): q = sum(obs)/sum(pred)
  //    over the fitted survey years (est_index_q == 7, "AnalyticalArith"). This
  //    reproduces the AMAK BTS q_bts = mean(ob_bts)/mean(eb_bts) and is the q
  //    form paired with the MVN covariance survey likelihood (Index_loglike =
  //    "MVN"). Uses index_hat before it is scaled by q just below.
  vector<Type> index_obs_sum(n_flt); index_obs_sum.setZero();
  vector<Type> index_hat_sum(n_flt); index_hat_sum.setZero();
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){
    index = index_ctl(index_ind, 0) - 1;
    flt_yr = index_ctl(index_ind, 2);
    if(flt_yr > 0){
      flt_yr = flt_yr - styr;
      if(flt_yr < nyrs_hind){
        index_obs_sum(index) += index_obs(index_ind, 0);
        index_hat_sum(index) += index_hat(index_ind);
      }
    }
  }
  for(index = 0; index < n_flt; index++){
    if(est_index_q(index) == 7){
      Type q_arith = index_obs_sum(index) / index_hat_sum(index);
      for(yr = 0; yr < nyrs_hind; yr++){
        index_q(index, yr) = q_arith;
      }
    }
  }


  // -- 8.3. Survey Biomass - multiply by q
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){

    index = index_ctl(index_ind, 0) - 1;            // Survey index
    flt_yr = index_ctl(index_ind, 2);      // Data year

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;
    }
    else if(flt_yr < 0){
      flt_yr = -flt_yr - styr;
    }

    // Hindcast or projection indexing for q
    int yr_ind = (flt_yr < nyrs_hind) ? flt_yr : (nyrs_hind - 1);

    index_hat(index_ind) = index_q(index, yr_ind) * index_hat(index_ind); // pow(index_hat(index_ind), (1 + index_q_pow(index)));
  }



  // --8.4. Calculate analytical sigma following Walters and Ludwig 1994
  // The positive-observation guard matches the catch estimator below; it is
  // belt-and-braces, data_check() having already refused a non-positive index.
  index_n_obs.setZero();
  index_analytical_sd.setZero();
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){

    index = index_ctl(index_ind, 0) - 1;            // Survey index
    flt_yr = index_ctl(index_ind, 2);      // Data year

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;

      if((flt_yr < nyrs_hind) && (index_obs(index_ind, 0) > 0)){
        index_n_obs(index) += 1;
        index_analytical_sd(index) += square( log(index_obs(index_ind, 0)) - log(index_hat(index_ind)));
      }
    }
  }

  for(index = 0 ; index < n_flt; index ++){
    if(index_n_obs(index) > 0){
      index_analytical_sd(index) = concentrated_lognormal_sd(index_analytical_sd(index) / index_n_obs(index), bias_adjust_obs);
    }
  }


  /** ------------------------------------------------------------------------ //
   // 9. FISHERY COMPONENTS EQUATIONS                                          //
   * ------------------------------------------------------------------------- */
  // 9.1. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (mt or thousands of fish)
  for(fsh_ind = 0; fsh_ind < catch_ctl.rows(); fsh_ind++){

    flt = catch_ctl(fsh_ind, 0) - 1;     // Fleet index (0-based)
    sp = catch_ctl(fsh_ind, 1) - 1;      // Species
    flt_yr = catch_ctl(fsh_ind, 2);      // Data year
    mo = catch_n(fsh_ind, 0);            // Month
    wt_idx_flt = nspp * 2 + flt;         // Weight index for fleet

    catch_hat(fsh_ind) = 0.0;
    max_catch_hat(fsh_ind) = 0.0;

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;
    }
    else if(flt_yr < 0){
      flt_yr = -flt_yr - styr;
    }

    for(sex = 0; sex < nsex(sp); sex ++){
      for(age = 0; age < nages(sp); age++) {

        // By weight
        if(flt_units(flt) == 1){
          catch_hat(fsh_ind) += F_flt_age(flt, sex, age, flt_yr) / Z_at_age(sp, sex, age, flt_yr) * (1.0 - exp(-Z_at_age(sp, sex, age, flt_yr))) * N_at_age(sp, sex, age, flt_yr) * weight_hat( wt_idx_flt, sex, age, flt_yr ); // 5.5.
          max_catch_hat(fsh_ind) += N_at_age(sp, sex, age, flt_yr) * weight_hat( wt_idx_flt, sex, age, flt_yr ) * sel_at_age(flt, sex, age, flt_yr) * proj_F_prop(flt);
        }

        // By numbers
        if(flt_units(flt) == 2){
          catch_hat(fsh_ind) += F_flt_age(flt, sex, age, flt_yr) / Z_at_age(sp, sex, age, flt_yr) * (1.0 - exp(-Z_at_age(sp, sex, age, flt_yr))) * N_at_age(sp, sex, age, flt_yr);
          max_catch_hat(fsh_ind) += N_at_age(sp, sex, age, flt_yr) * sel_at_age(flt, sex, age, flt_yr) * proj_F_prop(flt);
        }
      }
    }
  }

  // -- 9.1a. Initial equilibrium catch, in weight or numbers per the fleet's
  //    Observation_units, as catch_hat is:
  //   C_eq = sum_a Finit * s_a * w_a * N_eq_a * (1 - exp(-Z_a)) / Z_a
  //   Z_a  = M1_a + Finit * s_a
  // Baranov on the deviation-free equilibrium age structure, at the fleet's own
  // selectivity and weight in the first hindcast year; SS3's Equil_catch
  // (SS_popdyn.tpl, Do_Equil_Calc). R supplies these rows only under initMode 6,
  // the one mode that builds N_eq at this same Finit * selectivity.
  vector<Type> equil_catch_hat(equil_catch_ctl.rows()); equil_catch_hat.setZero();
  for(int eq_ind = 0; eq_ind < equil_catch_ctl.rows(); eq_ind++){
    flt = equil_catch_ctl(eq_ind, 0) - 1;
    sp  = equil_catch_ctl(eq_ind, 1) - 1;
    wt_idx_flt = nspp * 2 + flt;
    for(sex = 0; sex < nsex(sp); sex++){
      for(age = 0; age < nages(sp); age++){
        Type s_a = sel_at_age(flt, sex, age, 0);
        Type Z_a = M1_at_age(sp, sex, age, 0) + Finit(sp) * s_a;
        Type baranov = Finit(sp) * s_a / Z_a * (1.0 - exp(-Z_a)) * N_eq(sp, sex, age);
        if(flt_units(flt) == 1){          // by weight
          equil_catch_hat(eq_ind) += baranov * weight_hat( wt_idx_flt, sex, age, 0 );
        }
        if(flt_units(flt) == 2){          // by numbers
          equil_catch_hat(eq_ind) += baranov;
        }
      }
    }
  }

  // -- 9.1b. Analytical catch sigma, following Walters and Ludwig (1994)
  // The sd that minimises the catch density: squared log residuals accumulated
  // here, concentrated in concentrated_lognormal_sd(). The row predicate matches
  // the likelihood's own -- a fitted year, a fishery, a positive observation --
  // so the sd is estimated over exactly the rows it is applied to. A zero-catch
  // year is legal input and log(0) is -inf, which would poison the fleet's sd,
  // hence the positive-observation guard.
  //
  // It reads catch_obs, not obsvec, so oneStepPredict()'s candidate value never
  // reaches it: the sd enters an OSA sweep as a fixed full-data number, as an
  // estimated sd held at its MLE would.
  catch_n_obs.setZero();
  catch_analytical_sd.setZero();
  for(fsh_ind = 0; fsh_ind < catch_ctl.rows(); fsh_ind++){

    flt = catch_ctl(fsh_ind, 0) - 1;
    flt_yr = catch_ctl(fsh_ind, 2);

    if((flt_yr > 0) && (flt_yr <= endyr) && (flt_type(flt) == 1) && (catch_obs(fsh_ind, 0) > 0)){
      catch_n_obs(flt) += 1;
      catch_analytical_sd(flt) += square( log(catch_obs(fsh_ind, 0)) - log(catch_hat(fsh_ind)) );
    }
  }

  for(flt = 0; flt < n_flt; flt++){
    if(catch_n_obs(flt) > 0){
      catch_analytical_sd(flt) = concentrated_lognormal_sd(catch_analytical_sd(flt) / catch_n_obs(flt), bias_adjust_obs);
    }
  }


  // 9.2 Exploitable biomass ----
  exploitable_biomass.setZero();

  for(flt = 0; flt < n_flt; flt++) {

    sp = flt_spp(flt);
    wt_idx_flt = nspp * 2 + flt;

    if(flt_type(flt) == 1){ //Fishery only
      for(yr = 0; yr < nyrs; yr++) {
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            exploitable_biomass(sp, yr) += N_at_age(sp, sex, age, yr) * weight_hat( wt_idx_flt, sex, age, yr ) * sel_at_age(flt, sex, age, yr) * proj_F_prop(flt);
          }
        }
      }
    }
  }


  /** ------------------------------------------------------------------------ //
   * 10. COMPOSITION EQUATIONS                                                  //
   * ------------------------------------------------------------------------- */

  // -- 10.1. Marginal age/length composition
  age_hat.setZero();     // Age-composition
  age_obs_hat.setZero(); // Age-composition with ageing error
  comp_hat.setZero();    // Age- or length-comp following data
  for(comp_ind = 0; comp_ind < comp_hat.rows(); comp_ind++){

    flt = comp_ctl(comp_ind, 0) - 1;            // Fleet index (0-based)
    sp = comp_ctl(comp_ind, 1) - 1;             // Species
    flt_sex = comp_ctl(comp_ind, 2);            // Comp sex: 0 = combined, 1 = female, 2 = male
    comp_type = comp_ctl(comp_ind, 3);          // Comp type: 0 = age, 1 = length, 2 = CAAL
    yr = comp_ctl(comp_ind, 4);                 // Data year
    mo = (growth_model(sp) > 0) ? flt_month(flt) : comp_n(comp_ind, 0); // Month (use fleet month if estimating growth, otherwise use from data)
    int wtind = nspp * 2 + flt;
    n_hat(comp_ind) = 0.0;

    // Likelihood (yr > 0) vs prediction (yr < 0)
    if(yr > 0){
      yr = yr - styr;
    }
    else if(yr < 0){
      yr = -yr - styr;
    }

    // Hindcast
    Type Frate = 0.0;
    if(yr < nyrs_hind){
      yr_ind = yr;
      Frate = exp(log_F(flt, yr));
    }

    // Projection
    if(yr >= nyrs_hind){
      yr_ind = nyrs_hind - 1;
      Frate = proj_F_prop(flt) * proj_F(sp, yr);
    }

    // 10.1.1. Estimated growth
    if(growth_model(sp) > 0){
      for(age = 0; age < nages(sp); age++) {
        for(ln = 0; ln < nlengths(sp); ln++) {
          for(sex = 0; sex < nsex(sp); sex ++){

            switch(flt_type(flt)){
            case 1: // - Fishery
              if(flt_sel_dim(flt) == 1){ // Length based
                {
                  // sel(L) * P(L|age) over this data bin's population bins (2.3c).
                  Type ps_ln = 0.0;
                  for(int lp = pop_bin_lo(sp, ln); lp <= pop_bin_hi(sp, ln); lp++){
                    ps_ln += sel_at_length(flt, sex, lp, yr) * growth_matrix(wtind, sex, age, lp, yr);
                  }
                  pred_CAAL(flt, sex, age, ln, yr) = ps_ln * (Frate / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr));
                }
              }
              break;


            case 2: // - Survey
              if(flt_sel_dim(flt) == 1){ // Length based
                {
                  // sel(L) * P(L|age) over this data bin's population bins (2.3c).
                  Type ps_ln = 0.0;
                  for(int lp = pop_bin_lo(sp, ln); lp <= pop_bin_hi(sp, ln); lp++){
                    ps_ln += sel_at_length(flt, sex, lp, yr) * growth_matrix(wtind, sex, age, lp, yr);
                  }
                  pred_CAAL(flt, sex, age, ln, yr) = ps_ln * (N_at_age(sp, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr)));
                }
              }
              break;
            }
            // Comps read this, not pred_CAAL, which CAAL reads and stays length-based. A fleet with
            // selectivity at age (Fixed is always at age) has no sel(L), so sel(a) weights P(L|age).
            comp_pred(sex, age, ln) = pred_CAAL(flt, sex, age, ln, yr);
            if(flt_sel_dim(flt) != 1 || flt_sel_type(flt) == 0){
              Type p_ln = 0.0;
              for(int lp = pop_bin_lo(sp, ln); lp <= pop_bin_hi(sp, ln); lp++){
                p_ln += growth_matrix(wtind, sex, age, lp, yr);
              }
              switch(flt_type(flt)){
              case 1: // - Fishery
                comp_pred(sex, age, ln) = p_ln * (F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr));
                break;

              case 2: // - Survey
                comp_pred(sex, age, ln) = p_ln * (N_at_age(sp, sex, age, yr) * sel_at_age(flt, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr)));
                break;
              }
            }
          }
        }
      }

      // Assign CAAL to marginal age or length composition
      for(age = 0; age < nages(sp); age++) {
        for(ln = 0; ln < nlengths(sp); ln++) {

          // Age composition
          // - Catch-at-age
          if(comp_type == 0){
            switch(flt_sex){
            case 0: // Sexes combined or 1 sex assessment
              for(sex = 0; sex < nsex(sp); sex ++){
                age_hat(comp_ind, age ) += comp_pred(sex, age, ln);
              }
              break;

            case 1: case 2: // Sex-specific composition data
              sex = flt_sex - 1;
              age_hat(comp_ind, age ) += comp_pred(sex, age, ln);
              break;

            case 3: // Joint composition data
              for(sex = 0; sex < nsex(sp); sex ++){
                age_hat(comp_ind, age + nages(sp) * sex ) += comp_pred(sex, age, ln);
              }
              break;
            }
          }


          // Length composition
          // - Catch-at-length
          if(comp_type == 1){
            switch(flt_sex){
            case 0: // Sexes combined or 1 sex assessment
              for(sex = 0; sex < nsex(sp); sex ++){
                comp_hat(comp_ind, ln ) += comp_pred(sex, age, ln);
              }
              break;

            case 1: case 2: // Sex-specific composition data
              sex = flt_sex - 1;
              comp_hat(comp_ind, ln ) += comp_pred(sex, age, ln);
              break;

            case 3: // Joint composition data
              for(sex = 0; sex < nsex(sp); sex ++){
                comp_hat(comp_ind, ln + nlengths(sp) * sex ) += comp_pred(sex, age, ln) ;
              }
              break;
            }
          }
        }
      }

      // Adjust for aging error
      for(int obs_age = 0; obs_age < nages(sp); obs_age++) {
        for(int true_age = 0; true_age < nages(sp); true_age++) {
          age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(flt_ageing_error_index(flt), true_age, obs_age);
        }
      }

      // Adjust for aging error for joint data
      if(flt_sex == 3){
        for(int obs_age = nages(sp); obs_age < nages(sp) * 2; obs_age++) {
          for(int true_age = nages(sp); true_age < nages(sp) * 2; true_age++) {

            // Adjust indexing for joint age/length comp
            int true_age_tmp = true_age - nages(sp);
            int obs_age_tmp = obs_age - nages(sp);

            age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(flt_ageing_error_index(flt), true_age_tmp, obs_age_tmp);
          }
        }
      }

      // Assign age-comp to estimated comp
      if( comp_type == 0) {
        joint_adjust(comp_ind) = 1;
        if(flt_sex == 3){
          joint_adjust(comp_ind) = 2;
        }
        for(age = 0; age < nages(sp) * joint_adjust(comp_ind); age++) {
          comp_hat(comp_ind, age ) = age_obs_hat(comp_ind, age );
        }
      }
    } // End estimated growth loop


    // 10.1.2. Empirical weight-at-age
    if(growth_model(sp) == 0){
      for(age = 0; age < nages(sp); age++) {

        switch(flt_type(flt)){
        case 1: // - Fishery
          switch(flt_sex){
          case 0: // Sexes combined or 1 sex assessment
            for(sex = 0; sex < nsex(sp); sex ++){
              // Catch-at-age
              age_hat(comp_ind, age ) += F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            }
            break;

          case 1: case 2: // Sex-specific composition data
            sex = flt_sex - 1;
            // Catch-at-age
            age_hat(comp_ind, age ) = F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            break;

          case 3: // Joint composition data
            for(sex = 0; sex < nsex(sp); sex ++){
              // Catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            }
            break;
          }
          break;


        case 2: // - Survey
          switch(flt_sex){
          case 0: // Sexes combined or 1 sex assessment
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age ) += N_at_age(sp, sex, age, yr)  * sel_at_age(flt, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
            }
            break;

          case 1: case 2: // Sex-specific composition data
            sex = flt_sex - 1;
            // Survey catch-at-age
            age_hat(comp_ind, age ) = N_at_age(sp, sex, age, yr)  * sel_at_age(flt, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
            break;

          case 3: // Joint composition data
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = N_at_age(sp, sex, age, yr)  * sel_at_age(flt, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
            }
            break;
          }
        }
      }


      // Adjust for aging error
      for(int obs_age = 0; obs_age < nages(sp); obs_age++) {
        for(int true_age = 0; true_age < nages(sp); true_age++) {
          age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(flt_ageing_error_index(flt), true_age, obs_age);
        }
      }

      // Adjust for aging error for joint data
      if(flt_sex == 3){
        for(int obs_age = nages(sp); obs_age < nages(sp) * 2; obs_age++) {
          for(int true_age = nages(sp); true_age < nages(sp) * 2; true_age++) {

            // Adjust indexing for joint age/length comp
            int true_age_tmp = true_age - nages(sp);
            int obs_age_tmp = obs_age - nages(sp);

            age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(flt_ageing_error_index(flt), true_age_tmp, obs_age_tmp);
          }
        }
      }

      // Assign age-comp to estimated comp
      if( comp_type == 0) {
        joint_adjust(comp_ind) = 1;
        if(flt_sex == 3){
          joint_adjust(comp_ind) = 2;
        }
        for(age = 0; age < nages(sp) * joint_adjust(comp_ind); age++) {
          comp_hat(comp_ind, age ) = age_obs_hat(comp_ind, age );
        }
      }

      // Convert from catch-at-age to catch-at-length, reading age_hat, the TRUE ages, not
      // age_obs_hat: the age-length key is P(length | true age), and a measured length
      // carries no otolith reading, so ageing error has no part in a length composition.
      // SS3 likewise applies its age_age matrix only to age and CAAL data.
      if( comp_type == 1) {
        for(ln = 0; ln < nlengths(sp); ln++) {
          for(age = 0; age < nages(sp); age++) {

            sex = 0;

            // Adjust sex for males/females
            if((flt_sex > 0) & (flt_sex < 3)){
              sex = flt_sex - 1;
            }
            comp_hat(comp_ind, ln ) += age_hat(comp_ind, age) * age_trans_matrix(flt_age_transition_index(flt), sex, age, ln );
          }
        }


        // Convert from catch-at-age to catch-at-length for joint comp data
        if( flt_sex == 3) {
          for(ln = nlengths(sp); ln < nlengths(sp) * 2; ln++) {
            for(age = nages(sp); age < nages(sp) * 2; age++) {

              // Adjust indexing for joint age/length comp
              int obs_log_tmp = ln;
              int obs_age_tmp = age;

              obs_age_tmp = age - nages(sp);
              obs_log_tmp = ln - nlengths(sp);
              sex = 1;

              comp_hat(comp_ind, ln ) += age_hat(comp_ind, age) * age_trans_matrix(flt_age_transition_index(flt), sex, obs_age_tmp, obs_log_tmp );
            }
          }
        }
      }
    } // End empirical weight-at-age loop

    // Standardize to sum to 1
    Type penalty = 0.0;
    comp_hat.row(comp_ind) /= posfun(comp_hat.row(comp_ind).sum(), Type(0.001), penalty);
    (void) penalty;
  }


  // -- 10.2. CAAL
  // -- length-based selectivity was converted above
  caal_hat.setZero();
  for(int caal_ind = 0; caal_ind < caal_hat.rows(); caal_ind++){

    flt = caal_ctl(caal_ind, 0) - 1;            // Fleet index (0-based)
    sp = caal_ctl(caal_ind, 1) - 1;             // Species
    flt_sex = caal_ctl(caal_ind, 2);            // CAAL sex: 0 = combined, 1 = female, 2 = male
    yr = caal_ctl(caal_ind, 3);                 // Data year
    ln = caal_ctl(caal_ind, 4) - 1;             // Length bin
    mo = flt_month(flt);                        // Month (use fleet month b/c estimating growth)

    if(growth_model(sp) == 0) continue; // Skip empirical weight-at-age species

    if(flt_sex > 0){
      sex = flt_sex - 1;
    } else {
      sex = flt_sex;
    }

    int wtind = nspp * 2 + flt;

    // Likelihood (yr > 0) vs prediction (yr < 0)
    if(yr > 0){
      yr = yr - styr;
    }
    else if(yr < 0){
      yr = -yr - styr;
    }
    // Hindcast
    Type Frate = 0.0;
    if(yr < nyrs_hind){
      yr_ind = yr;
      Frate = exp(log_F(flt, yr));
    }

    // Projection
    if(yr >= nyrs_hind){
      yr_ind = nyrs_hind - 1;
      Frate = proj_F_prop(flt) * proj_F(sp, yr);
    }

    // 10.2. Calculate CAAL (predicted conditional age-at-length)
    for(age = 0; age < nages(sp); age++) {

      switch(flt_type(flt)){
      case 1: // - Fishery
        if(flt_sel_dim(flt) == 1){
          {
                  // sel(L) * P(L|age) over this data bin's population bins (2.3c).
                  Type ps_ln = 0.0;
                  for(int lp = pop_bin_lo(sp, ln); lp <= pop_bin_hi(sp, ln); lp++){
                    ps_ln += sel_at_length(flt, sex, lp, yr) * growth_matrix(wtind, sex, age, lp, yr);
                  }
                  pred_CAAL(flt, sex, age, ln, yr) = ps_ln * (Frate / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr));
                }
        }
        break;


      case 2: // - Survey
        if(flt_sel_dim(flt) == 1){
          {
                  // sel(L) * P(L|age) over this data bin's population bins (2.3c).
                  Type ps_ln = 0.0;
                  for(int lp = pop_bin_lo(sp, ln); lp <= pop_bin_hi(sp, ln); lp++){
                    ps_ln += sel_at_length(flt, sex, lp, yr) * growth_matrix(wtind, sex, age, lp, yr);
                  }
                  pred_CAAL(flt, sex, age, ln, yr) = ps_ln * (N_at_age(sp, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr)) );
                }
        }
        break;
      }
    }

    // Adjust for aging error
    for(int obs_age = 0; obs_age < nages(sp); obs_age++) {
      for(int true_age = 0; true_age < nages(sp); true_age++) {
        caal_hat(caal_ind, obs_age) += pred_CAAL(flt, sex, true_age, ln, yr) * age_error(flt_ageing_error_index(flt), true_age, obs_age);
      }
    }

    //  Observed CAAL standardize to sum to 1
    Type penalty = 0.0;
    caal_hat.row(caal_ind) /= posfun(caal_hat.row(caal_ind).sum(), Type(0.0001), penalty);
    (void) penalty;
  }



  /** ------------------------------------------------------------------------ //
   * 11. PREDICTED STOMACH CONTENT                                             //
   * ------------------------------------------------------------------------- */

  // Predict stomach content
  if((msmMode > 2) || (imax(suitMode) > 0)) {
    predict_stomach_content(nspp, nyrs_hind, styr, suit_styr, suit_endyr,
                            minage, nsex, nages, diet_obs, diet_ctl,
                            avgN_at_age, diet_prop_hat, diet_hat);
  }


  /** ------------------------------------------------------------------------ //
   // 12. DERIVED QUANTITIES
   * ------------------------------------------------------------------------- */


  // 12.1. Depletion
  biomass_depletion.setZero();
  ssb_depletion.setZero();
  for(sp = 0; sp < nspp; sp++){
    for(yr = 0; yr < nyrs; yr++){

      if(DynamicHCR == 0){
        biomass_depletion(sp, yr) = biomass(sp, yr)/B0(sp, nyrs-1);
        ssb_depletion(sp, yr) = ssb(sp, yr)/SB0(sp, nyrs-1);
      }

      if(DynamicHCR == 1){
        biomass_depletion(sp, yr) = biomass(sp, yr)/DynamicB0(sp, yr);
        ssb_depletion(sp, yr) = ssb(sp, yr)/DynamicSB0(sp, yr);
      }

      // Multi-species, no HCR: the projection runs unfished, so SSB in the last
      // projection year is multispecies SB0 once `projyr` is far enough out to
      // equilibrate. Excluded under DynamicHCR, where DynamicSB0 is the
      // reference and this would overwrite it.
      if((HCR == 0) & (msmMode > 0) & (DynamicHCR == 0)){
        biomass_depletion(sp, yr) = biomass(sp, yr)/biomass(sp, nyrs-1);
        ssb_depletion(sp, yr) = ssb(sp, yr)/ssb(sp, nyrs-1);
      }
    }
  }




  /** ------------------------------------------------------------------------ //
   // 13. LIKELIHOOD EQUATIONS
   * ------------------------------------------------------------------------- */

  // 13.0. OBJECTIVE FUNCTION
  // Data likelihoods come first (13.1-13.6), then priors and penalties (13.7-13.14).
  // Each SIMULATE block sits beside the density it draws from, in the same order,
  // so sim_mod() draws what the likelihood scores.
  // dnorm(..., true) and similar return a log density, which is subtracted; TMB's
  // MVNORM, AR1, SEPARABLE and GMRF return a negative log density, which is added.
  // Data terms read observations from obsvec, gated by keep so oneStepPredict() can
  // hold each one out for OSA residuals; in ordinary fitting keep = 1.
  int n_col = std::max(n_flt, nspp);
  // Rows of jnll_comp. Keep in sync with the row names in R/6-rename_output.R and
  // .JNLL_ROW_AXIS in R/9-profile.R; test-schema-jnll-rows.R checks all three.
  enum JnllRow {
    JNLL_INDEX          = 0,   // Index data
    JNLL_CATCH          = 1,   // Catch data
    JNLL_COMP           = 2,   // Composition data
    JNLL_CAAL           = 3,   // CAAL data
    JNLL_SEL_NONPARAM   = 4,   // Non-parametric selectivity
    JNLL_SEL_DEV        = 5,   // Selectivity deviates
    JNLL_Q_PRIOR        = 6,   // Catchability prior
    JNLL_Q_DEV          = 7,   // Catchability deviates
    JNLL_SRR_PRIOR      = 8,   // Stock-recruit prior
    JNLL_INIT_DEV       = 9,   // Initial abundance deviates
    JNLL_REC_DEV        = 10,  // Recruitment deviates
    JNLL_SRR_PENALTY    = 11,  // Stock-recruit penalty
    JNLL_REFPT_PENALTY  = 12,  // Reference point penalties
    JNLL_ZERO_N_PENALTY = 13,  // Zero n-at-age penalty
    JNLL_M_PRIOR        = 14,  // M prior
    JNLL_M_RE           = 15,  // M random effects
    JNLL_RATION         = 16,  // Ration
    JNLL_RATION_PENALTY = 17,  // Ration penalties
    JNLL_STOMACH        = 18,  // Stomach content data
    JNLL_LINKAGE_PRIOR  = 19,  // Linkage-table priors (per-row)
    JNLL_LINKAGE_RE     = 20,  // Linkage random effects
    JNLL_EQUIL_CATCH    = 21,  // Initial equilibrium catch
    JNLL_N_ROWS         = 22   // total row count (for dimensioning)
  };
  matrix<Type> jnll_comp(JNLL_N_ROWS, n_col); jnll_comp.setZero();  // negative log-likelihood components
  matrix<Type> unweighted_jnll_comp(JNLL_N_ROWS, n_col); unweighted_jnll_comp.setZero();  // see 13.15

  // A column is a fleet on rows JNLL_INDEX-JNLL_Q_DEV and JNLL_EQUIL_CATCH, a species
  // on JNLL_SRR_PRIOR-JNLL_LINKAGE_PRIOR, and neither on JNLL_LINKAGE_RE (one
  // model-wide density), so never sum down a column. Nothing writes the two ration rows.


  // 13.1. INDEX OF ABUNDANCE (JNLL_INDEX)
  // Index_distribution: 0 Lognormal; 1 MVN and 2 MVNORM, correlated across years
  // (fitted per fleet below); 3 Normal and 4 TruncatedNormal, both with an absolute
  // sd on the natural scale.
  Type index_std_dev = 0;
  for(index_ind = 0; index_ind < index_obs.rows(); index_ind++) {

    index = index_ctl(index_ind, 0) - 1;          // Fleet
    sp = index_ctl(index_ind, 1) - 1;             // Species
    flt_yr = index_ctl(index_ind, 2);             // Year
    int family = index_ll_type(index);

    switch (est_sigma_index(index)) {
    case 0: index_std_dev = index_obs(index_ind, 1); break;          // input, per row
    case 1: index_std_dev = exp(index_log_sd(index)); break;         // estimated
    case 2: index_std_dev = index_analytical_sd(index); break;       // Walters and Ludwig (1994), section 8.4
    default: error("Invalid 'Estimate_sigma_index'");
    }
    index_sd(index_ind) = index_std_dev;

    // Fit hindcast years with an observation.
    bool fit_row = (flt_yr > 0) && (flt_yr <= endyr) && (flt_type(index) > 0) && (index_obs(index_ind, 0) > 0);

    // Lognormal, read from obsvec (= log index) so OSA can gate each row with keep.
    if(fit_row && (family == 0)){
      int pos = index_obsvec_idx(index_ind);
      if(pos >= 0){
        Type index_mu = log(index_hat(index_ind)) - bias_adjust_obs*square(index_std_dev)/2.0;
        jnll_comp(JNLL_INDEX, index) -= keep(pos) * dnorm(obsvec(pos), index_mu, index_std_dev, true);
        if(osa_mode == 2){
          jnll_comp(JNLL_INDEX, index) -= osa_norm_cdf_terms(
            obsvec(pos), index_mu, index_std_dev, keep.cdf_lower(pos), keep.cdf_upper(pos));
        }
      }
    }

    // Normal (3) matches AMAK's avo_like/cpue_like, 0.5 * (obs - q * pred)^2 / sd^2.
    // TruncatedNormal (4) renormalizes over obs > 0: log Phi(mu / sd), in [log 0.5, 0].
    if(fit_row && ((family == 3) || (family == 4))){
      Type log_Z = (family == 4) ? log(pnorm(index_hat(index_ind) / index_std_dev)) : Type(0.0);
      if(osa_mode == 0){
        jnll_comp(JNLL_INDEX, index) -= dnorm(index_obs(index_ind, 0), index_hat(index_ind), index_std_dev, true);
        jnll_comp(JNLL_INDEX, index) += log_Z;
      } else {
        // OSA: obsvec holds the untransformed index for these families.
        int pos = index_obsvec_idx(index_ind);
        if(pos >= 0){
          jnll_comp(JNLL_INDEX, index) -= keep(pos) * dnorm(obsvec(pos), index_hat(index_ind), index_std_dev, true);
          jnll_comp(JNLL_INDEX, index) += keep(pos) * log_Z;
          if(osa_mode == 2){
            if(family == 4){
              jnll_comp(JNLL_INDEX, index) -= osa_trunc_norm_cdf_terms(
                obsvec(pos), index_hat(index_ind), index_std_dev,
                keep.cdf_lower(pos), keep.cdf_upper(pos));
            } else {
              jnll_comp(JNLL_INDEX, index) -= osa_norm_cdf_terms(
                obsvec(pos), index_hat(index_ind), index_std_dev,
                keep.cdf_lower(pos), keep.cdf_upper(pos));
            }
          }
        }
      }
    }

    // -- Simulate (families 0, 3, 4; MVN is drawn per fleet below). Projection rows
    // are always drawn, as run_mse() reveals them; a missing hindcast year stays NA.
    SIMULATE {
      bool sim_row = (flt_yr <= 0) || (flt_yr > endyr) || (index_obs(index_ind, 0) > 0);
      if(sim_row && (flt_type(index) > 0) && (family == 0)){
        index_obs(index_ind, 0) = exp(rnorm(log(index_hat(index_ind)) - bias_adjust_obs*square(index_std_dev)/2.0, index_std_dev));
        int sim_pos = index_obsvec_idx(index_ind);
        if(sim_pos >= 0){
          obsvec(sim_pos) = log(index_obs(index_ind, 0));
        }
      }
      if(sim_row && (flt_type(index) > 0) && (family == 3)){
        // The untruncated normal can draw <= 0, which data_check() refuses, so redraw and
        // count rejections; sim_mod() warns when they matter. Family 4 has no such gap.
        Type sim_draw = rnorm(index_hat(index_ind), index_std_dev);
        int  tn_tries = 0;
        while((sim_draw <= Type(0)) && (tn_tries < index_trunc_max_tries)){
          index_trunc_rejects_sim(index_ind) += 1;
          sim_draw = rnorm(index_hat(index_ind), index_std_dev);
          tn_tries++;
        }
        index_trunc_tries_sim(index_ind) += tn_tries + 1;
        index_trunc_sd_sim(index_ind) = index_std_dev;
        index_trunc_budget_sim(index_ind) = index_trunc_max_tries;
        index_obs(index_ind, 0) = sim_draw;
      }
      if(sim_row && (flt_type(index) > 0) && (family == 4)){
        // Exact inverse-CDF draw from the normal restricted to (0, inf).
        Type a = pnorm(-index_hat(index_ind) / index_std_dev);
        Type u = a + (Type(1.0) - a) * runif(Type(0.0), Type(1.0));
        index_obs(index_ind, 0) = index_hat(index_ind) + index_std_dev * qnorm(u);
      }
      // Keep obsvec in step: log(obs) when fitting, the untransformed obs under OSA.
      if((flt_type(index) > 0) && ((family == 3) || (family == 4))){
        int sim_pos = index_obsvec_idx(index_ind);
        if(sim_pos >= 0){
          if(osa_mode == 0){
            if(index_obs(index_ind, 0) > Type(0)) obsvec(sim_pos) = log(index_obs(index_ind, 0));
          } else {
            obsvec(sim_pos) = index_obs(index_ind, 0);
          }
        }
      }
    }
  }

  // MVN (1) and MVNORM (2): the AMAK/ebswp DoCovBTS survey likelihood on the arithmetic
  // residual r = obs - q * pred, with the fleet's fixed covariance index_cov_mat. MVN
  // reports the bare quadratic form 0.5 r' Sigma^-1 r, MVNORM the full density; the fit
  // is the same. Pair with Catchability 7 for the AMAK arithmetic-mean q.
  for(index = 0; index < n_flt; index++){
    if((flt_type(index) > 0) && (index_ll_type(index) == 1 || index_ll_type(index) == 2)){
      int n_mvn = 0;
      for(index_ind = 0; index_ind < index_obs.rows(); index_ind++){
        if((index_ctl(index_ind, 0) - 1) == index){
          flt_yr = index_ctl(index_ind, 2);
          if((flt_yr > 0) && (flt_yr <= endyr) && (index_obs(index_ind, 0) > 0)) n_mvn++;
        }
      }
      if(n_mvn > 0){
        // Fitted rows in index_obs order, which is the row order of Sigma.
        vector<Type> resid(n_mvn);
        vector<Type> mu(n_mvn);
        vector<int>  posv(n_mvn);
        vector<int>  rowv(n_mvn);
        int k = 0;
        for(index_ind = 0; index_ind < index_obs.rows(); index_ind++){
          if((index_ctl(index_ind, 0) - 1) == index){
            flt_yr = index_ctl(index_ind, 2);
            if((flt_yr > 0) && (flt_yr <= endyr) && (index_obs(index_ind, 0) > 0)){
              resid(k) = index_obs(index_ind, 0) - index_hat(index_ind);
              mu(k)    = index_hat(index_ind);
              posv(k)  = index_obsvec_idx(index_ind);
              rowv(k)  = index_ind;
              k++;
            }
          }
        }
        if(osa_mode == 0){
          // MVNORM() returns 0.5 * (r' Sigma^-1 r + logdet(Sigma) + n log(2 pi));
          // index_cov_const holds the constant MVN drops.
          Type dens = MVNORM(index_cov_mat(index))(resid);
          jnll_comp(JNLL_INDEX, index) += (index_ll_type(index) == 2) ? dens : (dens - index_cov_const(index));
        } else {
          // OSA: build_osa_data() stores z = L^-1 obs, with Sigma = L L'. Each z_k is an
          // independent N((L^-1 mu)_k, 1) (Thygesen et al. 2017).
          Eigen::LLT<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > llt(index_cov_mat(index));
          matrix<Type> Lmat = llt.matrixL();
          matrix<Type> muc(n_mvn, 1);
          for(k = 0; k < n_mvn; k++) muc(k, 0) = mu(k);
          matrix<Type> mwhite = Lmat.template triangularView<Eigen::Lower>().solve(muc);  // L^-1 mu
          for(k = 0; k < n_mvn; k++){
            int pos = posv(k);
            if(pos >= 0){
              jnll_comp(JNLL_INDEX, index) -= keep(pos) * dnorm(obsvec(pos), mwhite(k, 0), Type(1.0), true);
              if(osa_mode == 2){
                jnll_comp(JNLL_INDEX, index) -= osa_norm_cdf_terms(
                  obsvec(pos), mwhite(k, 0), Type(1.0),
                  keep.cdf_lower(pos), keep.cdf_upper(pos));
              }
            }
          }
        }

        // -- Simulate the fitted rows only (Sigma covers no others), redrawing the vector
        // while any row is <= 0: exact for the positive-orthant truncation, which per-margin
        // truncation would decorrelate. A draw left when the budget runs out is reported.
        SIMULATE {
          vector<Type> sim_dev = MVNORM(index_cov_mat(index)).simulate();
          int mvn_tries = 0;
          while(true){
            bool any_nonpos = false;
            for(k = 0; k < n_mvn; k++){
              if(mu(k) + sim_dev(k) <= Type(0)) any_nonpos = true;
            }
            if(!any_nonpos) break;
            if(mvn_tries >= index_trunc_max_tries * n_mvn) break;
            for(k = 0; k < n_mvn; k++){
              if(mu(k) + sim_dev(k) <= Type(0)) index_trunc_rejects_sim(rowv(k)) += 1;
            }
            sim_dev = MVNORM(index_cov_mat(index)).simulate();
            mvn_tries++;
          }
          for(k = 0; k < n_mvn; k++){
            index_trunc_tries_sim(rowv(k)) += mvn_tries + 1;
            index_trunc_sd_sim(rowv(k)) = sqrt(index_cov_mat(index)(k, k));
            index_trunc_budget_sim(rowv(k)) = index_trunc_max_tries * n_mvn;
            index_obs(rowv(k), 0) = mu(k) + sim_dev(k);
            // obsvec is log(obs) when fitting and whitened under OSA; update only the former.
            int sim_pos = posv(k);
            if((sim_pos >= 0) && (osa_mode == 0) && (index_obs(rowv(k), 0) > Type(0))){
              obsvec(sim_pos) = log(index_obs(rowv(k), 0));
            }
          }
        }
      }
    }
  }

  // Reported under a _sim name -- see the naming rule in section 5.2b.
  SIMULATE {
    matrix<Type> index_obs_sim = index_obs;
    REPORT(index_obs_sim);
    REPORT(index_trunc_tries_sim);
    REPORT(index_trunc_rejects_sim);
    REPORT(index_trunc_sd_sim);
    REPORT(index_trunc_budget_sim);
  }


  // 13.2. CATCH (JNLL_CATCH)
  // Lognormal on total catch, fishery fleets and hindcast years only.
  Type fsh_std_dev = 0;
  for(fsh_ind = 0; fsh_ind < catch_obs.rows(); fsh_ind++) {

    flt = catch_ctl(fsh_ind, 0) - 1;            // Fleet
    sp = catch_ctl(fsh_ind, 1) - 1;             // Species
    flt_yr = catch_ctl(fsh_ind, 2);             // Year
    yr = flt_yr - styr;                         // Year index, 0-based

    switch (est_sigma_fsh(flt)) {
    case 0: fsh_std_dev = catch_obs(fsh_ind, 1); break;          // input, per row
    case 1: fsh_std_dev = exp(catch_log_sd(flt)); break;         // estimated
    case 2: fsh_std_dev = catch_analytical_sd(flt); break;       // Walters and Ludwig (1994), section 9.1b
    default: error("Invalid 'Estimate_sigma_catch'");
    }
    catch_sd(fsh_ind) = fsh_std_dev;

    // Read from obsvec (= log catch) so OSA can gate each row with keep.
    if((flt_yr > 0) && (flt_yr <= endyr) && (flt_type(flt) == 1) && (catch_obs(fsh_ind, 0) > 0)){
      int pos = catch_obsvec_idx(fsh_ind);
      if(pos >= 0){
        Type catch_mu = log(catch_hat(fsh_ind)) - bias_adjust_obs*square(fsh_std_dev)/2.0;
        jnll_comp(JNLL_CATCH, flt) -= keep(pos) * dnorm(obsvec(pos), catch_mu, fsh_std_dev, true) ;
        if(osa_mode == 2){
          jnll_comp(JNLL_CATCH, flt) -= osa_norm_cdf_terms(
            obsvec(pos), catch_mu, fsh_std_dev, keep.cdf_lower(pos), keep.cdf_upper(pos));
        }
      }
    }

    // -- Simulate every row, with the dnorm's bias-corrected mean. Projection rows
    // (NA catch from clean_data()) are filled too; where catch_hat = 0 the draw is
    // exp(-inf) = 0, safe only on doubles. Do not share this mean with the dnorm (log 0
    // puts a NaN on the tape) or add a catch_hat > 0 guard (it changes the draw count).
    SIMULATE {
      catch_obs(fsh_ind, 0) = exp(rnorm(log(catch_hat(fsh_ind)) - bias_adjust_obs*square(fsh_std_dev)/2.0, fsh_std_dev));
      int sim_pos = catch_obsvec_idx(fsh_ind);
      if(sim_pos >= 0){
        obsvec(sim_pos) = log(catch_obs(fsh_ind, 0));
      }
    }
  }


  // 13.3. INITIAL EQUILIBRIUM CATCH (JNLL_EQUIL_CATCH)
  // Lognormal on its own Log_sd and its own row, as SS3's equ_catch_like. No OSA: one
  // observation has no sequence to condition on.
  for(int eq_ind = 0; eq_ind < equil_catch_ctl.rows(); eq_ind++){
    int eq_flt = equil_catch_ctl(eq_ind, 0) - 1;
    if((equil_catch_obs(eq_ind, 0) > 0) && (flt_type(eq_flt) == 1)){
      Type eq_sd = equil_catch_obs(eq_ind, 1);
      Type eq_mu = log(equil_catch_hat(eq_ind)) - bias_adjust_obs*square(eq_sd)/2.0;
      jnll_comp(JNLL_EQUIL_CATCH, eq_flt) -=
        dnorm(log(equil_catch_obs(eq_ind, 0)), eq_mu, eq_sd, true);

      // -- Simulate. Under initMode 6 this is the only observation on Finit, so a
      // self_test() replicate must redraw it.
      SIMULATE {
        equil_catch_obs(eq_ind, 0) = exp(rnorm(eq_mu, eq_sd));
      }
    }
  }

  // Reported under _sim names -- see the naming rule in section 5.2b. obsvec_sim is for
  // inspection: only its catch and index entries are redrawn.
  SIMULATE {
    matrix<Type> catch_obs_sim = catch_obs;
    matrix<Type> equil_catch_obs_sim = equil_catch_obs;
    vector<Type> obsvec_sim = obsvec;
    REPORT(catch_obs_sim);
    REPORT(equil_catch_obs_sim);
    REPORT(obsvec_sim);
  }


  // 13.4. AGE AND LENGTH COMPOSITION (JNLL_COMP)
  // Likelihood types are in comp_row_nll() (likelihood.hpp). comp_offset (default
  // 1e-5, from fit_control()) keeps log(0) out, and build_osa_data() uses the same.
  // FIXME: the multinomial (type 0) is renormalized after the offset, so its reported
  // NLL carries a per-row constant N_eff * log(1 + n_bins * offset); fit and MLE unchanged.
  Type comp_prop_offset = comp_offset;
  matrix<Type> comp_sim = comp_obs;   // simulated comps; see the _sim rule in section 5.2b

  for(comp_ind = 0; comp_ind < comp_obs.rows(); comp_ind++) {

    flt = comp_ctl(comp_ind, 0) - 1;        // Fleet
    sp = comp_ctl(comp_ind, 1) - 1;         // Species
    flt_sex = comp_ctl(comp_ind, 2);        // 0 = combined, 1 = female, 2 = male, 3 = joint
    comp_type = comp_ctl(comp_ind, 3);      // 0 = age, 1 = length
    yr = comp_ctl(comp_ind, 4);             // Year

    // A joint-sex row holds both sexes, so twice the bins.
    joint_adjust(comp_ind) = 1;
    if(flt_sex == 3){
      joint_adjust(comp_ind) = 2;
    }
    int nbins_blk = (comp_type == 0) ? nages(sp) : nlengths(sp);
    int n_comp = 0;   // a comp type other than age or length scores no bins
    if((comp_type == 0) || (comp_type == 1)) n_comp = nbins_blk * joint_adjust(comp_ind);

    vector<Type> comp_obs_tmp = comp_obs.row(comp_ind).segment(0, n_comp); // Observed proportion
    vector<Type> comp_hat_tmp = comp_hat.row(comp_ind).segment(0, n_comp); // Expected proportion

    // Tail accumulation (Comp_accum_young / Comp_accum_old). OSA reads an obsvec that
    // build_osa_data() folds the same way.
    fold_comp_tails(comp_obs_tmp, comp_hat_tmp, n_comp, nbins_blk, joint_adjust(comp_ind),
                    comp_accum_young(flt), comp_accum_old(flt));

    if((yr <= endyr) && (yr > 0) && (flt_type(flt) > 0) && (comp_n(comp_ind, 1) > 0)){
      comp_row_nll(comp_ll_type(flt), osa_mode, comp_obs_tmp, comp_hat_tmp,
                   comp_prop_offset, comp_n(comp_ind, 1), comp_weights(flt), DM_pars_comp(flt),
                   obsvec, keep, comp_obsvec_idx(comp_ind), "comp_ll_type",
                   jnll_comp(JNLL_COMP, flt), unweighted_jnll_comp(JNLL_COMP, flt));
    }

    // -- Simulate in raw bins, before the fold and offset (both families are closed
    // under merging bins), projection rows included, as run_mse() reveals them.
    SIMULATE {
      int nbin_raw = ((comp_type == 0) ? nages(sp) : nlengths(sp)) * joint_adjust(comp_ind);
      vector<Type> sim_p = comp_hat.row(comp_ind).segment(0, nbin_raw);
      Type sim_tot = sim_p.sum();
      Type n_nom = comp_n(comp_ind, 1);
      // A row with nothing predicted or nothing sampled comes back EMPTY, never as the
      // noise-free prediction: run_mse() reads an all-zero row as "no catch, no comp".
      vector<Type> sim_out(nbin_raw); sim_out.setZero();

      if((sim_tot > Type(0)) && (n_nom > Type(0))){
        vector<Type> p = sim_p / sim_tot;
        if(comp_ll_type(flt) == 1){
          // DM: theta already carries the weight, so draw at the nominal N.
          vector<Type> sim_alphas = p * n_nom * DM_pars_comp(flt);
          sim_out = rdirmultinom_rce(n_nom, sim_alphas);
        } else {
          // Multinomial: the weight scales the log-likelihood, so it is an effective
          // sample size and the draw uses it too.
          sim_out = rmultinom_rce(n_nom * comp_weights(flt), p);
        }
      }
      for(int b = 0; b < nbin_raw; b++){
        comp_sim(comp_ind, b) = sim_out(b);
      }
      // Zero the padding past this row's bins, which rearrange_data() would otherwise
      // normalize into the real ones.
      for(int b = nbin_raw; b < comp_sim.cols(); b++){
        comp_sim(comp_ind, b) = Type(0);
      }
    }
  }

  SIMULATE {
    matrix<Type> comp_obs_sim = comp_sim;
    REPORT(comp_obs_sim);
  }


  // 13.5. CONDITIONAL AGE-AT-LENGTH (JNLL_CAAL)
  // One row per fleet, year, sex and length bin, scored like a composition row.
  matrix<Type> caal_sim = caal_obs;   // simulated CAAL; see the _sim rule in section 5.2b

  for(int caal_ind = 0; caal_ind < caal_hat.rows(); caal_ind++){

    flt = caal_ctl(caal_ind, 0) - 1;            // Fleet
    sp = caal_ctl(caal_ind, 1) - 1;             // Species
    flt_sex = caal_ctl(caal_ind, 2);            // 0 = combined, 1 = female, 2 = male, 3 = joint
    yr = caal_ctl(caal_ind, 3);                 // Year
    ln = caal_ctl(caal_ind, 4) - 1;             // Length bin

    int n_caal = nages(sp);

    vector<Type> caal_obs_tmp = caal_obs.row(caal_ind).segment(0, n_caal); // Observed proportion
    vector<Type> caal_hat_tmp = caal_hat.row(caal_ind).segment(0, n_caal); // Expected proportion

    if((yr <= endyr) && (yr > 0) && (flt_type(flt) > 0) && (caal_n(caal_ind, 0) > 0)){
      comp_row_nll(caal_ll_type(flt), osa_mode, caal_obs_tmp, caal_hat_tmp,
                   comp_prop_offset, caal_n(caal_ind, 0), caal_weights(flt), DM_pars_caal(flt),
                   obsvec, keep, caal_obsvec_idx(caal_ind), "caal_ll_type",
                   jnll_comp(JNLL_CAAL, flt), unweighted_jnll_comp(JNLL_CAAL, flt));
    }

    // -- Simulate, as for the compositions above (CAAL has no tail fold).
    SIMULATE {
      vector<Type> sim_p = caal_hat.row(caal_ind).segment(0, n_caal);
      Type sim_tot = sim_p.sum();
      Type n_nom = caal_n(caal_ind, 0);
      vector<Type> sim_out(n_caal); sim_out.setZero();

      if((sim_tot > Type(0)) && (n_nom > Type(0))){
        vector<Type> p = sim_p / sim_tot;
        if(caal_ll_type(flt) == 1){
          vector<Type> sim_alphas = p * n_nom * DM_pars_caal(flt);
          sim_out = rdirmultinom_rce(n_nom, sim_alphas);
        } else {
          sim_out = rmultinom_rce(n_nom * caal_weights(flt), p);
        }
      }
      for(int b = 0; b < n_caal; b++){
        caal_sim(caal_ind, b) = sim_out(b);
      }
      for(int b = n_caal; b < caal_sim.cols(); b++){
        caal_sim(caal_ind, b) = Type(0);
      }
    }
  }

  SIMULATE {
    matrix<Type> caal_obs_sim = caal_sim;
    REPORT(caal_obs_sim);
  }


  // 13.6. STOMACH CONTENTS (JNLL_STOMACH)
  // Each stomach's prey proportions plus an "other prey" balance, scored per predator
  // when suitability is estimated.
  matrix<Type> diet_sim = diet_obs;   // simulated diet; see the _sim rule in section 5.2b

  if((msmMode > 2) || (imax(suitMode) > 0)) {

    int current_j = 0; // Track position in diet_ctl

    for (int i = 0; i < n_stomach_obs; ++i) {

      // This stomach's prey are the run of diet_ctl rows with stomach_id == i. The scan
      // only moves forward, so ids must be sorted 0, 1, 2, ... with no gaps; clean_data()
      // builds them that way and data_check() refuses a table that is not.
      int start_j = current_j;
      while((current_j < diet_ctl.rows()) && (stomach_id(current_j) == i)) {
        current_j++;
      }
      int n_prey = current_j - start_j;
      if (n_prey == 0) continue;   // before reading row start_j, which may be past the end

      int rsp = diet_ctl(start_j, 0) - 1;   // predator species
      Type N_s = diet_obs(start_j, 0);      // stomach sample size

      if (suitMode(rsp) <= 0) continue;

      // Prey proportions, with "other prey" as the last bin.
      vector<Type> obs_diet_prop(n_prey + 1); obs_diet_prop.setZero();
      vector<Type> pred_diet_prop(n_prey + 1); pred_diet_prop.setZero();

      for (int k = 0; k < n_prey; ++k) {
        obs_diet_prop(k) = diet_obs(start_j + k, 1);
        pred_diet_prop(k) = diet_hat(start_j + k, 1);
      }

      Type sum_obs_p = obs_diet_prop.head(n_prey).sum();
      if (sum_obs_p > 1.0) sum_obs_p = 1.0;
      obs_diet_prop(n_prey) = 1.0 - sum_obs_p;

      Type sum_est_p = pred_diet_prop.head(n_prey).sum();
      pred_diet_prop(n_prey) = posfun(1.0 - sum_est_p, Type(0.00001), penalty); // smooth floor at 0

      // Same offset as the compositions, then renormalize.
      obs_diet_prop += comp_offset;
      pred_diet_prop += comp_offset;

      obs_diet_prop /= obs_diet_prop.sum();
      pred_diet_prop /= pred_diet_prop.sum();

      vector<Type> obs_diet_content = obs_diet_prop * N_s;   // observed counts

      // Not comp_row_nll(): the stomach multinomial is TMB's dmultinom() on proportions
      // renormalized after the offset, a different number from the composition form.

      Type stomach_log_likelihood = 0;
      Type unweighted_stomach_log_likelihood = 0;
      Type DM_diet_par = DM_diet_pars(rsp);

      vector<Type> diet_alphas = pred_diet_prop * N_s * DM_diet_par;
      vector<Type> unweighted_diet_alphas = pred_diet_prop * N_s;

      if(osa_mode == 0){
        switch(diet_ll_type(rsp)){

        case 0:  // Multinomial; the weight multiplies the log-likelihood
          stomach_log_likelihood = dmultinom(obs_diet_content, pred_diet_prop, true);

          unweighted_jnll_comp(JNLL_STOMACH, rsp) -= stomach_log_likelihood;
          jnll_comp(JNLL_STOMACH, rsp) -= diet_comp_weights(rsp) * stomach_log_likelihood;
          break;
        case 1:  // Dirichlet-multinomial; theta = exp(weight)
          stomach_log_likelihood = ddirmultinom(obs_diet_content, diet_alphas, true);
          unweighted_stomach_log_likelihood = ddirmultinom(obs_diet_content, unweighted_diet_alphas, true);

          unweighted_jnll_comp(JNLL_STOMACH, rsp) -= unweighted_stomach_log_likelihood;
          jnll_comp(JNLL_STOMACH, rsp) -= stomach_log_likelihood;
          break;

        default:
          error("Invalid 'diet_ll_type'");
        }
      } else {
        // OSA: unweighted and keep-gated, read from obsvec over the n_prey + 1 bins.
        int start = diet_obsvec_idx(i);
        if(start >= 0){
          vector<Type> osa_x = obsvec.segment(start, n_prey + 1);
          if(diet_ll_type(rsp) == 1){   // Dirichlet-multinomial: no CDF term (comp_osa.hpp)
            jnll_comp(JNLL_STOMACH, rsp) -= ddirmultinom_osa(osa_x, diet_alphas, keep.segment(start, n_prey + 1), 1, 1);
          } else {
            jnll_comp(JNLL_STOMACH, rsp) -= dmultinom_osa(osa_x, pred_diet_prop, keep.segment(start, n_prey + 1), 1, 1, osa_mode == 2);
          }
        }
      }

      // -- Simulate from the same family, from the prediction before the offset. Only
      // the prey bins are written back; "other prey" is recomputed on the next fit.
      SIMULATE {
        vector<Type> sim_p(n_prey + 1);
        Type prey_sum = 0.0;
        for(int k = 0; k < n_prey; k++){
          sim_p(k) = diet_hat(start_j + k, 1);
          prey_sum += sim_p(k);
        }
        sim_p(n_prey) = (prey_sum < Type(1.0)) ? Type(1.0) - prey_sum : Type(0.0);
        Type sim_tot = sim_p.sum();
        if(sim_tot > Type(0)) sim_p /= sim_tot;

        if((N_s > Type(0)) && (sim_tot > Type(0))){
          vector<Type> sim_counts(n_prey + 1);
          switch(diet_ll_type(rsp)){
          case 0: {
            // The weight on a multinomial log-likelihood is an effective sample size.
            Type n_eff = N_s * diet_comp_weights(rsp);
            sim_counts = rmultinom_rce(n_eff, sim_p);
            break;
          }
          case 1: {
            vector<Type> sim_alphas = sim_p * N_s * DM_diet_par;
            sim_counts = rdirmultinom_rce(N_s, sim_alphas);
            break;
          }
          default:
            error("Invalid 'diet_ll_type'");
          }
          // Divide by the count actually drawn (rmultinom_rce rounds N), so the prey
          // proportions cannot sum past 1. A stomach that draws nothing comes back empty.
          Type n_drawn = sim_counts.sum();
          for(int k = 0; k < n_prey; k++){
            diet_sim(start_j + k, 1) = (n_drawn > Type(0))
              ? sim_counts(k) / n_drawn
              : Type(0);
          }
        }
      }
    }
  }

  // Reported under a _sim name -- see the naming rule in section 5.2b.
  SIMULATE {
    matrix<Type> diet_obs_sim = diet_sim;
    REPORT(diet_obs_sim);
  }


  // 13.7. SELECTIVITY PENALTIES (JNLL_SEL_NONPARAM, JNLL_SEL_DEV)
  // Charged once per Selectivity_index, on the lead fleet. Bins are on the fleet's
  // selectivity dimension (ages, or population length bins).
  // FIXME: penalize every deviation, not a sub-range, and the bin and year range controls
  // (Sel_pen_first_bin, Sel_pen_last_bin, Sel_cap_bin, Sel_start_year) could go.
  for(flt = 0; flt < n_flt; flt++){
    sp = flt_spp(flt);
    bool sel_is_length = (flt_sel_dim(flt) == 1);
    int  nbins = sel_is_length ? nlengths_pop(sp) : nages(sp);

    if(flt_type(flt) > 0 && flt_sel_lead(flt) == 1){

      // NonParametric (2): Ianelli/AMAK penalties on log selectivity centred to mean 1
      // (arrowtooth ADMB model), each year's curve when it varies (Atka mackerel).
      // Refused under random_sel = TRUE; form 13 is the integrable version.
      if(flt_sel_type(flt) == 2) {
        int nyrs_tmp = 1;
        if((flt_varying_sel(flt) == 1) || (flt_varying_sel(flt) == 4)){
          nyrs_tmp = nyrs_hind;
        }

        for(yr = 0; yr < nyrs_tmp; yr++){

          // 1. Decreasing selectivity: (|d| + d) / 2 = max(d, 0), d > 0 if decreasing.
          // FIXME: AMAK starts this penalty at bin nbins/2; here it starts at bin 0.
          for(sex = 0; sex < nsex(sp); sex++){
            for(age = 0; age < (nbins - 1); age++) {
              Type sel_ratio_tmp = log_non_par_sel(flt, sex, age, yr) - log_non_par_sel(flt, sex, age + 1, yr);
              jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 0) * square( (CppAD::abs(sel_ratio_tmp) + sel_ratio_tmp)/2.0);
            }
          }

          // 2. Curvature: squared second differences.
          for(sex = 0; sex < nsex(sp); sex++){
            vector<Type> sel_tmp(nbins); sel_tmp.setZero();
            for(age = 0; age < nbins; age++) {
              sel_tmp(age) = log_non_par_sel(flt, sex, age, yr);
            }
            vector<Type> sel_d2 = first_difference( first_difference( sel_tmp ) );
            for(int a2 = 0; a2 < sel_d2.size(); a2++) {
              jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 1) * sel_d2(a2) * sel_d2(a2);
            }
          }

          // 3. Time variation. RandomWalk (4): the year-to-year change in realized log
          // selectivity (Atka mackerel). IID (1): each deviate, over only the bins
          // build_map() estimates, as a fixed bin would pull the sd toward zero.
          if(flt_varying_sel(flt) == 4){
            if(yr > 0){
              for(sex = 0; sex < nsex(sp); sex++){
                for(age = 0; age < (nbins - 1); age++) {
                  jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_non_par_sel(flt, sex, age, yr), log_non_par_sel(flt, sex, age, yr - 1), sel_dev_sd(flt), true);
                }
              }
            }
          } else if(flt_varying_sel(flt) == 1){
            for(sex = 0; sex < nsex(sp); sex++){
              for(int bin = bin_first_selected(flt); bin < flt_n_sel_bins(flt); bin++) {
                jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_coff_dev(flt, sex, bin, yr), Type(0.0), sel_dev_sd(flt), true);
              }
            }
          }

          // 4. Average selectivity: holds each year's mean, deviates included.
          for(sex = 0; sex < nsex(sp); sex++){
            jnll_comp(JNLL_SEL_NONPARAM, flt) += 2.0 * square(avg_sel(flt, sex, yr));
          }
        }
      }

      // NonParametricPM (9): ADMB's Selectivity_Likelihood. JNLL_SEL_NONPARAM is
      // sel_like(1): shape, Sel_curve_pen1 = ctrl_flag(13); plus the optional avgsel (5).
      // JNLL_SEL_DEV is sel_like_dev(1): curvature, Sel_curve_pen2 = ctrl_flag(11)/nch; a
      // random walk with no normalizing constant; increments, Sel_curve_pen3 =
      // ctrl_flag(10)/group_num.
      if(flt_sel_type(flt) == 9) {

        int nyrs_tmp = 1;
        if(flt_varying_sel(flt) == 4){
          nyrs_tmp = nyrs_hind;
        }

        // Years before the fleet's selectivity start year (0-based) are skipped, as in
        // ADMB, whose base-year terms are anchored at styr.
        int start_yr = flt_sel_start_yr(flt);
        int shape_a0 = flt_sel_pen_first_bin(flt);              // first left bin of the shape pairs
        if(shape_a0 < 0) shape_a0 = bin_first_selected(flt);
        int shape_a1 = flt_sel_pen_last_bin(flt);               // last left bin of the shape pairs
        if(shape_a1 < 0) shape_a1 = nbins - 2;
        // Pairs are (age, age + 1), so the last left bin is nbins - 2.
        if(shape_a1 > nbins - 2) shape_a1 = nbins - 2;
        int shape_mode = flt_sel_shape_mode(flt);               // 0 = directional, 1 = smooth
        for(sex = 0; sex < nsex(sp); sex++){
          for(yr = start_yr; yr < nyrs_tmp; yr++){

            // (1) Shape. Directional (ADMB/AMAK): Sel_curve_pen1 >= 0 penalizes
            //     decreasing pairs, < 0 increasing ones. Smooth (RTMB "rpm"): two-sided d^2.
            for(age = shape_a0; age <= shape_a1; age++) {
              Type d = log_non_par_sel(flt, sex, age, yr) - log_non_par_sel(flt, sex, age + 1, yr); // > 0 if decreasing
              if(shape_mode == 1)
                jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 0) * d * d;
              else if(sel_curve_pen(flt, 0) >= 0)
                jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 0)  * square( (CppAD::abs(d) + d)/2.0 );
              else
                jnll_comp(JNLL_SEL_NONPARAM, flt) += -sel_curve_pen(flt, 0) * square( (CppAD::abs(d) - d)/2.0 );
            }

            // (2) Curvature; the start year counts only for a fleet starting at styr.
            if((yr > start_yr) || (start_yr == 0)){
              vector<Type> ls(nbins); ls.setZero();
              for(age = 0; age < nbins; age++) ls(age) = log_non_par_sel(flt, sex, age, yr);
              vector<Type> d2 = first_difference( first_difference( ls ) );
              for(int a2 = 0; a2 < d2.size(); a2++) jnll_comp(JNLL_SEL_DEV, flt) += sel_curve_pen(flt, 1) * d2(a2) * d2(a2);
            }

            // (3) Random walk, a bare sum of squares.
            if(yr > start_yr){
              for(age = 0; age < nbins; age++) {
                Type dd = log_non_par_sel(flt, sex, age, yr) - log_non_par_sel(flt, sex, age, yr - 1);
                jnll_comp(JNLL_SEL_DEV, flt) += dd * dd / (2.0 * sel_dev_sd(flt) * sel_dev_sd(flt));
              }
            }
          }

          // (4) Increment size: sum of squared sel_coff_dev, the random-walk increments.
          for(int bin = 0; bin < flt_n_sel_bins(flt); bin++){
            for(yr = start_yr; yr < nyrs_tmp; yr++){
              jnll_comp(JNLL_SEL_DEV, flt) += sel_curve_pen(flt, 2) * sel_coff_dev(flt, sex, bin, yr) * sel_coff_dev(flt, sex, bin, yr);
            }
          }

          // (5) AMAK avgsel: Sel_avgsel_pen * log(mean(exp(base coefficients)))^2 over
          //     the estimated bins; AMAK uses 10. Off (0) by default.
          if(flt_sel_avgsel_pen(flt) > 0){
            int n_base = flt_n_sel_bins(flt) - bin_first_selected(flt);
            vector<Type> base(n_base);
            for(int bin = 0; bin < n_base; bin++){
              base(bin) = sel_coff(flt, sex, bin + bin_first_selected(flt));
            }
            Type avgsel = log_mean_exp(base);
            jnll_comp(JNLL_SEL_NONPARAM, flt) += flt_sel_avgsel_pen(flt) * avgsel * avgsel;
          }
        }
      }

      // NonParametricIntegrable (13): the NonParametric shape penalties, charged once on
      // the base curve sel_coff, and a proper normal density on the deviates, so
      // random_sel = TRUE integrates a complete density. With no deviates it equals form 2.
      if(flt_sel_type(flt) == 13) {
        int n_sel_bins = flt_n_sel_bins(flt);
        int start_yr   = flt_sel_start_yr(flt);
        for(sex = 0; sex < nsex(sp); sex++){

          // The base curve as selectivity.hpp builds a year with no deviates: the
          // estimated bins, the last held through the tail, centred to a mean exp() of 1.
          vector<Type> base(nbins);
          for(int bin = 0; bin < nbins; bin++){
            base(bin) = sel_coff(flt, sex, (bin < n_sel_bins ? bin : n_sel_bins - 1));
          }
          vector<Type> est(n_sel_bins);
          for(int bin = 0; bin < n_sel_bins; bin++) est(bin) = base(bin);
          Type avg_base = log_mean_exp(est);
          Type ctr = log_mean_exp(base);
          for(int bin = 0; bin < nbins; bin++) base(bin) -= ctr;

          // 1. Decreasing selectivity.
          for(int bin = 0; bin < (nbins - 1); bin++) {
            Type d = base(bin) - base(bin + 1);
            jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 0) * square( (CppAD::abs(d) + d)/2.0 );
          }
          // 2. Curvature.
          vector<Type> d2 = first_difference( first_difference( base ) );
          for(int a2 = 0; a2 < d2.size(); a2++) {
            jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 1) * d2(a2) * d2(a2);
          }
          // 4. Average selectivity of the estimated coefficients.
          jnll_comp(JNLL_SEL_NONPARAM, flt) += 2.0 * square(avg_base);

          // 3. Deviates, over the estimated bins only: IID in every hindcast year, or
          //    random-walk increments after the start year (whose increment is fixed at 0).
          bool scored = (flt_varying_sel(flt) == 1) || (flt_varying_sel(flt) == 4);
          int yr_lo = (flt_varying_sel(flt) == 4) ? start_yr + 1 : 0;
          if(scored){
            for(yr = yr_lo; yr < nyrs_hind; yr++){
              for(int bin = bin_first_selected(flt); bin < n_sel_bins; bin++) {
                jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_coff_dev(flt, sex, bin, yr), Type(0.0), sel_dev_sd(flt), true);
              }
            }
          }
        }
      }

      // Logistic-family selectivity deviates: logistic (1), double logistic (3),
      // descending logistic (4), double normal (8). Each limb's inflection has sd and log
      // slope 4 * sd; the ascending random walk swaps them, to match the GOA pollock
      // assessment. Double normal: slot 0 is the peak and ascending width, slot 1 the
      // right-floor logit and descending width.
      bool asc_limb  = (flt_sel_type(flt) == 1) || (flt_sel_type(flt) == 3) || (flt_sel_type(flt) == 8);
      bool desc_limb = (flt_sel_type(flt) == 3) || (flt_sel_type(flt) == 4) || (flt_sel_type(flt) == 8);

      // Time_varying_sel IID (1) or 2, scored the same way.
      if((flt_varying_sel(flt) == 1) || (flt_varying_sel(flt) == 2)){
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = 0; yr < nyrs_hind; yr++){
            if(asc_limb){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_inf_dev(0, flt, sex, yr), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_sel_slp_dev(0, flt, sex, yr), Type(0.0), 4 * sel_dev_sd(flt), true);
            }
            if(desc_limb){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_inf_dev(1, flt, sex, yr), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_sel_slp_dev(1, flt, sex, yr), Type(0.0), 4 * sel_dev_sd(flt), true);
            }
          }
        }
      }

      // Random walk on both limbs (4) or the ascending limb only (5). Under 5 the
      // descending deviates are not estimated, so scoring them would bias the sd down.
      if((flt_varying_sel(flt) == 4) || (flt_varying_sel(flt) == 5)){
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = 1; yr < nyrs_hind; yr++){
            if(asc_limb){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_sel_slp_dev(0, flt, sex, yr) - log_sel_slp_dev(0, flt, sex, yr-1), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_inf_dev(0, flt, sex, yr) - sel_inf_dev(0, flt, sex, yr-1), Type(0.0), 4 * sel_dev_sd(flt), true);
            }
            if((flt_varying_sel(flt) == 4) && desc_limb){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_inf_dev(1, flt, sex, yr) - sel_inf_dev(1, flt, sex, yr-1), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_sel_slp_dev(1, flt, sex, yr) - log_sel_slp_dev(1, flt, sex, yr-1), Type(0.0), sel_dev_sd(flt) * 4, true);
            }
          }
        }
      }

      // LogisticPM (11): the AMAK "pm" survey penalties (ctrl_flag(19) > 0), from the year
      // after the selectivity start year. Sel_curve_pen1 (ctrl_flag(26)) on the change in
      // log selectivity over [Sel_norm_bin, Sel_norm_bin_upper], and Sel_curve_pen3 (ADMB's
      // 8) on the change in the free age-1 deviate. Sel_curve_pen2 is unused.
      if(flt_sel_type(flt) == 11){
        int start_yr = flt_sel_start_yr(flt);
        int alo = sel_norm_bin1(flt); int ahi = sel_norm_bin2(flt);
        if(alo < 0) alo = bin_first_selected(flt);            // default: whole selected range
        if(ahi < 0) ahi = nbins - 1;
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = start_yr + 1; yr < nyrs_hind; yr++){
            for(age = alo; age <= ahi; age++){
              Type s_now  = sel_is_length ? sel_at_length(flt, sex, age, yr)     : sel_at_age(flt, sex, age, yr);
              Type s_prev = sel_is_length ? sel_at_length(flt, sex, age, yr - 1) : sel_at_age(flt, sex, age, yr - 1);
              Type d = log(s_now) - log(s_prev);
              jnll_comp(JNLL_SEL_DEV, flt) += sel_curve_pen(flt, 0) * d * d;
            }
            Type da1 = sel_inf_dev(1, flt, sex, yr) - sel_inf_dev(1, flt, sex, yr - 1);
            jnll_comp(JNLL_SEL_DEV, flt) += sel_curve_pen(flt, 2) * da1 * da1;
          }
        }
      }

      // Hake non-parametric (5): IID deviates (Taylor et al. 2014), over all N_sel_bins,
      // including any bin build_map() holds at 0.
      if(((flt_varying_sel(flt) == 1) || (flt_varying_sel(flt) == 2)) && (flt_sel_type(flt) == 5)){
        for(int bin = 0; bin < flt_n_sel_bins(flt); bin++){
          for(sex = 0; sex < nsex(sp); sex++){
            for(yr = 0; yr < nyrs_hind; yr++){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_coff_dev(flt, sex, bin, yr), Type(0.0), sel_dev_sd(flt), true);
            }
          }
        }
      }

      // 2D AR1 (6): separable AR1 over year and bin. rho_trans() maps Sel_curve_pen1 and 2
      // to the bin and year correlations; sel_dev_sd is the conditional sd.
      if(flt_sel_type(flt) == 6){
        int n_sel_bins = flt_n_sel_bins(flt);
        array<Type> tmp_AR2(n_sel_bins, nyrs_hind);

        for(sex = 0; sex < nsex(sp); sex ++){
          for (int bin = 0; bin < n_sel_bins; bin++) {
            for (int yr = 0; yr < nyrs_hind; yr++) {
              tmp_AR2(bin, yr) =  sel_coff_dev(flt, sex, bin, yr);
            }
          }

          Type rho_a = rho_trans(sel_curve_pen(flt, 0));
          Type rho_y = rho_trans(sel_curve_pen(flt, 1));

          Type Sigma_sig_sel = pow(pow(sel_dev_sd(flt),2) / ((1-pow(rho_y,2))*(1-pow(rho_a,2))),0.5);
          // SEPARABLE(f, g): f correlates the last array index, g the first; tmp_AR2 is
          // (bin, year), so the year correlation goes first.
          jnll_comp(JNLL_SEL_DEV, flt) += SCALE(SEPARABLE(AR1(rho_y),AR1(rho_a)), Sigma_sig_sel)(tmp_AR2);
        }
      }

      // 3D AR1 (7): GMRF over bin, year and cohort (Cheng et al. 2024). Sel_curve_pen1-3
      // carry the bin, year and cohort correlations.
      if(flt_sel_type(flt) == 7){
        int n_sel_bins = flt_n_sel_bins(flt);
        int num_rows = n_sel_bins * nyrs_hind;
        matrix<Type> ay_index(num_rows, 2);

        for(sex = 0; sex < nsex(sp); sex ++){

          int row = 0;
          for (int year = 1; year <= nyrs_hind; ++year) {
            for (int bin = 1; bin <= n_sel_bins; ++bin) {
              ay_index(row, 0) = bin;
              ay_index(row, 1) = year;
              row++;
            }
          }

          // Arguments run bin, year, cohort. construct_Q()'s parameter names read the
          // other way round, but ay_index is also filled the other way round from WHAM's,
          // so the two cancel. Do not "correct" the call to match the names.
          Eigen::SparseMatrix<Type> Q_sparse(num_rows, num_rows);
          Q_sparse = construct_Q(nyrs_hind, n_sel_bins, ay_index,
                                 sel_curve_pen(flt, 0), sel_curve_pen(flt, 1), sel_curve_pen(flt, 2), // natural scale
                                 log(square(sel_dev_sd(flt))), 0); // conditional variance

          array<Type> tmp_AR2(n_sel_bins, nyrs_hind);
          for (int bin = 0; bin < n_sel_bins; bin++) {
            for (int yr = 0; yr < nyrs_hind; yr++) {
              tmp_AR2(bin, yr) =  sel_coff_dev(flt, sex, bin, yr);
            }
          }
          jnll_comp(JNLL_SEL_DEV, flt) += GMRF(Q_sparse)(tmp_AR2);
        }
      }
    }
  }


  // 13.8. CATCHABILITY (JNLL_Q_PRIOR, JNLL_Q_DEV)
  // Charged once per Catchability_index, on the lead fleet.
  for(flt = 0; flt < n_flt; flt++){
   if(flt_q_lead(flt) == 1){

    // Lognormal prior (Catchability 2, Estimated-with-prior). Catchability_init is the
    // mean of q when bias_adjust_proc = 1, the median when 0.
    if( est_index_q(flt) == 2){
      jnll_comp(JNLL_Q_PRIOR, flt) -= dnorm(index_log_q(flt), index_log_q_prior(flt) - bias_adjust_proc*square(index_q_sd(flt))/2.0, index_q_sd(flt), true);
    }

    // Deviates, on an estimated q (Catchability 1 or 2) of a survey or CPUE fleet.
    bool q_devs = (flt_type(flt) > 0) && ((est_index_q(flt) == 1) || (est_index_q(flt) == 2));

    // IID (1); the retired AR1 (2) is scored as IID.
    if(q_devs && ((index_varying_q(flt) == 1) || (index_varying_q(flt) == 2))){
      for(yr = 0; yr < nyrs_hind; yr++){
        jnll_comp(JNLL_Q_DEV, flt) -= dnorm(index_q_dev(flt, yr), Type(0.0), index_q_dev_sd(flt), true );
      }
    }

    // RandomWalk (4)
    if(q_devs && (index_varying_q(flt) == 4)){
      for(yr = 1; yr < nyrs_hind; yr++){
        jnll_comp(JNLL_Q_DEV, flt) -= dnorm(index_q_dev(flt, yr) - index_q_dev(flt, yr-1), Type(0.0), index_q_dev_sd(flt), true );
      }
    }
   }
  }


  // 13.9. RECRUITMENT (JNLL_SRR_PRIOR, JNLL_INIT_DEV, JNLL_REC_DEV, JNLL_SRR_PENALTY)
  // Lognormal deviates are centred at -sigma^2/2 when bias_adjust_proc = 1, so the
  // expected recruitment is the curve's (mean-unbiased); at 0 the curve is the median.
  for(sp = 0; sp < nspp; sp++) {
    penalty = 0.0;
    // Input numbers-at-age (estDynamics > 0) have no stock-recruit curve to penalize.
    int srr_terms_on = (estDynamics(sp) == 0);

    // Beverton-Holt steepness prior: lognormal (srr_prior is its mean or median, as above)...
    if(srr_terms_on && (srr_est_mode == 2) & ((srr_pred_fun == 2) | (srr_pred_fun == 3))){
      jnll_comp(JNLL_SRR_PRIOR, sp) -= dnorm(log(steepness(sp, 0)), log(srr_prior(sp)) - bias_adjust_proc*square(srr_prior_sd(sp))/2.0, srr_prior_sd(sp), true);
    }

    // ...or beta, from the prior mean and sd.
    if(srr_terms_on && (srr_est_mode == 3) & ((srr_pred_fun == 2) | (srr_pred_fun == 3))){
      Type beta_alpha = ((1 - srr_prior(sp))/ square(srr_prior_sd(sp)) - 1/srr_prior(sp)) * square(srr_prior(sp));
      Type beta_beta = beta_alpha * (1/srr_prior(sp) - 1);
      jnll_comp(JNLL_SRR_PRIOR, sp) -= dbeta(steepness(sp, 0), beta_alpha, beta_beta, true);
    }

    // Ricker alpha prior, lognormal.
    if(srr_terms_on && (srr_est_mode == 2) & ((srr_pred_fun == 4) | (srr_pred_fun == 5))){
      jnll_comp(JNLL_SRR_PRIOR, sp) -= dnorm(rec_pars(sp, 1), log(srr_prior(sp)) - bias_adjust_proc*square(srr_prior_sd(sp))/2.0, srr_prior_sd(sp), true);
    }

    // Ricker: penalty when 1 / beta (spawning biomass at peak recruitment) exceeds
    // Bmsy_lim / 1e6.
    if(srr_terms_on && (!isNA(Bmsy_lim(sp))) && ((srr_pred_fun == 4) || (srr_pred_fun == 5))){
      Type bmsy = 1.0/exp(rec_pars(sp, 2));
      bmsy =  posfun(Bmsy_lim(sp)/Type(1000000.0) - bmsy, Type(0.001), penalty);
      jnll_comp(JNLL_SRR_PRIOR, sp) += 100 * penalty;
    }

    // Initial-abundance deviates; none under FreeParams (0), Equilibrium (1) or
    // OffsetEquilibrium (5).
    if(initMode > 1 && initMode != 5){
      for(age = 1; age < nages(sp); age++) {
        jnll_comp(JNLL_INIT_DEV, sp) -= dnorm( init_dev(sp, age - 1), -bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp), true);
      }
    }

    // Recruitment deviates.
    for(yr = 0; yr < nyrs_hind; yr++) {
      jnll_comp(JNLL_REC_DEV, sp) -= dnorm( rec_dev(sp, yr),  -bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp), true);
    }

    // Stock-recruit curve penalty (AMAK/Ianelli): recruitment is estimated about a mean
    // (srr_fun 0) and a curve is fitted to it (srr_pred_fun > 0). R_hat is the mean of R
    // when bias_adjust_proc = 1 (Ianelli pm.tpl, Dorn 2002), the median (AMAK) when 0.
    if(srr_terms_on && (srr_fun == 0) & (srr_pred_fun  > 0)){
      for(yr = srr_hat_styr; yr <= srr_hat_endyr; yr++) {
        jnll_comp(JNLL_SRR_PENALTY, sp) -= dnorm( log(R(sp, yr)), log(R_hat(sp,yr)) - bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp), true);
      }
    }
  }


  // 13.10. REFERENCE POINT PENALTIES (JNLL_REFPT_PENALTY)
  // Each finds the F that meets a harvest-control-rule target by penalizing the miss.

  // HCR 1 (CMSY): maximize total terminal-year catch across species, subject to
  // SSB depletion staying above Plimit.
  Type CMSY = 0;
  if(HCR == 1){
    for(fsh_ind = 0; fsh_ind < catch_ctl.rows(); fsh_ind++){
      flt = catch_ctl(fsh_ind, 0) - 1;
      sp = catch_ctl(fsh_ind, 1) - 1;
      flt_yr = catch_ctl(fsh_ind, 2);
      if((flt_type(flt) == 1) && (forecast(sp) == 1) && (estDynamics(sp) == 0)){
        if(flt_yr == projyr){
          CMSY  += catch_hat(fsh_ind);
        }
      }
    }

    jnll_comp(JNLL_REFPT_PENALTY, 0) = - square(CMSY/1000000.0);

    for(sp = 0; sp < nspp; sp++) {
      if((forecast(sp) == 1) && (estDynamics(sp) == 0)){
        penalty = 0.0;
        Type nothing_useful =  posfun( (ssb_depletion(sp, nyrs-1) - Plimit(sp)), Type(0.0001), penalty); (void) nothing_useful;
        jnll_comp(JNLL_REFPT_PENALTY, sp) += 500.0 * square(CMSY/1000.0) * penalty;
      }
    }
  }

  for(sp = 0; sp < nspp; sp++) {
    // Single species: Flimit and Ftarget are SPR fractions, except under HCR 3, where
    // Ftarget is a fraction of SB0.
    if((forecast(sp) == 1) && (msmMode == 0) && (estDynamics(sp) == 0) && ((DynamicHCR == 0) || (DynamicHCR == 1))){

      // ConstantF (2): Flimit is still estimated, for reference.
      if(HCR == 2){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
      }

      // ConstantFSSB (3): SSB at Ftarget of SB0 in the last projection year (static),
      // or every year against dynamic SB0.
      if(HCR == 3){
        if(DynamicHCR == 0){
          jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Ftarget_percent(sp));
        } else {
          for(yr = 1; yr < nyrs; yr++){ // No initial abundance
            jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((DynamicSBF(sp, yr)/DynamicSB0(sp, yr))-Ftarget_percent(sp));
          }
        }
      }

      // ConstantFSPR, NPFMC, PFMC, SESSF (4-7): SPR-based limit and target.
      if(HCR > 3){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((SPRtarget(sp)/SPR0(sp))-Ftarget_percent(sp));
      }
    }

    // Multispecies: per-recruit quantities are undefined under predation, so each
    // target is a fraction of SB0 in the last projection year (2 and 6 Flimit, 3 Ftarget).
    if((forecast(sp) == 1) && (msmMode > 0) && (estDynamics(sp) == 0)){
      if(HCR == 2){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Flimit_percent(sp));
      }
      if(HCR == 3){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Ftarget_percent(sp));
      }
      if(HCR == 6){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Flimit_percent(sp));
      }
    }
  }


  // 13.11. FLOOR PENALTIES (JNLL_ZERO_N_PENALTY): what posfun() adds in the dynamics
  // when it holds numbers-at-age, or recruitment terms, above zero.
  for(sp = 0; sp < nspp; sp++){
    if(estDynamics(sp) == 0){
      jnll_comp(JNLL_ZERO_N_PENALTY, sp) += zero_N_pen(sp);
    }
  }


  // 13.12. NATURAL MORTALITY (JNLL_M_PRIOR, JNLL_M_RE)
  for(sp = 0; sp < nspp; sp++) {

    // Lognormal prior: M_prior is its mean when bias_adjust_proc = 1, its median when 0.
    Type M_prior_mean = log(M_prior(sp)) - bias_adjust_proc*square(M_prior_sd(sp)) / 2.0;

    // On the estimated M1 only: one value (M1_model 1), one per sex (2), or one per sex
    // and age (3).
    if( (M1_use_prior(sp) == 1) && (M2_use_prior(sp) == 0) ) {
      int nsex_tmp = (M1_model(sp) > 1) ? nsex(sp) : 1;
      int nage_tmp = (M1_model(sp) == 3) ? nages(sp) : 1;

      for(int sex = 0; sex < nsex_tmp; sex++) {
        for(int age = 0; age < nage_tmp; age++) {
          jnll_comp(JNLL_M_PRIOR, sp) -= dnorm(log_M1(sp, sex, age), M_prior_mean, M_prior_sd(sp), true);
        }
      }
    }

    // On total M = M1 + M2, every sex, age and hindcast year.
    // FIXME: sex-specific prior?
    if( M2_use_prior(sp) == 1 ) {
      for(int sex = 0; sex < nsex(sp); sex++) {
        for(int age = 0; age < nages(sp); age++) {
          for(int yr = 0; yr < nyrs_hind; yr++) {
            jnll_comp(JNLL_M_PRIOR, sp) -= dnorm(log(M_at_age(sp, sex, age, yr)), M_prior_mean, M_prior_sd(sp), true);
          }
        }
      }
    }

    // M1 random effects (M1_re): 1 IID by age, 2 IID by year, 3 IID by age and year,
    // 4 AR1 by age, 5 AR1 by year, 6 2D AR1; the IID codes are AR1 with rho mapped at
    // its start. sigma (the conditional sd) and rho are per species, shared by sex;
    // the deviates are per sex under M1_model 2. SCALE() takes the marginal sd.
    int num_re_sexes = (M1_model(sp) == 2 && nsex(sp) > 1) ? 2 : 1;

    // By age, constant over years.
    if( (M1_re(sp) == 1) || (M1_re(sp) == 4) ) {
      Type sigma_M = exp(M1_dev_log_sd(sp, 0));
      Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
      Type Sigma_M = pow(pow(sigma_M, 2) / (1.0 - pow(rho_M_a, 2)), 0.5);
      vector<Type> M_re_age(nages(sp));

      for(int sex = 0; sex < num_re_sexes; sex++) {
        for(int age = 0; age < nages(sp); age++) {
          M_re_age(age) = log_M1_dev(sp, sex, age, 0);
        }
        jnll_comp(JNLL_M_RE, sp) += SCALE(AR1(rho_M_a), Sigma_M)(M_re_age);
      }
    }

    // By year, constant over ages.
    if( (M1_re(sp) == 2) || (M1_re(sp) == 5) ) {
      Type sigma_M = exp(M1_dev_log_sd(sp, 0));
      Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
      Type Sigma_M = pow(pow(sigma_M, 2) / (1.0 - pow(rho_M_y, 2)), 0.5);
      vector<Type> M_re_yr(nyrs_hind);

      for(int sex = 0; sex < num_re_sexes; sex++) {
        for(int yr = 0; yr < nyrs_hind; yr++) {
          M_re_yr(yr) = log_M1_dev(sp, sex, 0, yr);
        }
        jnll_comp(JNLL_M_RE, sp) += SCALE(AR1(rho_M_y), Sigma_M)(M_re_yr);
      }
    }

    // By age and year.
    if( (M1_re(sp) == 3) || (M1_re(sp) == 6) ) {
      Type sigma_M = exp(M1_dev_log_sd(sp, 0));
      Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
      Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
      Type Sigma_M = pow(pow(sigma_M, 2) / ((1.0 - pow(rho_M_y, 2)) * (1.0 - pow(rho_M_a, 2))), 0.5);
      array<Type> M_re_a_yr(nages(sp), nyrs_hind);

      for(int sex = 0; sex < num_re_sexes; sex++) {
        for(int age = 0; age < nages(sp); age++) {
          for(int yr = 0; yr < nyrs_hind; yr++) {
            M_re_a_yr(age, yr) = log_M1_dev(sp, sex, age, yr);
          }
        }
        // SEPARABLE(f, g): f correlates the last array index, g the first; the array is
        // (age, year), so the year correlation goes first.
        jnll_comp(JNLL_M_RE, sp) += SCALE(SEPARABLE(AR1(rho_M_y), AR1(rho_M_a)), Sigma_M)(M_re_a_yr);
      }
    }
  }


  // 13.13. LINKAGE PRIORS (JNLL_LINKAGE_PRIOR)
  // Families are in linkage_log_prior() (likelihood.hpp). A prior on a log-stored base (R0,
  // alpha, M1, K, q, slope, apical) is on its natural scale; one on a sel inflection or
  // DoubleNormalSS3 parameter is on its stored scale. Most (Intercept) rows are held at 0 (2.3b), so
  // the prior reads the base parameter they stand in for; a slope row reads its coefficient.
  for (int i = 0; i < beta_linkage.size(); ++i) {
    int fam = linkage_prior_family(i);
    if (fam == 0) continue;   // 0 = no prior
    Type p1  = linkage_prior_p1(i);
    Type p2  = linkage_prior_p2(i);
    int slot_col = linkage_species(i) - 1;   // jnll column: the row's species, or 0 if it names none
    if (slot_col < 0) slot_col = 0;        // all species => column 0
    if (slot_col >= n_col) slot_col = 0;

    Type b = beta_linkage(i);
    bool base_is_log = true;   // most base parameters are stored as logs
    if (linkage_is_intercept(i) == 1) {
      // Species/sex/age 0 means "all"; the prior reads the first level's base parameter.
      int sp_in = linkage_species(i);
      int sx_in = linkage_sex(i);
      int ab_in = linkage_age_bin(i);
      int sp_idx = (sp_in == 0) ? 0 : (sp_in - 1);
      int sx_idx = (sx_in == 0) ? 0 : (sx_in - 1);
      int ab_idx = (ab_in == 0) ? 0 : (ab_in - 1);

      int proc  = linkage_process(i);
      int param = linkage_param(i);
      if (proc == RCEATTLE_PROC_RECRUIT) {
        // rec_pars columns: log R0, log alpha, log Beta (stock-recruit). R_init has no base parameter;
        // its intercept is the level, so its prior stays on beta_linkage(i).
        if (param != RCEATTLE_REC_R_INIT) {
          b = rec_pars(sp_idx, param);
        }
      } else if (proc == RCEATTLE_PROC_M) {
        b = log_M1(sp_idx, sx_idx, ab_idx);
      } else if (proc == RCEATTLE_PROC_GROWTH) {
        // K, L1, Linf, m in log_growth_pars; then sd_L1, sd_Linf in growth_log_sd.
        if (param < RCEATTLE_N_GROWTH_PARAMS) {
          b = log_growth_pars(sp_idx, sx_idx, param);
        } else {
          b = growth_log_sd(sp_idx, sx_idx, param - RCEATTLE_N_GROWTH_PARAMS);
        }
      } else if (proc == RCEATTLE_PROC_Q) {
        // Fleet-indexed: log q.
        int fl_in  = linkage_fleet(i);
        int fl_idx = (fl_in == 0) ? 0 : (fl_in - 1);
        b = index_log_q(fl_idx);
      } else if (proc == RCEATTLE_PROC_COMP) {
        // Dirichlet-multinomial log theta: comp and CAAL by fleet, diet by predator.
        int fl_in  = linkage_fleet(i);
        int fl_idx = (fl_in == 0) ? 0 : (fl_in - 1);
        if (param == RCEATTLE_COMP_THETA_COMP)      b = comp_weights(fl_idx);
        else if (param == RCEATTLE_COMP_THETA_CAAL) b = caal_weights(fl_idx);
        else if (param == RCEATTLE_COMP_THETA_DIET) b = diet_comp_weights(sp_idx);
      } else if (proc == RCEATTLE_PROC_SEL) {
        // Parametric selectivity, element 0 ascending, 1 descending. The (Intercept)
        // adds no offset, so the prior regularizes the base parameter (as in Cole
        // Monnahan's GOA pollock priors). RCEATTLE_SEL_COFF has no scalar base and is
        // refused by .check_sel_linkage_support().
        int fl_in  = linkage_fleet(i);
        int fl_idx = (fl_in == 0) ? 0 : (fl_in - 1);
        if (param == RCEATTLE_SEL_SLP_ASC || param == RCEATTLE_SEL_SLP_DESC) {
          b = log_sel_slp(param - RCEATTLE_SEL_SLP_ASC, fl_idx, sx_idx);   // log slope
        } else if (param == RCEATTLE_SEL_INF_ASC || param == RCEATTLE_SEL_INF_DESC) {
          // As stored: an inflection, DoubleNormal's peak or logit floor, or LogisticPM's
          // log first-bin selectivity.
          b = sel_inf(param - RCEATTLE_SEL_INF_ASC, fl_idx, sx_idx);
          base_is_log = false;
        } else if (param == RCEATTLE_SEL_APICAL) {
          b = log_sel_apical(fl_idx, sx_idx);             // log apical multiplier; 0 = a multiplier of 1
        } else if (param >= RCEATTLE_SEL_DN_PEAK && param <= RCEATTLE_SEL_END_LOGIT) {
          // DoubleNormalSS3, on SS3's own scale (peak in bins or cm, log widths, logits),
          // where an SS3 control file puts its priors.
          b = sel_dn6(param - RCEATTLE_SEL_DN_PEAK, fl_idx, sx_idx);
          base_is_log = false;
        }
      }
    }

    Type log_prior = linkage_log_prior(fam, b, (linkage_is_intercept(i) == 1) && base_is_log,
                                       p1, p2, bias_adjust_proc);
    jnll_comp(JNLL_LINKAGE_PRIOR, slot_col)            -= log_prior;
    unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, slot_col) -= log_prior;
  }

  // Priors on each random-effect group's sd, from linkage_spec(priors = list(sigma = )).
  // FIXME(jacobian): normal, gamma and beta read the sd while log_sigma is estimated,
  // with no Jacobian -- a penalty, not a proper prior. Add log_sigma for a density.
  if (log_sigma_linkage.size() > 0) {
    for (int g = 0; g < log_sigma_linkage.size(); ++g) {
      int fam = linkage_re_sigma_prior_family(g);
      if (fam == 0) continue;   // 0 = no prior
      Type sd = exp(log_sigma_linkage(g));
      Type log_prior = linkage_log_prior(fam, sd, false, linkage_re_sigma_prior_p1(g),
                                         linkage_re_sigma_prior_p2(g), bias_adjust_proc);
      jnll_comp(JNLL_LINKAGE_PRIOR, 0)            -= log_prior;
      unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, 0) -= log_prior;
    }
  }

  // Priors on each AR1 group's correlation rho in (-1, 1), from priors = list(rho = ):
  // normal on rho, or beta on (rho + 1) / 2.
  // FIXME(jacobian): as for sigma, no rho_trans Jacobian.
  if (trans_rho_linkage.size() > 0) {
    for (int g = 0; g < linkage_re_rho.size(); ++g) {
      int fam = linkage_re_rho_prior_family(g);
      if ((fam != 1) && (fam != 4)) continue;   // only normal (1) and beta (4) are read for rho
      Type rho = rho_trans(trans_rho_linkage(linkage_re_rho(g)));
      Type x = (fam == 4) ? (rho + Type(1)) / Type(2) : rho;
      Type log_prior = linkage_log_prior(fam, x, false, linkage_re_rho_prior_p1(g),
                                         linkage_re_rho_prior_p2(g), bias_adjust_proc);
      jnll_comp(JNLL_LINKAGE_PRIOR, 0)            -= log_prior;
      unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, 0) -= log_prior;
    }
  }


  // 13.14. LINKAGE RANDOM EFFECTS (JNLL_LINKAGE_RE)
  // Group g's deviations u, in time order: IID u_t ~ N(0, sigma_g); random walk
  // u_t - u_{t-1} ~ N(0, sigma_g); or stationary AR1, marginal sd sigma_g
  // (rceattle_re_group_nll()). Integrated and penalized deviations share this density.
  for (int grp = 0; grp < log_sigma_linkage.size(); ++grp) {
    vector<int> slots = rceattle_group_slots(linkage_re_sigma, grp);
    int len = slots.size();
    if (len == 0) continue;
    vector<Type> re(len);
    for (int t = 0; t < len; ++t) re(t) = beta_linkage_re_all(slots(t));

    int st = linkage_re_struct(grp);
    Type sigma = exp(log_sigma_linkage(grp));
    Type rho = (st == RCEATTLE_RE_AR1) ? rho_trans(trans_rho_linkage(linkage_re_rho(grp))) : Type(0);
    rceattle_re_group_nll(st, re, sigma, rho,
                          jnll_comp(JNLL_LINKAGE_RE, 0), unweighted_jnll_comp(JNLL_LINKAGE_RE, 0));

    // Rogers QAR1: covariate x(t) ~ N(re(t), obs_sd), obs_sd fixed unless obs_sd_est =
    // TRUE. Tying re to x identifies beta_linkage_obs, which multiplies re (5.2b).
    if (linkage_re_obs(grp) >= 0) {
      Type osd = exp(log_obs_sd_linkage(linkage_re_obs(grp)));
      for (int t = 0; t < len; ++t) {
        int s = slots(t);
        if (linkage_re_obs_mask(s) == 0) continue;   // no observation that year
        if (osa_mode == 0) {
          Type log_dens = dnorm(linkage_re_obs_value(s), re(t), osd, true);
          jnll_comp(JNLL_LINKAGE_RE, 0)            -= log_dens;
          unweighted_jnll_comp(JNLL_LINKAGE_RE, 0) -= log_dens;
        } else {
          // OSA (osa_mode 1-2): the same x(t), from obsvec so oneStepPredict() can hold
          // it out through keep(pos); residualized against re(t), as WHAM's Ecov.
          int pos = linkage_re_obsvec_idx(s);
          if (pos >= 0) {
            jnll_comp(JNLL_LINKAGE_RE, 0)            -= keep(pos) * dnorm(obsvec(pos), re(t), osd, true);
            unweighted_jnll_comp(JNLL_LINKAGE_RE, 0) -= keep(pos) * dnorm(obsvec(pos), re(t), osd, true);
            // The CDF term is not a likelihood, so only jnll_comp gets it here.
            if (osa_mode == 2) {
              jnll_comp(JNLL_LINKAGE_RE, 0) -= osa_norm_cdf_terms(
                obsvec(pos), re(t), osd, keep.cdf_lower(pos), keep.cdf_upper(pos));
            }
          }
        }
      }
    }
  }


  // 13.15. UNWEIGHTED COMPONENTS
  // jnll_comp without data weights. In ordinary fitting the composition, CAAL and
  // stomach rows were written unweighted above (0 under OSA), and the linkage rows to
  // both; every other row carries no data weight, so it is copied.
  for(flt = 0; flt < n_col; flt++){
    for(int ind = 0; ind < JNLL_N_ROWS; ind++){
      bool written_above = (ind == JNLL_COMP) || (ind == JNLL_CAAL) || (ind == JNLL_STOMACH) ||
                           (ind == JNLL_LINKAGE_PRIOR) || (ind == JNLL_LINKAGE_RE);
      if(!written_above){
        unweighted_jnll_comp(ind, flt) = jnll_comp(ind, flt);
      }
    }
  }


  /** ------------------------------------------------------------------------ //
   // 14. REPORT SECTION                                                        //
   * ------------------------------------------------------------------------- */

  // 14.1. Population components
  REPORT( pop_scalar );
  ADREPORT( pop_scalar );
  REPORT( avg_R );
  REPORT( Z_at_age );
  REPORT( N_at_age );
  REPORT( avgN_at_age );
  REPORT( biomass_at_age );
  REPORT( biomass_depletion );
  REPORT( ssb_depletion );
  REPORT( biomass );
  REPORT( ssb );
  REPORT( exploitable_biomass );
  REPORT(R_sd);
  REPORT( R0 );
  REPORT( R_init );
  REPORT( R );
  REPORT( R_hat );
  REPORT( M_at_age );
  //REPORT( mature_females );
  //REPORT( sex_ratio_hat );
  //REPORT( N_at_age_dB0);
  //REPORT( N_at_age_dBF);
  // ADREPORT( M_at_age );            // no SEs: see ?plot_mortality
  REPORT( weight_hat );
  REPORT( growth_matrix );
  REPORT( length_hat );

  ADREPORT( biomass );
  ADREPORT( ssb );

  // -- Secondary display series, so add_ci = TRUE can draw an interval at all.
  // Natural scale ONLY: exploitable_biomass sums over fisheries alone, so it is
  // exactly 0 for a survey-only species or proj_F_prop = 0, and the depletions
  // divide by B0 / SB0. log() of either puts -Inf on the tape and NaNs the whole
  // sdreport; the plotters recover sd(log x) = sd(x)/x instead.
  ADREPORT( exploitable_biomass );
  ADREPORT( biomass_depletion );
  ADREPORT( ssb_depletion );

  // -- Log-scale derived series.
  // sdreport linearizes once about the MLE, so its SD is exact only for a linear function
  // of the parameters. These are built multiplicatively (R = R0 * exp(rec_dev); n-at-age
  // is a product of survivals), so log(x) is near-linear and exp(log(x) +/- z * sd) is
  // right-skewed and cannot cross zero. The log SD is also the CV, the form the ABC / OFL
  // buffer wants. The natural-scale ADREPORTs stay so existing sdrep$value callers work.
  matrix<Type>  log_biomass = biomass;  log_biomass = log(biomass.array());
  matrix<Type>  log_ssb     = ssb;      log_ssb     = log(ssb.array());
  matrix<Type>  log_R       = R;        log_R       = log(R.array());
  REPORT( log_biomass );
  REPORT( log_ssb );
  REPORT( log_R );
  ADREPORT( log_biomass );
  ADREPORT( log_ssb );
  ADREPORT( log_R );
  ADREPORT( R_sd );
  ADREPORT( R );

  // -- Log selectivity-at-age, for a delta-method interval on the curve. Log and not
  // logit, because the non-parametric forms normalize to a MEAN of one, so sel_at_age
  // routinely exceeds one. Gated because the delta method costs every reported value
  // against every parameter; see ?fit_control.
  //
  // Only estimated, age-based lead fleets, from their first selected bin: those cells
  // cannot be a structural zero, and one log(0) = -Inf on the tape turns EVERY sdreport
  // quantity to NaN, biomass and SSB included. The four routes to a structural zero are
  // enumerated in ?fit_control and all read from the data, so the reported SET never
  // depends on a parameter value. A reported VALUE can still underflow -- a narrow dome
  // puts exp(-x^2/2) below 1e-308 many bins from its peak -- so the log is floored with
  // max2(), exact above the floor. log_sel_at_age_index carries each row's fleet, sex,
  // absolute age and calendar year, since which cells are reported depends on the model.
  if(adreport_sel == 1){
    // First bin and bin count per fleet. data_check() bounds Bin_first_selected
    // by nages, but a 0 in the workbook still arrives here as -1, and the model
    // builds safebounds = FALSE, so both ends are clamped before indexing.
    vector<int> sel_first_bin(n_flt); sel_first_bin.setZero();
    vector<int> sel_rep_bins(n_flt);  sel_rep_bins.setZero();
    for(int flt = 0; flt < n_flt; flt++){
      if((flt_type(flt) <= 0) || (flt_sel_lead(flt) != 1) ||
         (flt_sel_type(flt) == 0) || (flt_sel_dim(flt) == 1)) continue;
      int n_age = nages(flt_spp(flt));
      int first = bin_first_selected(flt);
      if(first < 0) first = 0;
      if(first > n_age) first = n_age;
      sel_first_bin(flt) = first;
      sel_rep_bins(flt)  = n_age - first;
    }

    int n_sel = 0;
    for(int flt = 0; flt < n_flt; flt++){
      n_sel += nsex(flt_spp(flt)) * sel_rep_bins(flt) * nyrs_hind;
    }

    vector<Type> log_sel_at_age(n_sel); log_sel_at_age.setZero();
    matrix<Type> log_sel_at_age_index(n_sel, 4); log_sel_at_age_index.setZero();

    int i_sel = 0;
    for(int flt = 0; flt < n_flt; flt++){
      if(sel_rep_bins(flt) == 0) continue;
      int sel_sp = flt_spp(flt);
      for(int sel_sex = 0; sel_sex < nsex(sel_sp); sel_sex++){
        for(int sel_age = sel_first_bin(flt); sel_age < nages(sel_sp); sel_age++){
          for(int sel_yr = 0; sel_yr < nyrs_hind; sel_yr++){
            log_sel_at_age(i_sel) = log(max2(sel_at_age(flt, sel_sex, sel_age, sel_yr),
                                             Type(1e-300)));
            log_sel_at_age_index(i_sel, 0) = Type(flt + 1);
            log_sel_at_age_index(i_sel, 1) = Type(sel_sex + 1);
            log_sel_at_age_index(i_sel, 2) = Type(minage(sel_sp) + sel_age);
            log_sel_at_age_index(i_sel, 3) = Type(styr + sel_yr);
            i_sel++;
          }
        }
      }
    }

    REPORT( log_sel_at_age );
    REPORT( log_sel_at_age_index );
    ADREPORT( log_sel_at_age );
  }


  // -- 14.2. Biological reference points
  REPORT( NByageF );
  REPORT( NByage0);
  REPORT( NbyageSPR);
  REPORT( B0 );
  REPORT( SB0 );
  REPORT( SBF );
  REPORT( DynamicB0 );
  REPORT( DynamicSB0 );
  REPORT( DynamicSBF );
  REPORT( SPR0 );
  REPORT( SPRFinit );
  REPORT( SPRlimit );
  REPORT( SPRtarget );
  REPORT( steepness );
  REPORT( proj_F );
  REPORT( Ftarget );
  REPORT( Flimit );

  // -- 14.3. Selectivity
  REPORT( sel_at_age );
  REPORT( sel_at_length );


  // -- 14.4. Survey components
  REPORT( index_hat );
  REPORT( index_sd );
  REPORT( index_q );
  vector<Type>  log_index_hat = index_hat;  log_index_hat = log(index_hat.array());// Log predicted index
  REPORT( log_index_hat );
  ADREPORT( log_index_hat );

  // -- 14.4b. Formula-linkage effect sizes: the estimated coefficients, the random-effect
  // deviations and the Rogers (2024) QAR1 effect size, so a user can read the *effect*
  // and, through ADREPORT, its uncertainty -- not just the realized parameter (index_q
  // above). All may be length 0 (no linkages), in which case they report nothing.
  REPORT( beta_linkage );
  ADREPORT( beta_linkage );
  REPORT( beta_linkage_re );
  REPORT( beta_linkage_re_pen );
  // Slot-ordered deviations across both vectors, so a caller can read a group's
  // walk without knowing how it is stored.
  REPORT( beta_linkage_re_all );
  REPORT( beta_linkage_obs );
  ADREPORT( beta_linkage_obs );


  // -- 14.5. Fishery components
  REPORT( F_spp );
  REPORT( F_flt );
  REPORT( F_flt_age );
  REPORT( F_at_age );
  REPORT( catch_hat );
  REPORT( max_catch_hat );
  REPORT( catch_sd );
  REPORT( catch_analytical_sd );


  // 14.6. Age/length composition
  REPORT( age_obs_hat );
  REPORT( age_hat );
  REPORT( comp_obs );
  REPORT( comp_hat );
  REPORT( caal_hat );
  REPORT( equil_catch_hat );
  REPORT( caal_obs );
  REPORT( pred_CAAL );


  // -- 14.7. Likelihood components
  REPORT( jnll_comp );
  REPORT( unweighted_jnll_comp );


  // -- 14.8. Ration components
  REPORT( consumption_at_age );
  REPORT( fT );

  // -- 14.9. Predation components
  REPORT( M2_prop );
  REPORT( diet_obs );
  REPORT( suitability );
  REPORT( M1_at_age );
  REPORT( M2_at_age );
  REPORT( B_eaten );
  REPORT( B_eaten_as_prey );
  // ADREPORT( B_eaten_as_prey );    // no SEs: see ?plot_b_eaten_prop
  REPORT( vulnerability );
  REPORT( vulnerability_other );
  REPORT( gam_a );
  REPORT( gam_b );
  REPORT( diet_hat );
  REPORT( suit_other );
  REPORT( stom_div_bio );
  REPORT( avail_food );


  // -- 14.10. Initial-state cumulative mortality (feeds the initial age structure)
  REPORT(mort_sum);

  /** ------------------------------------------------------------------------ //
   // 15. END MODEL                                                             //
   * ------------------------------------------------------------------------- */

  Type jnll = 0;

  // Modes 0-3 return the real objective. Mode 3 builds the TMB object without
  // running the optimizer, so obj$fn()/obj$gr() are usable for diagnostics.
  if(estimateMode < 4) {
    jnll = jnll_comp.sum();
  }

  // Mode 4: build_map() maps out every hindcast parameter, leaving `dummy` as
  // the only free parameter. Placeholder objective so nlminb has something
  // valid to minimize -- a plumbing smoke test, not a likelihood.
  if(estimateMode > 3) {
    jnll = dummy * dummy;
  }


  REPORT( jnll );
  return jnll;
}

