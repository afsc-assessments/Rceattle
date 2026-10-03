#ifndef GROWTH_HPP
#define GROWTH_HPP

/**
 * @brief SD of length-at-age (cm) for one age, year and sex.
 *
 * The variability term v is \f$e^{sd_0}\f$ at or below `age_L1`, pinned to
 * \f$e^{sd_1}\f$ for the plus group under `sd_style == 1` (WHAM; also SS3 when
 * Growth_Age_for_L2 = 999), and otherwise interpolated linearly in length
 * between \f$(l_1, e^{sd_0})\f$ and \f$(L_\infty, e^{sd_1})\f$. `sd_form` says what v is:
 * 1 = an SD in cm (SS3 CV_Growth_Pattern 2), 2 = a CV, so SD = v * L
 * (SS3 CV_Growth_Pattern 0).
 */
template<class Type>
Type length_sd_at_age(Type current_age, Type age_L1, bool plus_group,
                      int sd_style, int sd_form, Type l1, Type linf, Type len,
                      Type log_sd0, Type log_sd1) {
  Type sd0 = exp(log_sd0);
  Type sd1 = exp(log_sd1);
  Type v;
  if(current_age <= age_L1) {
    v = sd0;
  } else if(plus_group && sd_style == 1) {
    v = sd1;
  } else {
    v = sd0 + (sd1 - sd0) / (linf - l1) * (len - l1);
  }
  if(sd_form == 2) v = v * len;
  return v;
}


/**
 * @brief Age-length key, weight-at-age and mature weight-at-age for one cell.
 *
 * Length is normal with mean `mu` and SD `sd` (cm). Probabilities are taken on
 * the population length bins `lengths_pop` (lower edges): the first bin is a
 * minus group below the second edge and the last a plus group above the last
 * edge (SS3's convention). They stay on that grid; the consumers that need data
 * bins sum each bin's run of population bins themselves.
 *
 * Weight-at-age is \f$\sum_l P(l|a) \alpha L_{mid}^\beta\f$ at population-bin
 * midpoints (kg); the last bin's midpoint sits half a bin width above its lower
 * edge. When the species has maturity-at-length, mature weight-at-age is the
 * same sum with the logistic maturity \f$1/(1+e^{-s(L_{mid} - L_{50})})\f$ inside:
 * SS3's fecundity-at-age when fecundity equals body weight.
 */
template<class Type>
void fill_age_length_key(int wtind, int sp, int sex, int age, int yr,
                         Type mu, Type sd,
                         const vector<int>& nlengths_pop,
                         matrix<Type>& lengths_pop,
                         matrix<Type>& weight_length_pars,
                         const vector<int>& mat_len_use,
                         matrix<Type>& mat_len_pars,
                         array<Type>& growth_matrix,
                         array<Type>& weight_hat,
                         array<Type>& mat_weight_hat) {
  int np = nlengths_pop(sp);
  // Kept on the POPULATION bins. Its readers multiply it by
  // selectivity-at-length, and SS3 forms that product at population resolution
  // before binning, so aggregating here would apply one bin-average selectivity
  // across a data bin. ceattle.cpp sums each data bin's run of population bins
  // (section 2.3c) where it needs data bins.
  for(int lp = 0; lp < np; lp++) growth_matrix(wtind, sex, age, lp, yr) = Type(0.0);

  Type expected_weight = 0.0;
  Type expected_mat_weight = 0.0;
  for(int lp = 0; lp < np; lp++) {
    Type prob;
    if(np < 2) {
      // One bin spans the whole length distribution, so it holds all of it.
      // Tested before the lp == 0 branch, which would read the second edge.
      prob = Type(1.0);
    } else if(lp == 0) {
      prob = pnorm((lengths_pop(sp, 1) - mu) / sd);
    } else if(lp == (np - 1)) {
      prob = 1.0 - pnorm((lengths_pop(sp, np - 1) - mu) / sd);
    } else {
      prob = pnorm((lengths_pop(sp, lp + 1) - mu) / sd) - pnorm((lengths_pop(sp, lp) - mu) / sd);
    }
    growth_matrix(wtind, sex, age, lp, yr) = prob;

    // Midpoint (cm) of the bin: lower edge plus half its own width, the last bin
    // taking the previous bin's width. A single-bin grid has no neighbour to take
    // a width from, so it keeps the edge itself -- as selectivity.hpp does.
    Type lenmid;
    if(lp < np - 1) {
      lenmid = (lengths_pop(sp, lp) + lengths_pop(sp, lp + 1)) / Type(2.0);
    } else if(np < 2) {
      lenmid = lengths_pop(sp, np - 1);
    } else {
      lenmid = lengths_pop(sp, np - 1) + (lengths_pop(sp, np - 1) - lengths_pop(sp, np - 2)) / Type(2.0);
    }
    Type wt_len = weight_length_pars(sp, 0) * pow(lenmid, weight_length_pars(sp, 1));
    expected_weight += prob * wt_len;
    if(mat_len_use(sp) == 1) {
      expected_mat_weight += prob * wt_len / (Type(1.0) + exp(-mat_len_pars(sp, 1) * (lenmid - mat_len_pars(sp, 0))));
    }
  }
  weight_hat(wtind, sex, age, yr) = expected_weight;
  if(mat_len_use(sp) == 1) mat_weight_hat(wtind, sex, age, yr) = expected_mat_weight;
}

/**
 * @brief Integrated Growth, Size-Transition, and Weight-at-Age Module for Month = 0.
 *
 * Computes Jan-1 mean length-at-age, the age->length probability matrix, and
 * integrated weight-at-age for a single species.
 *
 * @section math_models Mathematical Models:
 * 1. Mean Length-at-Age (\f$L_a\f$) — anchored at \f$L(a_{L1}) = l_1\f$:
 *    - Von Bertalanffy (Model 1): \f$L_a = L_{\infty} + (l_1 - L_{\infty}) \cdot e^{-K(a - a_{L1})}\f$
 *    - Richards (Model 2): \f$L_a = [L_{\infty}^m + (l_1^m - L_{\infty}^m) \cdot e^{-K(a - a_{L1})}]^{1/m}\f$
 *    - For year > 0, a cohort recursion advances \f$L_{a, y}\f$ from \f$L_{a-1, y-1}\f$
 *      using lag-year parameters. At the cohort boundary (current_age ==
 *      age_L1_ceil) the closed-form anchor at \f$l_1\f$ is used; `age_L1_safe`
 *      keeps this anchor at \f$l_1\f$ for both minage > 0 and minage = 0.
 *    - Linear ramp from Lmin_sp at age 0 to \f$l_1\f$ at age_L1 applies when
 *      current_age <= age_L1 (only reachable for minage > 0).
 *
 * 2. Weight-at-Age (\f$W_a\f$):
 *    Integrated across length to account for Jensen's Inequality:
 *    \f$W_a = \sum_{ln} P(ln | a) \cdot \alpha \cdot L_{mid}^{\beta}\f$
 *    Bin midpoints are computed per-bin to support non-uniform length bins;
 *    the length plus-group is extended by half the final interior bin width.
 *
 * @section logic Biological Logic:
 * - **Temporal Resolution**: Jan-1 (month = 0).
 * - **Plus-Age Group**: The oldest age's mean length is set by
 *   `growth_plus_length`: an M1-weighted mean toward \f$L_{\infty}\f$ (1, the
 *   default), no adjustment (2), or SS3's Linf_decay forms (3 = -999, 4 = a
 *   decay rate). See section 2 of the body.
 * - **SD-at-Age**: length_sd_at_age(); `growth_sd_form` makes the two
 *   endpoints SDs (cm) or CVs.
 * - **Size Transition**: fill_age_length_key(), on the population length
 *   bins; consumers bin the product with selectivity, not the key itself.
 *
 * @param wtind Weight index slot to write into.
 * @param sp Species index.
 * @param nyrs Number of years in the simulation/model.
 * @param nsex/nages/nlengths Vectors of species-specific counts per dimension.
 * @param minage Vector of minimum modeled ages (0 supported).
 * @param growth_model Per-species growth model selector (1 = VBGF, 2 = Richards).
 * @param lengths Matrix of length bin boundaries (nspp x length_bins). Bin
 *                widths may be non-uniform.
 * @param growth_parameters 4D array of time-varying growth parameters (K, L1, Linf, m).
 * @param growth_log_sd 3D array of log-scale SDs (index 0 = SD at L1, index 1 = SD at Linf).
 * @param weight_length_pars Matrix of length-weight parameters (\f$\alpha\f$, \f$\beta\f$).
 *
 * @note Outputs are passed by reference.
 * @param length_hat [Output] 4D array filled with mean length-at-age.
 * @param growth_matrix [Output] 5D array filled with the growth transition matrix.
 * @param weight_hat [Output] 4D array filled with integrated weight-at-age.
 */
template<class Type>
void estimate_growth(
    int wtind,
    int sp,
    int nspp,
    int nyrs,
    const vector<int>&  nsex,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  minage,
    const vector<Type>& growth_age_L1,
    const vector<int>&  growth_model,
    const vector<int>&  growth_sd_style,   // Plus-group SD-at-age: 1 = WHAM (pin to exp(sd_Linf)), 2 = interpolate by length
    const vector<int>&  growth_sd_form,    // 1 = SD endpoints in cm, 2 = CV endpoints (SD = CV * L)
    const vector<int>&  growth_plus_length,// Plus-group mean length: 1 = M1-weighted, 2 = none, 3 = SS3.24, 4 = decay
    const vector<Type>& plus_group_decay,  // Decay rate (per year) for growth_plus_length == 4
    matrix<Type>& lengths,
    const vector<int>&  nlengths_pop,
    matrix<Type>& lengths_pop,             // Population length bins, lower edges (cm)
    array<Type>& growth_parameters,
    array<Type>& growth_log_sd,
    matrix<Type>& weight_length_pars,
    const vector<int>&  mat_len_use,       // 1 = maturity-at-length for this species
    matrix<Type>& mat_len_pars,            // [sp, 0] = L50 (cm), [sp, 1] = logistic slope (per cm)
    array<Type>& log_M1,         // Base natural mortality at age [nspp, nsex, nages], log scale
    array<Type> &length_hat,     // Modified by reference
    array<Type> &growth_matrix,  // Modified by reference
    array<Type> &weight_hat,     // Modified by reference
    array<Type> &mat_weight_hat  // Modified by reference (maturity-at-length species only)
) {

  // Initialize output and temporary storage


  // Calculate mean-length, SD, and growth matrix, for all years:
  // lengths is vector with lengths mm (2, 4, 6, 8, etc)
  Type b_len, last_linear, current_age;

  // The linear ramp below age_L1 starts from the lowest population length
  // edge at age 0 (SS3's len_bins(1)); with no population grid that is the
  // lowest data bin.
  Type Lmin_sp = lengths_pop(sp, 0);
  // age_L1 is the VB anchor age (= age at which `l1` is the length). Read
  // from data_list$growth_age_L1[sp] (= SS3's Growth_Age_for_L1 ctl input).
  // R-side fit_mod() resolves the default to max(0.5, minage[sp]) so models
  // with minage = 0 get an SS3-style half-year anchor and minage >= 1 stays
  // backwards-compatible.
  Type age_L1 = growth_age_L1(sp);
  // age_L1_ceil is compared to current_age (the slot age in years) to detect
  // the youngest VB-relevant slot. That slot needs the closed-form boundary
  // formula (not the cohort recursion which would index age-1 = -1 at slot 0
  // when minage = 0, or use a stale value at slot 1 when the anchor is at
  // age 0.5). For minage = 0, the youngest VB slot is C++ age = 1
  // (current_age = 1, anchor at 0.5). For minage >= 1, slot 0 itself
  // (current_age = minage). Either way: age_L1_ceil = max(1, minage).
  Type age_L1_ceil = (Type(minage(sp)) >= Type(1)) ? Type(minage(sp)) : Type(1);

  // When minage = 0, the linear ramp from Lmin_sp at age 0 to l1 at age_L1
  // becomes degenerate (divide by zero in b_len) AND the ramp branch should
  // never execute (current_age >= 1 always > age_L1 = 0). Replace the ramp
  // anchor with a small positive value so b_len is finite; the value is
  // unused because the if(current_age <= age_L1) branch is unreachable.
  Type age_L1_safe = age_L1;
  if (age_L1 <= Type(0)) age_L1_safe = Type(1);  // safe denominator only

  for(int sex = 0; sex < nsex(sp); sex++) {
    for(int yr = 0; yr < nyrs; yr++) {
      for(int age = 0; age < nages(sp); age++) {

        // Parameters for parametric growth
        Type kappa = growth_parameters(sp, sex, yr, 0);
        Type l1 = growth_parameters(sp, sex, yr, 1);
        Type linf = growth_parameters(sp, sex, yr, 2);
        Type m = growth_parameters(sp, sex, yr, 3);
        // current_age is the slot's actual age (start of year). Slot index
        // age (0-based) holds the (age + minage) cohort. The historical
        // `age + 1.0` shifted everything by one year at any minage != 1
        // (gave end-of-year length instead of start). Note this is the SLOT
        // age, NOT the VB anchor (which is age_L1, defined above).
        current_age = Type(age) + Type(minage(sp));

        // 1. Calculate Mean Length at Age ---
        switch(growth_model(sp)) {
        case 1: // Von Bertalanffy

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 when minage > 0; == 1 when minage = 0 (ramp branch unreachable)

          // Age < minage
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else {
            if(yr == 0) {
              length_hat(wtind,  sex, age, yr) = linf + (l1 - linf) * exp(-kappa * (current_age - age_L1));
            } else { // Yr > 0
              if(current_age == age_L1_ceil) {
                // age_L1_safe ensures last_linear == l1 in both minage > 0
                // (where age_L1_safe == age_L1) and minage = 0 (where
                // age_L1 = 0 would otherwise collapse this to Lmin_sp).
                last_linear = Lmin_sp + b_len * age_L1_safe;
                length_hat(wtind,  sex, age, yr) = last_linear + (last_linear - linf) * (exp(-kappa * (current_age - age_L1)) - 1.0);
              } else {
                // Lag 1-year parameters
                Type lagkappa = growth_parameters(sp, sex, yr - 1, 0);
                Type laglinf = growth_parameters(sp, sex, yr - 1, 2);
                length_hat(wtind,  sex, age, yr) = length_hat(wtind,  sex, age - 1, yr - 1) + (length_hat(wtind,  sex, age - 1, yr - 1) - laglinf) * (exp(-lagkappa) - 1.0);
              }
            }
          }
          break;

        case 2: // Richards

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 when minage > 0; == 1 when minage = 0 (ramp branch unreachable)

          if((current_age) <= age_L1) {
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else {
            if(yr == 0) {
              length_hat(wtind,  sex, age, yr) = pow(pow(linf, m) + (pow(l1, m) - pow(linf, m)) * exp(-kappa * (current_age - age_L1)), 1 / m);
            } else { // Yr > 0
              if(current_age == age_L1_ceil) {
                // age_L1_safe ensures last_linear == l1 in both minage > 0
                // and minage = 0; see VBGF branch above for full reasoning.
                last_linear = Lmin_sp + b_len * age_L1_safe;
                length_hat(wtind,  sex, age, yr) = pow(pow(last_linear, m) + (pow(last_linear, m) - pow(linf, m)) * (exp(-kappa * (current_age - age_L1)) - 1.0), 1 / m);
              } else {
                // Lag 1-year parameters
                Type lagkappa = growth_parameters(sp, sex, yr - 1, 0);
                Type laglinf = growth_parameters(sp, sex, yr - 1, 2);
                Type lagm = growth_parameters(sp, sex, yr - 1, 3);
                length_hat(wtind,  sex, age, yr) = pow(pow(length_hat(wtind,  sex, age - 1, yr - 1), lagm) + (pow(length_hat(wtind,  sex, age - 1, yr - 1), lagm) - pow(laglinf, lagm)) * (exp(-lagkappa) - 1.0), 1 / lagm);
              }
            }
          }
          break;

        case 3: // Non-parametric (Free parameters)
          error("Non-parametric growth not yet implemented");
          // length_hat(wtind,  sex, age, yr) = exp(length_par(sp, sex, age) + length_par_re(sp, sex, age, yr));
          break;

        default:
          error("Invalid 'growth_model");
        } // Growth_model switch


        // 2. Plus-Group Mean Length (Oldest Age Only) ---
        // The plus group holds fish older than the oldest age, so its mean length
        // (cm) sits between L(oldest age) and L-infinity. growth_plus_length picks
        // the form: 1 = M1-weighted, 2/3/4 = SS3's Linf_decay -998/-999/d. See
        // ?build_growth.
        if(growth_model(sp) < 3 && age == (nages(sp) - 1)) {
          Type current_size = length_hat(wtind,  sex, age, yr);
          Type diff = linf - current_size;
          int oldest_age = minage(sp) + nages(sp) - 1;
          if(growth_plus_length(sp) == 1) {
            Type temp_n = 0, temp_sum = 0, weight_a = 1.0;
            Type surv = exp(-exp(log_M1(sp, sex, nages(sp) - 1)));
            for(int a = 0; a <= nages(sp); a++) {
              temp_sum += weight_a * (current_size + (Type(a) / Type(nages(sp))) * diff);
              temp_n += weight_a;
              weight_a *= surv;
            }
            length_hat(wtind,  sex, age, yr) = temp_sum / temp_n;
          } else if(growth_plus_length(sp) == 3) {
            Type temp_n = 0, temp_sum = 0;
            // Survival weights at a fixed 0.2/yr, the nominal mortality SS3 uses
            // for Linf_decay = -999 -- not this species' estimated M1, which is
            // what form 1 above weights by.
            for(int a = 0; a <= oldest_age; a++) {
              Type weight_a = exp(Type(-0.2) * Type(a));
              temp_sum += weight_a * (current_size + (Type(a) / Type(oldest_age)) * diff);
              temp_n += weight_a;
            }
            length_hat(wtind,  sex, age, yr) = temp_sum / temp_n;
          } else if(growth_plus_length(sp) == 4) {
            Type size = current_size, temp_sum = current_size, temp_n = 1.0, weight_a = 1.0;
            // Averaged over twice the oldest age, far enough past it that the
            // decay weights leave no meaningful length behind (SS3's own span).
            for(int a = 1; a <= 2 * oldest_age; a++) {
              weight_a *= exp(-plus_group_decay(sp));
              size += (linf - size) * (Type(1.0) - exp(-kappa));
              temp_sum += weight_a * size;
              temp_n += weight_a;
            }
            length_hat(wtind,  sex, age, yr) = temp_sum / temp_n;
          }
        }

        // 3. SD of length-at-age, then age-length key and weight-at-age ---
        if(growth_model(sp) < 3) {
          Type len = length_hat(wtind,  sex, age, yr);
          Type sd = length_sd_at_age(current_age, age_L1, age == (nages(sp) - 1),
                                     growth_sd_style(sp), growth_sd_form(sp), l1, linf, len,
                                     growth_log_sd(sp, sex, 0), growth_log_sd(sp, sex, 1));
          fill_age_length_key(wtind, sp, sex, age, yr, len, sd, nlengths_pop,
                              lengths_pop, weight_length_pars,
                              mat_len_use, mat_len_pars, growth_matrix, weight_hat, mat_weight_hat);
        }
      } // age
    } // yr
  } // sex
}



/**
 * @brief Integrated Predator Growth, Size-Transition, and Weight-at-Age Module at Month X.
 *
 * @section math_models_within Mathematical Models:
 * 1. Mean Length-at-Age (\f$L_a\f$):
 *    Advances the Jan-1 length stored at id_pop forward by fracyr of growth.
 *    - Von Bertalanffy (Model 1): \f$L_a = L_{\infty} + (L_{a, jan1} - L_{\infty}) \cdot e^{-K \cdot fracyr}\f$
 *    - Richards (Model 2): \f$L_a = [L_{\infty}^m + (L_{a, jan1}^m - L_{\infty}^m) \cdot e^{-K \cdot fracyr}]^{1/m}\f$
 *    - Linear growth is applied for ages < minage (current_age <= age_L1).
 *
 * 2. Weight-at-Age (\f$W_a\f$):
 *    Calculated via integration across the length distribution to account for Jensen's Inequality:
 *    \f$W_a = \sum_{ln} P(ln | a) \cdot \alpha \cdot L_{mid}^{\beta}\f$
 *
 * @section logic_within Biological Logic:
 * - **Temporal Resolution**: Incorporates `fracyr` to allow within-year (seasonal)
 *   growth and differentiability for time-varying parameters.
 * - **Plus Group**: Starts from the Jan-1 plus-group length that
 *   `estimate_growth()` set (`id_pop`). Under `growth_plus_length == 1` it keeps
 *   that length through the year; under the SS3 forms (2-4) it grows within the
 *   year like every other age, as in SS3.
 * - **SD-at-Age**: length_sd_at_age(), as in `estimate_growth()`.
 * - **Size Transition**: fill_age_length_key(), on the population length bins;
 *   consumers bin the product with selectivity, not the key itself.
 *
 * @param wtind Weight index for population/fleet.
 * @param id_pop Index for population (Jan-1) weight-at-age; the within-year
 *               growth advance is anchored on `length_hat(id_pop, ...)`.
 * @param sp Species index.
 * @param fracyr Fraction of the year (0.0 - 1.0) for within-year growth.
 * @param nyrs Number of years in the simulation/model.
 * @param nsex/nages/nlengths Vectors of species-specific counts per dimension.
 * @param lengths Matrix of length bin boundaries (nspp x length_bins). Bin
 *                widths may be non-uniform.
 * @param growth_parameters 4D array of time-varying growth parameters (K, L1, Linf, m).
 * @param growth_log_sd 3D array of log-scale SDs (index 0 = SD at L1, index 1 = SD at Linf).
 * @param weight_length_pars Matrix of length-weight parameters (\f$\alpha\f$, \f$\beta\f$).
 *
 * @note Outputs are passed by reference.
 * @param length_hat [Output] 4D array filled with mean length-at-age.
 * @param growth_matrix [Output] 5D array filled with the growth transition matrix.
 * @param weight_hat [Output] 4D array filled with integrated weight-at-age.
 */
template<class Type>
void estimate_growth_within_yr(
    int wtind,
    int id_pop,
    int sp,
    Type fracyr,
    int nspp,
    int nyrs,
    const vector<int>&  nsex,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  minage,
    const vector<Type>& growth_age_L1,
    const vector<int>&  growth_model,
    const vector<int>&  growth_sd_style,   // Plus-group SD-at-age: 1 = WHAM (pin to exp(sd_Linf)), 2 = interpolate by length
    const vector<int>&  growth_sd_form,    // 1 = SD endpoints in cm, 2 = CV endpoints (SD = CV * L)
    const vector<int>&  growth_plus_length,// 1 = plus group held at its Jan-1 length within the year; 2-4 (SS3 forms) = it grows
    matrix<Type>& lengths,
    const vector<int>&  nlengths_pop,
    matrix<Type>& lengths_pop,             // Population length bins, lower edges (cm)
    array<Type>& growth_parameters,
    array<Type>& growth_log_sd,
    matrix<Type>& weight_length_pars,
    const vector<int>&  mat_len_use,       // 1 = maturity-at-length for this species
    matrix<Type>& mat_len_pars,            // [sp, 0] = L50 (cm), [sp, 1] = logistic slope (per cm)
    array<Type> &length_hat,     // Modified by reference
    array<Type> &growth_matrix,  // Modified by reference
    array<Type> &weight_hat,     // Modified by reference
    array<Type> &mat_weight_hat  // Modified by reference (maturity-at-length species only)
) {

  // Initialize output and temporary storage


  // Calculate mean-length, SD, and growth matrix, for all years:
  // lengths is vector with lengths mm (2, 4, 6, 8, etc)
  Type b_len, current_age, last_linear;

  // The linear ramp below age_L1 starts from the lowest population length
  // edge at age 0 (SS3's len_bins(1)); with no population grid that is the
  // lowest data bin.
  Type Lmin_sp = lengths_pop(sp, 0);
  // age_L1 is the VB anchor age (= age at which `l1` is the length). Read
  // from data_list$growth_age_L1[sp] (= SS3's Growth_Age_for_L1 ctl input).
  // R-side fit_mod() resolves the default to max(0.5, minage[sp]) so models
  // with minage = 0 get an SS3-style half-year anchor and minage >= 1 stays
  // backwards-compatible.
  Type age_L1 = growth_age_L1(sp);

  // Safe denominator for the linear ramp slope when minage = 0. See month=0
  // estimate_growth() for full context.
  Type age_L1_safe = age_L1;
  if (age_L1 <= Type(0)) age_L1_safe = Type(1);

  for(int sex = 0; sex < nsex(sp); sex++) {
    for(int yr = 0; yr < nyrs; yr++) {
      for(int age = 0; age < nages(sp); age++) {

        // Parameters for parametric growth
        Type kappa = growth_parameters(sp, sex, yr, 0);
        Type l1 = growth_parameters(sp, sex, yr, 1);
        Type linf = growth_parameters(sp, sex, yr, 2);
        Type m = growth_parameters(sp, sex, yr, 3);
        // See month=0 overload: slot k = age (k - 1 + minage); current_age is
        // slot age (NOT the VB anchor) offset by fracyr. Was
        // `age + 1.0 + fracyr`, which is off by one when minage != 1.
        current_age = Type(age) + Type(minage(sp)) + fracyr;

        // 1. Calculate Mean Length at Age ---
        switch(growth_model(sp)) {
        case 1: // Von Bertalanffy

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 when minage > 0; == 1 when minage = 0 (ramp branch unreachable)

          // Age < minage
          // The plus group is advanced through within-year growth like every
          // other age (SS3 convention). The recruitment-weighted-mean
          // correction for fish promoted into the plus group is reapplied at
          // every year boundary by estimate_growth(),
          // so id_pop already carries the corrected Jan-1 length.
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          // WHAM-style parameterization: for ages between the anchor and the
          // plus group, blend the linear ramp with the growth curve; the plus
          // group is pinned at its Jan-1 (id_pop) length rather than advanced by
          // within-year growth.
          // TODO(review): this branch (and its Richards mirror below) tests
          // `age + 1.0 < age_L1`, mixing the slot index with an age, while the
          // sibling branch above uses `current_age` (= age + minage + fracyr).
          // For minage != 1 the two disagree, so the interpolated length can be
          // off by one age. Unreachable for minage 0/1; verify for minage >= 2.
          }else if(age + 1.0 < age_L1){ // Linear + growth curve mixed
            last_linear = Lmin_sp + b_len * age_L1;
            length_hat(wtind,  sex, age, yr) = last_linear + (last_linear - linf) * (exp(-kappa * (current_age - age_L1)) - 1.0);
          }else if(age + 1.0 == nages(sp) && growth_plus_length(sp) == 1) { // Plus group held at its Jan-1 length; under the SS3 forms it grows like every other age
            length_hat(wtind,  sex, age, yr) = length_hat(id_pop,  sex, age, yr);
          }else { // Growth curve
            length_hat(wtind,  sex, age, yr) = length_hat(id_pop,  sex, age, yr) + (length_hat(id_pop,  sex, age, yr) - linf) * (exp(-kappa * fracyr) - 1.0); // Add fracyr growth
          }
          break;

        case 2: // Richards

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 when minage > 0; == 1 when minage = 0 (ramp branch unreachable)

          // Mirror the VBGF branch above with the Richards curve so the plus
          // group is treated identically regardless of growth family: a linear
          // ramp below the anchor, a mixed linear+growth branch, the plus group
          // pinned at its Jan-1 (id_pop) length, and within-year Richards growth
          // for the remaining ages. (`last_linear` is declared at function scope.)
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          }else if(age + 1.0 < age_L1){ // Linear + growth curve mixed
            last_linear = Lmin_sp + b_len * age_L1;
            length_hat(wtind,  sex, age, yr) = pow(pow(last_linear, m) + (pow(last_linear, m) - pow(linf, m)) * (exp(-kappa * (current_age - age_L1)) - 1.0), 1 / m);
          }else if(age + 1.0 == nages(sp) && growth_plus_length(sp) == 1) { // Plus group held at its Jan-1 length (see VBGF branch)
            length_hat(wtind,  sex, age, yr) = length_hat(id_pop,  sex, age, yr);
          } else {
            length_hat(wtind,  sex, age, yr) = pow(pow(length_hat(id_pop,  sex, age, yr), m) + (pow(length_hat(id_pop,  sex, age, yr), m) - pow(linf, m)) * (exp(-kappa * fracyr) - 1.0), 1 / m); // Add fracyr growth
          }
          break;

        case 3: // Non-parametric (Free parameters)
          error("Non-parametric growth not yet implemented");
          // length_hat(wtind,  sex, age, yr) = exp(length_par(sp, sex, age) + length_par_re(sp, sex, age, yr));
          break;

        default:
          error("Invalid 'growth_model");
        } // Growth_model switch

        // 2. SD of length-at-age, then age-length key and weight-at-age ---
        if(growth_model(sp) < 3) {
          Type len = length_hat(wtind,  sex, age, yr);
          Type sd = length_sd_at_age(current_age, age_L1, age == (nages(sp) - 1),
                                     growth_sd_style(sp), growth_sd_form(sp), l1, linf, len,
                                     growth_log_sd(sp, sex, 0), growth_log_sd(sp, sex, 1));
          fill_age_length_key(wtind, sp, sex, age, yr, len, sd, nlengths_pop,
                              lengths_pop, weight_length_pars,
                              mat_len_use, mat_len_pars, growth_matrix, weight_hat, mat_weight_hat);
        }
      } // age
    } // yr
  } // sex
}


/**
 * @brief Calculates population and fleet-specific weight-at-age.
 *
 * This function populates the weight_hat array based on either empirical data
 * (growth_model == 0) or estimated growth parameters (growth_model > 0).
 * It handles both hindcast and projection years by carrying over the last
 * hindcast year's empirical data.
 *
 * @param weight_hat [ref] 4D array to store calculated weights (wt_index, sex, age, year)
 * @param length_hat [ref] 4D array for estimated lengths
 * @param growth_matrix [ref] 5D array for age-length transition matrices
 * @param weight_obs Empirical weight-at-age input (wt_index, sex, age, year), kg
 * @param growth_model Integer vector indicating growth type (0=empirical, >0=estimated)
 * @param growth_sd_style Per species, plus-group SD-at-age: 1 = pinned to exp(sd_Linf) (WHAM), 2 = interpolated by length
 * @param growth_sd_form Per species, 1 = SD endpoints in cm, 2 = CV endpoints
 * @param growth_plus_length Per species, plus-group mean length form (see estimate_growth())
 * @param mat_weight_hat [ref] Mature weight-at-age (kg) for maturity-at-length species
 * @param nspp Number of species
 * @param nyrs Total number of years (hindcast + projection)
 * @param nyrs_hind Number of hindcast years
 * @param n_flt Number of fleets
 * @param flt_spp Vector mapping fleet to species index
 * @param flt_month Vector mapping fleet to month of operation
 * @param nsex Vector of number of sexes per species
 * @param nages Vector of number of ages per species
 * @param pop_wt_index Index for population biomass weights
 * @param ssb_wt_index Index for spawning stock biomass weights
 * @param flt_wt_index Index for fleet-specific weights
 * @param spawn_month Vector of spawning months per species
 * [Other parameters for estimate_growth: lengths, nlengths, minage,
 * growth_parameters, growth_log_sd, weight_length_pars]
 */
template <class Type>
void calculate_weight(
    array<Type> &weight_hat,   // Modified by reference
    array<Type> &length_hat,   // Modified by reference
    array<Type> &growth_matrix,// Modified by reference
    array<Type> &mat_weight_hat,// Modified by reference: mature weight-at-age, maturity-at-length species only
    array<Type>& weight_obs,
    const vector<int>&  growth_model,
    const vector<int>&  growth_sd_style,
    const vector<int>&  growth_sd_form,
    const vector<int>&  growth_plus_length,
    const vector<Type>& plus_group_decay,
    const int& nspp,
    const int& nyrs,
    const int& nyrs_hind,
    const int& n_flt,
    const vector<int>&  flt_spp,
    vector<Type> flt_month,
    const vector<int>&  nsex,
    const vector<int>&  minage,
    const vector<Type>& growth_age_L1,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  pop_wt_index,
    const vector<int>&  ssb_wt_index,
    const vector<int>&  flt_wt_index,
    vector<Type> spawn_month,
    matrix<Type>& lengths,
    const vector<int>&  nlengths_pop,
    matrix<Type>& lengths_pop,
    array<Type>& growth_parameters,
    array<Type>& growth_log_sd,
    matrix<Type> weight_length_pars,
    const vector<int>&  mat_len_use,
    matrix<Type>& mat_len_pars,
    array<Type>& log_M1
) {
  int yr_ind;
  int wt_idx_pop;
  int wt_idx_ssb;
  int wt_idx_flt;

  // 1. POPULATION WEIGHT-AT-AGE
  for (int sp = 0; sp < nspp; sp++) {

    wt_idx_pop = 2 * sp;
    wt_idx_ssb = 2 * sp + 1;

    // -- 1.1. Empirical weight-at-age
    if (growth_model(sp) == 0) {
      for (int sex = 0; sex < nsex(sp); sex++) {
        for (int age = 0; age < nages(sp); age++) {
          for (int yr = 0; yr < nyrs; yr++) {

            // Handle projection logic
            yr_ind = (yr < nyrs_hind) ? yr : (nyrs_hind - 1);

            // Biomass weight
            weight_hat(wt_idx_pop, sex, age, yr) = weight_obs(pop_wt_index(sp), sex, age, yr_ind);

            // SSB weight
            weight_hat(wt_idx_ssb, sex, age, yr) = weight_obs(ssb_wt_index(sp), sex, age, yr_ind);
          }
        }
      }
    }

    // -- 1.2. Estimated growth
    if (growth_model(sp) > 0) {
      // Biomass weight (beginning of year / month 0)
      estimate_growth(
        wt_idx_pop,
        sp,
        nspp,
        nyrs,
        nsex,
        nages,
        nlengths,
        minage,
        growth_age_L1,
        growth_model,
        growth_sd_style,
        growth_sd_form,
        growth_plus_length,
        plus_group_decay,
        lengths,
        nlengths_pop,
        lengths_pop,
        growth_parameters,
        growth_log_sd,
        weight_length_pars,
        mat_len_use,
        mat_len_pars,
        log_M1,
        length_hat,     // Pass by reference
        growth_matrix,  // Pass by reference
        weight_hat,     // Pass by reference
        mat_weight_hat  // Pass by reference
      );

      // SSB weight (at month of spawning)
      estimate_growth_within_yr(
        wt_idx_ssb,
        wt_idx_pop,
        sp,
        spawn_month(sp) / Type(12.0),
        nspp,
        nyrs,
        nsex,
        nages,
        nlengths,
        minage,
        growth_age_L1,
        growth_model,
        growth_sd_style,
        growth_sd_form,
        growth_plus_length,
        lengths,
        nlengths_pop,
        lengths_pop,
        growth_parameters,
        growth_log_sd,
        weight_length_pars,
        mat_len_use,
        mat_len_pars,
        length_hat,     // Pass by reference
        growth_matrix,  // Pass by reference
        weight_hat,     // Pass by reference
        mat_weight_hat  // Pass by reference
      );
    }
  }

  // 2. FLEET WEIGHT-AT-AGE
  for (int flt = 0; flt < n_flt; flt++) {
    int sp = flt_spp(flt);
    Type mo = flt_month(flt);
    wt_idx_pop = 2 * sp;
    wt_idx_flt = nspp * 2 + flt;

    // -- 2.1. Empirical weight-at-age
    if (growth_model(sp) == 0) {
      for (int sex = 0; sex < nsex(sp); sex++) {
        for (int age = 0; age < nages(sp); age++) {
          for (int yr = 0; yr < nyrs; yr++) {

            yr_ind = (yr < nyrs_hind) ? yr : (nyrs_hind - 1);
            weight_hat(wt_idx_flt, sex, age, yr) = weight_obs(flt_wt_index(flt), sex, age, yr_ind);
          }
        }
      }
    }

    // -- 2.2. Estimated growth
    if (growth_model(sp) > 0) {

      estimate_growth_within_yr(
        wt_idx_flt,
        wt_idx_pop,
        sp,
        mo / Type(12.0),
        nspp,
        nyrs,
        nsex,
        nages,
        nlengths,
        minage,
        growth_age_L1,
        growth_model,
        growth_sd_style,
        growth_sd_form,
        growth_plus_length,
        lengths,
        nlengths_pop,
        lengths_pop,
        growth_parameters,
        growth_log_sd,
        weight_length_pars,
        mat_len_use,
        mat_len_pars,
        length_hat,     // Pass by reference
        growth_matrix,  // Pass by reference
        weight_hat,     // Pass by reference
        mat_weight_hat  // Pass by reference
      );
    }
  }
}



#endif

