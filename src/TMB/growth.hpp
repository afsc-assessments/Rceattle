#ifndef GROWTH_HPP
#define GROWTH_HPP

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
 * - **Plus-Age Group**: Oldest age class is corrected via a survival-weighted
 *   mean of \f$[current\_size, ..., L_{\infty}]\f$ over \f$a = 0..nages\f$, with
 *   weights \f$\exp(-M_1 a)\f$ at the species' base natural mortality at the oldest
 *   age (`log_M1`). This is the static analogue of SS3's N-at-age-weighted
 *   recruitment correction at the season transition.
 * - **SD-at-Age**: For current_age <= age_L1, SD = \f$e^{sd_0}\f$. Otherwise
 *   linear interpolation in length between SD(\f$l_1\f$) = \f$e^{sd_0}\f$ and
 *   SD(\f$L_{\infty}\f$) = \f$e^{sd_1}\f$. The plus group is pinned to \f$e^{sd_1}\f$
 *   under `growth_sd_style == 1` (WHAM) and interpolated like every other age
 *   under `growth_sd_style == 2` (SS3), as in estimate_growth_within_yr().
 * - **Size Transition**: Converts mean length and SD into a probability
 *   matrix \f$P(\text{Length} | \text{Age})\f$ via `pnorm`. First length bin is
 *   a minus-group; last length bin is a plus-group.
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
    const vector<int>&  growth_sd_style,   // Plus-group SD-at-age: 1 = WHAM (pin to exp(sd_Linf)), 2 = SS3 (interpolate by length)
    matrix<Type>& lengths,
    array<Type>& growth_parameters,
    array<Type>& growth_log_sd,
    matrix<Type>& weight_length_pars,
    array<Type>& log_M1,         // Base natural mortality at age [nspp, nsex, nages], log scale
    array<Type> &length_hat,     // Modified by reference
    array<Type> &growth_matrix,  // Modified by reference
    array<Type> &weight_hat      // Modified by reference
) {

  // Initialize output and temporary storage
  array<Type> length_sd(nsex(sp), nages(sp), nyrs); length_sd.setZero();         // SD in length-at-age


  // Calculate mean-length, SD, and growth matrix, for all years:
  // lengths is vector with lengths mm (2, 4, 6, 8, etc)
  Type Fac1, Fac2, Slope, b_len, last_linear, current_age;

  Type Lmin_sp = lengths(sp, 0);
  Type Lmax_sp = lengths(sp, nlengths(sp) - 1);
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


        // 2. Plus-Group Correction (Oldest Age Only) ---
        // Ages pooled in the plus group, weighted by survival at the oldest-age base
        // M1 (not Z: F and predation are excluded), lengths interpolated to L-infinity.
        if(growth_model(sp) < 3 && age == (nages(sp) - 1)) {
          Type current_size = length_hat(wtind,  sex, age, yr);
          Type temp_n = 0, temp_sum = 0, weight_a = 1.0;
          Type diff = linf - current_size;
          Type surv = exp(-exp(log_M1(sp, sex, nages(sp) - 1)));
          for(int a = 0; a <= nages(sp); a++) {
            temp_sum += weight_a * (current_size + (Type(a) / Type(nages(sp))) * diff);
            temp_n += weight_a;
            weight_a *= surv;
          }
          length_hat(wtind,  sex, age, yr) = temp_sum / temp_n;
        }

        // 3. Calculate SD (Integrated) ---
        // Length-based linear interpolation: SD(l1) = sd0, SD(linf) = sd1.
        // SD-at-age: e^{sd0} up to age_L1, then linear-in-length interpolation
        // to e^{sd1} (the plus group pinned to the upper anchor e^{sd1}).
        if(growth_model(sp) < 3) {
          if((current_age) <= age_L1) {
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0));
          // Plus-group SD-at-age. WHAM (growth_sd_style == 1) pins it to the
          // upper anchor exp(sd_Linf) (SDAA plus group = SD_len(1), matching
          // estimate_growth_within_yr()); SS3 (growth_sd_style == 2) falls
          // through to the length interpolation below, treating the plus group
          // like any interior age.
          } else if(age == (nages(sp) - 1) && growth_sd_style(sp) == 1) {
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 1));
          } else {
            Slope = (exp(growth_log_sd(sp, sex, 1)) - exp(growth_log_sd(sp, sex, 0))) / (linf - l1);
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0)) + Slope * (length_hat(wtind,  sex, age, yr) - l1);
          }

          // Free parameters
          if(growth_model(sp) == 3) {
            // Slope = (exp(growth_log_sd(sp, sex, 1)) - exp(growth_log_sd(sp, sex, 0)))/(length_hat(wtind,  sex, nages(sp)-1, yr) - length_hat(wtind,  sex, 0, yr));
            // length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0) + Slope * (length_hat(wtind,  sex, age, yr) - length_hat(wtind,  sex, 0, yr));
          }
        }

        // 4. Build Growth Matrix & Weight-at-Age simultaneously ---
        Type expected_weight = 0.0;

        for(int ln = 0; ln < nlengths(sp); ln++) {
          Type prob;
          if(ln == 0) {
            Fac1 = (Lmin_sp + lengths(sp, 1) - lengths(sp, 0) - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            prob = pnorm(Fac1);
          } else if(ln == (nlengths(sp) - 1)) {
            Fac1 = (Lmax_sp - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            prob = 1.0 - pnorm(Fac1);
          } else {
            Fac1 = (lengths(sp, ln + 1) - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            Fac2 = (lengths(sp, ln) - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            prob = pnorm(Fac1) - pnorm(Fac2);
          }

          // Explicit assignment to avoid the 5D operator warning
          growth_matrix(wtind, sex, age, ln, yr) = prob;

          // Bin midpoint for weight-at-length. Supports non-uniform bins:
          // each interior bin uses its own [ln, ln+1] midpoint, and the
          // length plus-group extends by half the final bin's width.
          Type lenmid;
          if(ln < nlengths(sp) - 1) {
            lenmid = (lengths(sp, ln) + lengths(sp, ln + 1)) / Type(2.0);
          } else {
            Type last_width = lengths(sp, nlengths(sp) - 1) - lengths(sp, nlengths(sp) - 2);
            lenmid = lengths(sp, nlengths(sp) - 1) + last_width / Type(2.0);
          }

          // Weighted sum for Weight-at-Age
          expected_weight += prob * weight_length_pars(sp, 0) * pow(lenmid, weight_length_pars(sp, 1));
        }
        weight_hat(wtind, sex, age, yr) = expected_weight;
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
 * - **Plus Group**: Advanced by within-year growth identically to other ages
 *   (SS3 convention). The recruitment-weighted-mean correction for fish promoted
 *   into the plus group is applied at every year boundary by `estimate_growth()`
 *   (month 0), so `id_pop` already carries the corrected Jan-1 length.
 * - **SD-at-Age**: Length-based linear interpolation between `SD(l1)` and
 *   `SD(linf)` above `age_L1`. The plus group is pinned to `SD(linf)` =
 *   `exp(sd_Linf)` under `growth_sd_style == 1` (WHAM) and interpolated like
 *   every other age under `growth_sd_style == 2` (SS3), as in `estimate_growth()`.
 * - **Size Transition**: Converts mean length and SD into a probability matrix
 *   \f$P(\text{Length} | \text{Age})\f$ using a cumulative normal distribution (`pnorm`).
 *   The first length bin is a minus-group on length; the last is a plus-group.
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
    const vector<int>&  growth_sd_style,   // Plus-group SD-at-age: 1 = WHAM (pin to exp(sd_Linf)), 2 = SS3 (interpolate by length)
    matrix<Type>& lengths,
    array<Type>& growth_parameters,
    array<Type>& growth_log_sd,
    matrix<Type>& weight_length_pars,
    array<Type> &length_hat,     // Modified by reference
    array<Type> &growth_matrix,  // Modified by reference
    array<Type> &weight_hat      // Modified by reference
) {

  // Initialize output and temporary storage
  array<Type> length_sd(nsex(sp), nages(sp), nyrs); length_sd.setZero();         // SD in length-at-age


  // Calculate mean-length, SD, and growth matrix, for all years:
  // lengths is vector with lengths mm (2, 4, 6, 8, etc)
  Type Fac1, Fac2, Slope, b_len, current_age, last_linear;

  Type Lmin_sp = lengths(sp, 0);
  Type Lmax_sp = lengths(sp, nlengths(sp) - 1);
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
          }else if(age + 1.0 == nages(sp)) { // Plus group pinned at Jan-1. Comment out for SS3 style.
            length_hat(wtind,  sex, age, yr) = length_hat(id_pop,  sex, age, yr);
          }else { // Growth curve (excl. plus group)
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
          }else if(age + 1.0 == nages(sp)) { // Plus group pinned at Jan-1. Comment out for SS3 style.
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

        // 2. Calculate SD (Integrated) ---
        // SD-at-age: e^{sd0} up to age_L1, then linear-in-length interpolation
        // to e^{sd1} (the plus group pinned to the upper anchor e^{sd1}).
        if(growth_model(sp) < 3) {
          if((current_age) <= age_L1) {
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0));
          // Plus-group SD-at-age. WHAM (growth_sd_style == 1) pins it to the
          // upper anchor exp(sd_Linf) (SDAA plus group = SD_len(1), matching
          // estimate_growth()); SS3 (growth_sd_style == 2) falls through to the
          // length interpolation below, treating the plus group like any
          // interior age.
          } else if(age == (nages(sp) - 1) && growth_sd_style(sp) == 1) {
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 1));
          } else {
            Slope = (exp(growth_log_sd(sp, sex, 1)) - exp(growth_log_sd(sp, sex, 0))) / (linf - l1);
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0)) + Slope * (length_hat(wtind,  sex, age, yr) - l1);
          }

          // Free parameters
          if(growth_model(sp) == 3) {
            // Slope = (exp(growth_log_sd(sp, sex, 1)) - exp(growth_log_sd(sp, sex, 0)))/(length_hat(wtind,  sex, nages(sp)-1, yr) - length_hat(wtind,  sex, 0, yr));
            // length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0) + Slope * (length_hat(wtind,  sex, age, yr) - length_hat(wtind,  sex, 0, yr));
          }
        }

        // 4. Build Growth Matrix & Weight-at-Age simultaneously ---
        Type expected_weight = 0.0;

        for(int ln = 0; ln < nlengths(sp); ln++) {
          Type prob;
          if(ln == 0) {
            Fac1 = (Lmin_sp + lengths(sp, 1) - lengths(sp, 0) - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            prob = pnorm(Fac1);
          } else if(ln == (nlengths(sp) - 1)) {
            Fac1 = (Lmax_sp - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            prob = 1.0 - pnorm(Fac1);
          } else {
            Fac1 = (lengths(sp, ln + 1) - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            Fac2 = (lengths(sp, ln) - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            prob = pnorm(Fac1) - pnorm(Fac2);
          }

          // Explicit assignment to avoid the 5D operator warning
          growth_matrix(wtind, sex, age, ln, yr) = prob;

          // Bin midpoint for weight-at-length. Supports non-uniform bins:
          // each interior bin uses its own [ln, ln+1] midpoint, and the
          // length plus-group extends by half the final bin's width.
          Type lenmid;
          if(ln < nlengths(sp) - 1) {
            lenmid = (lengths(sp, ln) + lengths(sp, ln + 1)) / Type(2.0);
          } else {
            Type last_width = lengths(sp, nlengths(sp) - 1) - lengths(sp, nlengths(sp) - 2);
            lenmid = lengths(sp, nlengths(sp) - 1) + last_width / Type(2.0);
          }

          // Weighted sum for Weight-at-Age
          expected_weight += prob * weight_length_pars(sp, 0) * pow(lenmid, weight_length_pars(sp, 1));
        }
        weight_hat(wtind, sex, age, yr) = expected_weight;
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
 * @param growth_sd_style Per species, plus-group SD-at-age: 1 = pinned to exp(sd_Linf) (WHAM), 2 = interpolated by length (SS3)
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
    array<Type>& weight_obs,
    const vector<int>&  growth_model,
    const vector<int>&  growth_sd_style,
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
    array<Type>& growth_parameters,
    array<Type>& growth_log_sd,
    matrix<Type> weight_length_pars,
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
        lengths,
        growth_parameters,
        growth_log_sd,
        weight_length_pars,
        log_M1,
        length_hat,     // Pass by reference
        growth_matrix,  // Pass by reference
        weight_hat      // Pass by reference
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
        lengths,
        growth_parameters,
        growth_log_sd,
        weight_length_pars,
        length_hat,     // Pass by reference
        growth_matrix,  // Pass by reference
        weight_hat      // Pass by reference
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
        lengths,
        growth_parameters,
        growth_log_sd,
        weight_length_pars,
        length_hat,     // Pass by reference
        growth_matrix,  // Pass by reference
        weight_hat      // Pass by reference
      );
    }
  }
}



// ------------------------------------------------------------------------- //
// TODO                                                                      //
// ------------------------------------------------------------------------- //
// add growth parameters to build param and build map (AR1 and variance)
// define age_L1 and age_L1 ceiling

#endif

