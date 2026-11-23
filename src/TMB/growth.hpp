#ifndef GROWTH_HPP
#define GROWTH_HPP

/**
 * @brief SD of length-at-age (cm) for one age, year and sex.
 *
 * The variability term v is \f$e^{sd_0}\f$ at or below `age_L1`, pinned to
 * \f$e^{sd_1}\f$ for the plus group under `sd_style == 1` (WHAM; also SS3 with
 * Growth_Age_for_L2 = 999), and otherwise interpolated linearly in length between
 * \f$(l_1, e^{sd_0})\f$ and \f$(L_\infty, e^{sd_1})\f$. `sd_form` 1 makes v an SD in cm
 * (SS3 CV_Growth_Pattern 2), 2 a CV, so SD = v * L (SS3 CV_Growth_Pattern 0).
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
 * Length is normal with mean `mu` and SD `sd` (cm) on the population length bins
 * `lengths_pop` (lower edges): the first bin is a minus group, the last a plus group
 * (SS3's convention).
 *
 * Weight-at-age is \f$\sum_l P(l|a) \alpha L_{mid}^\beta\f$ at population-bin midpoints
 * (kg). With maturity-at-length, mature weight-at-age is the same sum with the logistic
 * maturity \f$1/(1+e^{-s(L_{mid} - L_{50})})\f$ inside: SS3's fecundity-at-age when
 * fecundity equals body weight.
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
  // Kept on the POPULATION bins: readers multiply by selectivity-at-length there, as SS3
  // does, and sum into data bins afterwards (ceattle.cpp section 2.3c). Aggregating here
  // would apply one bin-average selectivity across a data bin.
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

    // Bin midpoint (cm): lower edge plus half its own width, the last bin taking the
    // previous bin's width. A single-bin grid has no neighbour to take a width from,
    // so it keeps the edge itself -- as selectivity.hpp does.
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
 * @brief Jan-1 (month 0) mean length, age-length key and weight-at-age for one species.
 *
 * Mean length-at-age (cm) is anchored at \f$L(a_{L1}) = l_1\f$:
 * - von Bertalanffy (`growth_model` 1):
 *   \f$L_a = L_{\infty} + (l_1 - L_{\infty}) \cdot e^{-K(a - a_{L1})}\f$
 * - Richards (2):
 *   \f$L_a = [L_{\infty}^m + (l_1^m - L_{\infty}^m) \cdot e^{-K(a - a_{L1})}]^{1/m}\f$
 * - Below the anchor, a linear ramp from `Lmin_sp` at age 0 to \f$l_1\f$ at `age_L1`.
 * - After year 0 a cohort recursion advances \f$L_{a, y}\f$ from \f$L_{a-1, y-1}\f$ with
 *   the previous year's parameters; where that predecessor is still on the ramp, or at
 *   slot 0, the closed-form anchor at \f$l_1\f$ is used instead.
 * Adapted from Giancarlo Correa's mean length-at-age block in src/wham_v0.cpp, WHAM's
 * growth branch (timjmiller/wham at 16770d4, GPL-3).
 *
 * Weight-at-age (kg) is integrated across the length distribution rather than taken at
 * mean length, which Jensen's inequality would bias:
 * \f$W_a = \sum_{ln} P(ln | a) \cdot \alpha \cdot L_{mid}^{\beta}\f$.
 *
 * @param wtind Weight index slot to write into.
 * @param sp Species index.
 * @param minage Minimum modelled age per species; 0 is supported.
 * @param growth_model Per species, 1 = von Bertalanffy, 2 = Richards.
 * @param growth_parameters Time-varying growth parameters by species, sex, year and
 *   parameter: K (per year), L1 and Linf (cm), Richards m.
 * @param growth_log_sd Log-scale variability endpoints: index 0 at L1, index 1 at Linf.
 * @param weight_length_pars Allometric \f$\alpha\f$ and \f$\beta\f$, kg from cm.
 * @param length_hat [Output] mean length-at-age (cm).
 * @param growth_matrix [Output] age-length key on the population length bins.
 * @param weight_hat [Output] integrated weight-at-age (kg).
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
    const vector<int>&  growth_sd_style,   // Plus-group SD-at-age: 1 = WHAM pin, 2 = SS3 by length
    const vector<int>&  growth_sd_form,    // 1 = SD endpoints in cm, 2 = CV endpoints (SD = CV * L)
    const vector<int>&  growth_plus_length,// Plus-group mean length: 1 = M1-weighted,
                                           // 2 = no adjustment, 3/4 = SS3 Linf_decay
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

  Type b_len, last_linear, current_age;

  // The linear ramp below age_L1 starts from the lowest population length edge (cm) at
  // age 0 (SS3's len_bins(1)); with no population grid that is the lowest data bin.
  Type Lmin_sp = lengths_pop(sp, 0);
  // Anchor age, where mean length is `l1` (SS3's Growth_Age_for_L1). fit_mod() resolves
  // an unset value to max(0.5, minage[sp]): an SS3-style half-year anchor at minage = 0.
  Type age_L1 = growth_age_L1(sp);

  // An anchor at age 0 would divide by zero in the ramp slope b_len, so the slope and the
  // cohort-boundary test below both floor it at 1 year.
  Type age_L1_safe = age_L1;
  if (age_L1 <= Type(0)) age_L1_safe = Type(1);

  for(int sex = 0; sex < nsex(sp); sex++) {
    for(int yr = 0; yr < nyrs; yr++) {
      for(int age = 0; age < nages(sp); age++) {

        // This sex and year: K (per year), L1 and Linf (cm), Richards m.
        Type kappa = growth_parameters(sp, sex, yr, 0);
        Type l1 = growth_parameters(sp, sex, yr, 1);
        Type linf = growth_parameters(sp, sex, yr, 2);
        Type m = growth_parameters(sp, sex, yr, 3);
        // Age (years) at the start of the year: 0-based slot `age` holds the
        // (age + minage) cohort. This is the SLOT age, not the anchor age_L1 above, and
        // adding 1 rather than minage reads the wrong age at any minage != 1.
        current_age = Type(age) + Type(minage(sp));

        // 1. Calculate Mean Length at Age ---
        switch(growth_model(sp)) {
        case 1: // Von Bertalanffy

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 unless the anchor is <= 0

          // Linear ramp below the anchor age
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else {
            if(yr == 0) {
              length_hat(wtind,  sex, age, yr) = linf + (l1 - linf) * exp(-kappa * (current_age - age_L1));
            } else { // Yr > 0
              // The recursion below reads the same cohort one year younger, so it is
              // valid only once that slot is on the curve. Where it is still on the ramp
              // -- or at slot 0, which has no predecessor -- read off the anchor instead.
              if(current_age - Type(1) < age_L1_safe || age == 0) {
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
          b_len = (l1 - Lmin_sp) / age_L1_safe;

          // Linear ramp below the anchor age
          if((current_age) <= age_L1) {
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else {
            if(yr == 0) {
              length_hat(wtind,  sex, age, yr) = pow(pow(linf, m) + (pow(l1, m) - pow(linf, m)) * exp(-kappa * (current_age - age_L1)), 1 / m);
            } else { // Yr > 0
              // Same cohort-boundary test as the VBGF branch above.
              if(current_age - Type(1) < age_L1_safe || age == 0) {
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
          break;

        default:
          error("Invalid 'growth_model");
        } // Growth_model switch


        // 2. Plus-Group Mean Length (Oldest Age Only) ---
        // The plus group holds fish past the oldest age, so its mean length (cm) lies
        // between L(oldest age) and Linf. growth_plus_length: 1 = weighted by survival at
        // the oldest age's M1, 2/3/4 = SS3's Linf_decay -998/-999/d (?build_growth).
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
            // Survival weights at a fixed 0.2/yr, the nominal mortality SS3 uses for
            // Linf_decay = -999 -- not this species' M1, which is what form 1 above uses.
            for(int a = 0; a <= oldest_age; a++) {
              Type weight_a = exp(Type(-0.2) * Type(a));
              temp_sum += weight_a * (current_size + (Type(a) / Type(oldest_age)) * diff);
              temp_n += weight_a;
            }
            length_hat(wtind,  sex, age, yr) = temp_sum / temp_n;
          } else if(growth_plus_length(sp) == 4) {
            Type size = current_size, temp_sum = current_size, temp_n = 1.0, weight_a = 1.0;
            // Each further age is one more year along the von Bertalanffy curve, for
            // either growth family, weighted by exp(-plus_group_decay) per year. Averaged
            // over 2 x the oldest age, SS3's span, by when the weights are negligible.
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
 * @brief Mean length-at-age, age-length key and weight-at-age `fracyr` into the year.
 *
 * Advances the Jan-1 mean length stored at `id_pop` by `fracyr` of growth, for SSB
 * weight-at-age at the spawning month and each fleet's weight-at-age at its own month.
 * Adapted from Giancarlo Correa's pred_LAA() in the same WHAM growth branch:
 * - von Bertalanffy (`growth_model` 1):
 *   \f$L_a = L_{\infty} + (L_{a, jan1} - L_{\infty}) \cdot e^{-K \cdot fracyr}\f$
 * - Richards (2):
 *   \f$L_a = [L_{\infty}^m + (L_{a, jan1}^m - L_{\infty}^m) \cdot e^{-K \cdot fracyr}]^{1/m}\f$
 * - Below the anchor age (`current_age <= age_L1`) the linear ramp applies, as at Jan-1.
 *
 * Under `growth_plus_length == 1` the plus group keeps the Jan-1 mean length
 * `estimate_growth()` set for it; under the SS3 forms (2-4) it grows within the year like
 * every other age, as SS3 does.
 *
 * Weight-at-age (kg) is integrated over the length distribution as in estimate_growth().
 *
 * @param wtind Weight index slot to write into.
 * @param id_pop Weight index holding the Jan-1 lengths the advance starts from.
 * @param sp Species index.
 * @param fracyr Fraction of the year (0.0 - 1.0) to advance growth by.
 * @param growth_parameters,growth_log_sd,weight_length_pars As in estimate_growth().
 * @param length_hat,growth_matrix,weight_hat [Output] As in estimate_growth().
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
    const vector<int>&  growth_sd_style,   // Plus-group SD-at-age: 1 = WHAM pin, 2 = SS3 by length
    const vector<int>&  growth_sd_form,    // 1 = SD endpoints in cm, 2 = CV endpoints (SD = CV * L)
    const vector<int>&  growth_plus_length,// 1 = plus group held at Jan-1 length; 2-4 = it grows
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

  Type b_len, current_age, last_linear;

  // Ramp origin (lowest population length edge, cm, SS3's len_bins(1)) and anchor age
  // (SS3's Growth_Age_for_L1), both as in estimate_growth().
  Type Lmin_sp = lengths_pop(sp, 0);
  Type age_L1 = growth_age_L1(sp);

  // Floors the ramp slope's denominator at 1 year when the anchor is at age 0. This
  // routine advances each slot from its own Jan-1 length, so unlike estimate_growth() it
  // has no cohort-boundary test to floor.
  Type age_L1_safe = age_L1;
  if (age_L1 <= Type(0)) age_L1_safe = Type(1);

  for(int sex = 0; sex < nsex(sp); sex++) {
    for(int yr = 0; yr < nyrs; yr++) {
      for(int age = 0; age < nages(sp); age++) {

        // This sex and year: K (per year), L1 and Linf (cm), Richards m.
        Type kappa = growth_parameters(sp, sex, yr, 0);
        Type l1 = growth_parameters(sp, sex, yr, 1);
        Type linf = growth_parameters(sp, sex, yr, 2);
        Type m = growth_parameters(sp, sex, yr, 3);
        // As in estimate_growth(): the slot age, not the anchor age_L1, offset by fracyr.
        // Adding 1 rather than minage reads the wrong age at any minage != 1.
        current_age = Type(age) + Type(minage(sp)) + fracyr;

        // 1. Calculate Mean Length at Age ---
        switch(growth_model(sp)) {
        case 1: // Von Bertalanffy

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 unless the anchor is <= 0

          // Linear ramp below the anchor age
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          // Plus group PINNED at the Jan-1 length estimate_growth() wrote for it, so this
          // override comes BEFORE the anchor test below, which also matches it when the
          // whole age range sits under the anchor and would discard that pinned length.
          // `age + 1.0 == nages(sp)` is a slot ordinal (0-based `age`, `nages` a count).
          }else if(age + 1.0 == nages(sp) && growth_plus_length(sp) == 1) { // Pinned at Jan-1
            length_hat(wtind,  sex, age, yr) = length_hat(id_pop,  sex, age, yr);
          // The cohort that CROSSED the anchor during the year -- Jan-1 age below age_L1,
          // mid-year age above it -- is read off the growth curve from the anchor, not
          // advanced from its ramp length. `last_linear` is identically `l1` here, so the
          // body is linf + (l1 - linf) * exp(-kappa * (current_age - age_L1)). The test
          // is the JAN-1 age: `current_age` carries fracyr and the first branch took it.
          }else if(age + Type(minage(sp)) < age_L1){ // Linear + growth curve mixed
            last_linear = Lmin_sp + b_len * age_L1;
            length_hat(wtind,  sex, age, yr) = last_linear + (last_linear - linf) * (exp(-kappa * (current_age - age_L1)) - 1.0);
          }else { // Growth curve
            length_hat(wtind,  sex, age, yr) = length_hat(id_pop,  sex, age, yr) + (length_hat(id_pop,  sex, age, yr) - linf) * (exp(-kappa * fracyr) - 1.0); // Add fracyr growth
          }
          break;

        case 2: // Richards

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;

          // Same branch order as the VBGF case above, with the Richards curve, so the
          // plus group is treated identically whichever growth family a species uses.
          // Linear ramp below the anchor age
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          // Plus group PINNED at its Jan-1 length, ordered before the anchor test for the
          // reason given in the VBGF branch above.
          }else if(age + 1.0 == nages(sp) && growth_plus_length(sp) == 1) { // Pinned at Jan-1
            length_hat(wtind,  sex, age, yr) = length_hat(id_pop,  sex, age, yr);
          // As in the VBGF branch, but read off the Richards curve: `last_linear` is
          // identically `l1` here, so the body is
          // [linf^m + (l1^m - linf^m) * exp(-kappa * (current_age - age_L1))]^(1/m).
          }else if(age + Type(minage(sp)) < age_L1){ // Linear + growth curve mixed
            last_linear = Lmin_sp + b_len * age_L1;
            length_hat(wtind,  sex, age, yr) = pow(pow(last_linear, m) + (pow(last_linear, m) - pow(linf, m)) * (exp(-kappa * (current_age - age_L1)) - 1.0), 1 / m);
          } else {
            length_hat(wtind,  sex, age, yr) = pow(pow(length_hat(id_pop,  sex, age, yr), m) + (pow(length_hat(id_pop,  sex, age, yr), m) - pow(linf, m)) * (exp(-kappa * fracyr) - 1.0), 1 / m); // Add fracyr growth
          }
          break;

        case 3: // Non-parametric (Free parameters)
          error("Non-parametric growth not yet implemented");
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
 * @brief Population, SSB and fleet weight-at-age (kg) for every species and fleet.
 *
 * Empirical weight-at-age (`growth_model == 0`) is read straight from `weight_obs`, each
 * projection year carrying the last hindcast year forward. Under estimated growth
 * (`growth_model > 0`) population weight comes from estimate_growth() at Jan-1, and SSB
 * and fleet weight from estimate_growth_within_yr() at the spawning month and at each
 * fleet's own month. Arguments not listed below pass through to those two, except n_flt,
 * which this function's own fleet loop consumes.
 *
 * @param weight_hat [ref] weight-at-age (wt_index, sex, age, year), kg
 * @param length_hat [ref] mean length-at-age (cm)
 * @param growth_matrix [ref] age-length key on the population length bins
 * @param mat_weight_hat [ref] mature weight-at-age (kg), maturity-at-length species only
 * @param weight_obs Empirical weight-at-age input (wt_index, sex, age, year), kg
 * @param growth_model Per species, 0 = empirical, 1 = von Bertalanffy, 2 = Richards
 * @param nyrs Years, hindcast plus projection.
 * @param nyrs_hind Years in the hindcast.
 * @param flt_spp Species index of each fleet
 * @param flt_month Fleet observation month (0 = not specified); flt_month/12 of a year
 * @param pop_wt_index Weight-at-age input index for population biomass
 * @param ssb_wt_index Weight-at-age input index for spawning-stock biomass
 * @param flt_wt_index Weight-at-age input index per fleet
 * @param spawn_month Spawning month per species; spawn_month/12 of a year (0 = Jan-1)
 */
template <class Type>
void calculate_weight(
    array<Type> &weight_hat,   // Modified by reference
    array<Type> &length_hat,   // Modified by reference
    array<Type> &growth_matrix,// Modified by reference
    array<Type> &mat_weight_hat,// Modified by reference; maturity-at-length species only
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

            // Projection years reuse the last hindcast year's empirical weight-at-age.
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
        length_hat,
        growth_matrix,
        weight_hat,
        mat_weight_hat
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
        length_hat,
        growth_matrix,
        weight_hat,
        mat_weight_hat
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
        length_hat,
        growth_matrix,
        weight_hat,
        mat_weight_hat
      );
    }
  }
}



#endif

