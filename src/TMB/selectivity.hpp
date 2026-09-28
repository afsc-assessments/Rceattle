#ifndef SELECTIVITY_HPP
#define SELECTIVITY_HPP

/**
 * @brief Performs final normalization and projection of fishery selectivity (age/length) across all fleets.
 *
 * @details This function iterates through all fleets to process selectivity data.
 * It handles:
 * 1. Zeroing out selectivity for ages below bin_first_selected.
 * 2. Normalizing selectivity values based on a single age (sel_norm_bin1 >= 0),
 * an age range average (sel_norm_bin2 >= 0), or the maximum value
 * across all ages/sexes (sel_norm_bin1 between -1 and -499).
 * 3. Projecting the final year of hindcast selectivity into future projection years.
 *
 * @param flt fleet to iterate through.
 * @param nyrs_hind Total number of years in the hindcast period.
 * @param nyrs Total number of years including projection.
 * @param flt_spp Function/Array mapping fleet index to species index.
 * @param flt_sel_type Function/Array mapping fleet index to selectivity type.
 * @param flt_sel_dim Age or length based selectivity.
 * @param bin_first_selected Array/function returning the minimum age/length bin of selection for a fleet.
 * @param nages Array/function returning the number of ages for a species.
 * @param nlengths Array/function returning the number of length bins for a species.
 * @param nsex Array/function returning the number of sexes for a species.
 * @param sel_norm_bin1 Array/function returning the normalization age or control flag.
 * @param sel_norm_bin2 Array/function returning the upper bound for age-range normalization.
 * @param sel_norm_scope Per fleet: 0 = WithinSex (each sex divided by its own
 *   reference, so both reach 1), 1 = AcrossSexes (one reference pooled over the
 *   sexes, so the less-selected sex stays below 1). Orthogonal to
 *   sel_norm_bin1/2, which say WHERE the reference is taken.
 * @param selectivity 4D container (fleet, sex, age or length bin, year) modified in-place.
 */
template<class Type>
void normalize_and_project_selectivity(
    const int& flt,
    const int& nyrs_hind,
    const int& nyrs,
    const vector<int>&  flt_spp,
    const vector<int>&  flt_sel_type,
    const vector<int>&  flt_sel_dim,
    const vector<int>&  bin_first_selected,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  nsex,
    const vector<int>&  sel_norm_bin1,
    const vector<int>&  sel_norm_bin2,
    const vector<int>&  sel_norm_scope,
    array<Type> &selectivity // Modified by reference
) {
  int sp = flt_spp(flt);
  int sel_type = flt_sel_type(flt);
  int nbins = (flt_sel_dim(flt) == 0) ? nages(sp) : nlengths(sp);

  // Ages not selected
  for(int yr = 0; yr < nyrs_hind; yr++) {
    for(int bin = 0; bin < nbins; bin++){
      for(int sex = 0; sex < nsex(sp); sex++){
        if(bin < bin_first_selected(flt)){
          selectivity(flt, sex, bin, yr) = 0;
        }
      }
    }
  }

  // Normalization makes two orthogonal choices, carried in two inputs rather than
  // encoded in the sign of one:
  //   sel_norm_bin1/2 -- WHERE the reference is taken: a named bin, the mean over
  //     a bin range, or (< 0) the maximum over bins;
  //   sel_norm_scope  -- WHOSE scale it sets: WithinSex (0) gives each sex its own
  //     reference so both reach 1; AcrossSexes (1) pools one reference so the
  //     less-selected sex stays below 1 and relative sex selectivity survives.
  // Identical for a one-sex species. Hake (5/12) normalizes in its own year/sex
  // block above; LogisticPM (11) reuses sel_norm_bin1/2 as a penalty age-range.
  if((sel_norm_bin1(flt) > -500) && (sel_type != 5) && (sel_type != 12) && (sel_type != 11)) {
    bool at_bin     = (sel_norm_bin1(flt) >= 0);
    bool over_range = at_bin && (sel_norm_bin2(flt) >= 0);
    bool across     = (sel_norm_scope(flt) == 1);

    vector<Type> ref(nsex(sp));

    for(int yr = 0; yr < nyrs_hind; yr++) {

      // Reference value for each sex.
      for(int sex = 0; sex < nsex(sp); sex++){
        if(over_range){
          // Mean selectivity over the bin range
          ref(sex) = 0;
          for(int bin = sel_norm_bin1(flt); bin <= sel_norm_bin2(flt); bin++) {
            ref(sex) += selectivity(flt, sex, bin, yr)/(sel_norm_bin2(flt) - sel_norm_bin1(flt) + 1);
          }
        } else if(at_bin){
          // Single normalization bin, floored so the division stays finite
          ref(sex) = max2(Type(0.001), selectivity(flt, sex, sel_norm_bin1(flt), yr));
        } else {
          // Maximum over bins. Fold rather than branch: TMB tapes once, so an
          // `if` on an AD Type would freeze the argmax bin at its initial
          // position.
          ref(sex) = 0;
          for(int bin = 0; bin < nbins; bin++){
            ref(sex) = max2(ref(sex), selectivity(flt, sex, bin, yr));
          }
        }
      }

      // One divisor for every sex, so their relative levels survive.
      if(across){
        Type pooled = ref(0);
        for(int sex = 1; sex < nsex(sp); sex++) pooled = max2(pooled, ref(sex));
        for(int sex = 0; sex < nsex(sp); sex++) ref(sex) = pooled;
      }

      for(int bin = 0; bin < nbins; bin++){
        for(int sex = 0; sex < nsex(sp); sex++){
          selectivity(flt, sex, bin, yr) /= ref(sex);
        }
      }
    }
  }

  // Project outwards assuming final years selectivity
  for(int yr = nyrs_hind; yr < nyrs; yr++) {
    for(int bin = 0; bin < nbins; bin++){
      for(int sex = 0; sex < nsex(sp); sex++){
        selectivity(flt, sex, bin, yr) = selectivity(flt, sex, bin, nyrs_hind - 1);
      }
    }
  }
}


/**
 * @brief Converts length-based selectivity to age-based selectivity and stores it in a 4D array.
 * * @details This function maps selectivity from length bins to age classes by
 * integrating over a growth transition matrix (Prob(Length | Age)). The resulting
 * age-based selectivity is stored directly into the provided 'sel_at_age' array for the
 * specified fleet, sex, and year using pass by reference.
 *
 * @tparam Type The numeric type (e.g., double, TMB::adouble).
 * @param flt Index of the fleet/fishery.
 * @param sp Index of the species.
 * @param sex Index of the sex.
 * @param yr Index of the year.
 * @param wtind Index for the weight-at-age/growth transition matrix type.
 * @param nages Reference to vector containing number of ages per species.
 * @param nlengths Reference to vector containing number of length bins per species.
 * @param sel_length Reference to 4D array of selectivity-at-length [flt, sex, ln, yr].
 * @param growth_matrix Reference to 5D transition matrix [wtind, sex, age, ln, yr].
 * @param sel_at_age Reference to 4D output array [flt, sex, age, yr] to be updated in-place.
 */
template <class Type>
void convert_length_selectivity(
    const int& flt,
    const int& sp,
    const int& sex,
    const int& yr,
    const int& wtind,
    const vector<int>& nages,
    const vector<int>& nlengths,
    array<Type>& sel_length,
    array<Type>& growth_matrix,
    array<Type>& sel_at_age // Modified by reference
) {
  // Iterate through ages for the specific species
  for(int age = 0; age < nages(sp); age++) {

    // Initialize the specific cell in the 4D array to zero before accumulating
    sel_at_age(flt, sex, age, yr) = 0.0;

    for(int ln = 0; ln < nlengths(sp); ln++) {
      // sel_at_age = sum over lengths of (Prob(Length|Age) * sel_at_length)
      sel_at_age(flt, sex, age, yr) += growth_matrix(wtind, sex, age, ln, yr) * sel_length(flt, sex, ln, yr);
    }
  }
}


/**
 * @brief Unified function to calculate selectivity across fleets, sexes, bins (age or length), and years.
 * * @details This function calculates selectivity for the hindcast period,
 * normalizes it, projects it into the future, and—if the base calculation was length-based—uses a
 * growth transition matrix to output the equivalent age-based selectivity.
 * * @section functional_forms Selectivity Models:
 * The dimension (age vs length) is determined by `flt_sel_dim` (0 = age, 1 = length) and is shared
 * across all parametric cases; the age/length distinction is handled inside each case, not by a
 * separate case number. After calculation, length-based selectivity is converted to age-based via
 * `convert_length_selectivity()`, and all selectivities are normalized/projected via
 * `normalize_and_project_selectivity()`.
 *
 * - Case 0: Empirical — fixed selectivity read from `emp_sel` input; assumed age-based.
 *
 * - Case 1: Logistic — 2-parameter ascending logistic [Age or Length]
 *     sel(x) = 1 / (1 + exp(-slope * (x - inflection)))
 *     Parameters: log_sel_slp(0) = log(slope); sel_inf(0) = inflection point
 *
 * - Case 2: Non-parametric (Ianelli/AMAK style) [Age or Length]
 *     Bin-specific log-scale coefficients, normalized by mean; penalized for monotonicity and curvature.
 *     Parameters: sel_coff[flt, sex, bin] for bins 0..N_sel_bins-1
 *
 * - Case 3: Double Logistic (Dorn and Methot 1990) [Age or Length]
 *     sel(x) = asc(x) * (1 - desc(x)), where each limb is a 2-parameter logistic.
 *     Parameters: log_sel_slp(0/1) = log(ascending/descending slope); sel_inf(0/1) = ascending/descending inflection
 *
 * - Case 4: Descending Logistic [Age or Length]
 *     sel(x) = 1 - 1/(1 + exp(-slope * (x - inflection)))
 *     Parameters: log_sel_slp(1) = log(slope); sel_inf(1) = inflection point
 *
 * - Case 5: Non-parametric (Hake/Taylor style) [Age or Length]
 *     Cumulative sum of bin coefficients, exponentiated and normalized to max = 1.
 *     Parameters: sel_coff[flt, sex, bin] for bins bin_first_selected..N_sel_bins-1
 *
 * - Case 6: 2D AR1 (age x year) [Age or Length]
 *     Logistic-transformed bin-year deviates with AR1 correlation across age and year.
 *
 * - Case 7: 3D AR1 (age x cohort x year, Cheng et al 2024) [Age or Length]
 *     Conditional-variance parameterisation of a 3D random field across age, cohort and year.
 *
 * - Case 8: Double Normal [Age or Length]
 *     Four-parameter dome with a configurable right-tail floor (analogous to SS3 pattern
 *     24 / P6). Gaussian ascending limb left of peak; Gaussian descending limb right of
 *     peak that approaches `right_floor` instead of 0; blended smoothly via a steep
 *     logistic weight so the function is twice differentiable everywhere.
 *       w(x)         = 1 / (1 + exp(-20*(x - peak)))     [~0 left of peak, ~1 right]
 *       asc(x)       = exp(-0.5*((x-peak)/sigma_asc)^2)
 *       desc(x)      = right_floor + (1 - right_floor) * exp(-0.5*((x-peak)/sigma_desc)^2)
 *       sel(x)       = (1 - w(x)) * asc(x) + w(x) * desc(x)
 *       right_floor  = 1 / (1 + exp(-(sel_inf(1) + sel_inf_dev(1))))   [logit]
 *     Parameters (reuses existing arrays; all four support annual deviates):
 *       sel_inf(0)     = peak bin (mode);            TV deviate: sel_inf_dev(0)
 *       log_sel_slp(0) = log(sigma_ascending);       TV deviate: log_sel_slp_dev(0)
 *       log_sel_slp(1) = log(sigma_descending);      TV deviate: log_sel_slp_dev(1)
 *       sel_inf(1)     = logit(right_floor) (~SS3 P6/end_logit). logit -> -Inf gives a
 *                        fully dome-shaped curve; logit -> +Inf collapses to logistic
 *                        (ascending only). TV deviate: sel_inf_dev(1)
 * * @param nspp Number of species.
 * @param n_flt Number of fleets (fisheries and surveys).
 * @param nyrs Total years (including hindcast and projection).
 * @param nyrs_hind Number of hindcast years.
 * @param styr Start year (for indexing empirical data).
 * @param nsex Vector of number of sexes per species.
 * @param nages Vector of max age bins per species.
 * @param nlengths Vector of max length bins per species.
 * @param lengths Matrix of length bin boundaries (nspp x max_length_bins).
 * @param flt_spp Vector mapping fleet index to species index.
 * @param flt_sel_type Vector of selectivity form.
 * @param flt_sel_dim Age or length based selectivity
 * @param bin_first_selected Vector of minimum age/length bins selected per fleet.
 * @param flt_n_sel_bins Vector of maximum bins with estimated coefficients per fleet.
 * @param sel_norm_bin1 Vector of bins used for normalization (or control flags) per fleet.
 * @param sel_norm_bin2 Vector of upper bins for mean normalization per fleet.
 * @param emp_sel_obs Matrix of empirical selectivity observation values.
 * @param emp_sel_ctl Control matrix for empirical selectivity indexing.
 * @param log_sel_slp Matrix of log-scale logistic slope parameters.
 * @param log_sel_slp_dev 4D array of yearly slope deviations.
 * @param sel_inf Matrix of bin-at-inflection parameters.
 * @param sel_inf_dev 4D array of yearly inflection deviations.
 * @param sel_coff 3D array of bin-specific coefficients for non-parametric models.
 * @param sel_coff_dev 4D array of yearly coefficient deviations.
 * @param avg_sel [Output] 3D array of average selectivity per fleet/sex/year (modified by reference).
 * @param non_par_sel [Output] 4D array of unnormalized non-parametric selectivity (modified by reference).
 * @param log_non_par_sel [Output] The same curve on the log scale, which is the scale the
 *   non-parametric selectivity penalties are written on (modified by reference).
 * @param sel_at_length [Output] 4D array of final length-based selectivity (modified by reference).
 * @param sel_at_age [Output] 4D array of final age-based selectivity (modified by reference).
 * @param growth_matrix 5D array defining the length-at-age probability transition matrix.
 */
template<class Type>
void calculate_selectivity(
    const int& nspp,
    const int& n_flt,
    const int& nyrs,
    const int& nyrs_hind,
    const int& styr,
    const vector<int>&  nsex,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    matrix<Type>& lengths,
    const vector<int>&  flt_spp,
    const vector<int>&  flt_sel_type,
    const vector<int>&  flt_varying_sel,
    const vector<int>&  flt_sel_dim,
    const vector<int>&  bin_first_selected,
    const vector<int>&  flt_n_sel_bins,
    const vector<int>&  flt_sel_cap_bin,
    const vector<int>&  sel_norm_bin1,
    const vector<int>&  sel_norm_bin2,
    const vector<int>&  sel_norm_scope,
    const vector<int>&  flt_sel_start_yr,
    matrix<Type>& emp_sel_obs,
    matrix<int>& emp_sel_ctl,
    array<Type>& log_sel_slp,
    array<Type>& log_sel_slp_dev,
    array<Type>& sel_inf,
    array<Type>& sel_inf_dev,
    array<Type>& sel_coff,
    array<Type>& sel_coff_dev,
    array<Type> &avg_sel,
    array<Type> &non_par_sel,
    array<Type> &log_non_par_sel,
    array<Type> &sel_at_length,
    array<Type> &sel_at_age,
    array<Type>& growth_matrix,
    // Environmental-linkage offsets, one log-scale and one natural-scale tensor
    // per parameter family. All zero unless a selectivity linkage was supplied,
    // so a model without one is numerically unchanged.
    // Non-const references only because TMB's multi-index array operator() is
    // not const-qualified; these tensors are read, never written, here.
    array<Type>& sel_slp_off,      array<Type>& sel_slp_off_nat,
    array<Type>& sel_inf_off,      array<Type>& sel_inf_off_nat,
    array<Type>& sel_coff_off,     array<Type>& sel_coff_off_nat,
    // Per-sex apical height: exp(log_sel_apical + log offset) + natural offset
    // multiplies one sex's whole curve. Compiled in only when a selectivity
    // linkage names `apical` (sel_apical_on), so every other model keeps its
    // AD tape unchanged.
    const int& sel_apical_on,
    array<Type>& log_sel_apical,   // [n_flt, max_sex]
    array<Type>& sel_apical_off,   array<Type>& sel_apical_off_nat   // [n_flt, max_sex, nyrs]
) {
  sel_at_age.setZero();
  sel_at_length.setZero();
  Type max_sel = 0;
  Type avgsel_tmp = 0;

  // --- 1. EMPIRICAL SELECTIVITY (Assumed Age-based from input) ---
  for (int i = 0; i < emp_sel_obs.rows(); i++) {
    int flt = emp_sel_ctl(i, 0) - 1;
    int sp  = emp_sel_ctl(i, 1) - 1;
    int s_idx = emp_sel_ctl(i, 2);
    int yr = emp_sel_ctl(i, 3) - styr;

    if (flt_sel_type(flt) == 0 && yr < nyrs_hind) {
      int s_start = (s_idx == 0) ? 0 : s_idx - 1;
      int s_end   = (s_idx == 0 && nsex(sp) == 2) ? 1 : s_start;

      for (int sex = s_start; sex <= s_end; sex++) {
        for (int age = 0; age < nages(sp); age++) {
          if (!isNA(emp_sel_obs(i, age))) {
            sel_at_age(flt, sex, age, yr) = emp_sel_obs(i, age);
          }
        }
      }
    }
  }

  // --- 2. ESTIMATED SELECTIVITY (Unified Age & Length) ---
  for (int flt = 0; flt < n_flt; flt++) {
    int sp = flt_spp(flt);
    int sel_type = flt_sel_type(flt);
    if (sel_type == 0) continue;

    // NonParametricIntegrable (13) builds its curve two ways: a random walk on the
    // base coefficients when its deviations are a walk, otherwise the shared
    // non-parametric branch. Time_varying_sel picks which, so the switch below
    // dispatches on this rather than on the Selectivity code alone.
    // SEL_CASE_NP_INTEGRABLE_WALK is internal; it is not a Selectivity value.
    const int SEL_CASE_NP_INTEGRABLE_WALK = -13;
    int sel_case = (sel_type == 13 && flt_varying_sel(flt) == 4)
                     ? SEL_CASE_NP_INTEGRABLE_WALK : sel_type;

    bool is_length_based = flt_sel_dim(flt) == 1;
    int nbins =  is_length_based? nlengths(sp) : nages(sp);
    int n_sel_bins = flt_n_sel_bins(flt);
    Type binwidth = is_length_based ? (lengths(sp, 1) - lengths(sp, 0)) : Type(1.0);

    // Uncapped, per-year-centered log-selectivity, carried across years for the
    // NonParametricPM (type 9) random walk (the realized curve is then capped).
    array<Type> np_unc(nsex(sp), nbins, nyrs_hind); np_unc.setZero();

    for (int yr = 0; yr < nyrs_hind; yr++) {
      for (int sex = 0; sex < nsex(sp); sex++) {

        switch (sel_case) {
        case 1: // Logistic
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? (lengths(sp, bin) + 0.5 * binwidth) : Type(bin + 1);
            // slope stored on the log scale: log-link offset rides inside the
            // exp (multiplicative on the natural slope), natural-link offset
            // adds after. inflection stored natural: identity-link offset adds,
            // log-link offset multiplies. Both no-op when the offsets are zero.
            Type slope = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr)
                             + sel_slp_off(0, flt, sex, yr)) + sel_slp_off_nat(0, flt, sex, yr);
            Type inf   = (sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr)
                          + sel_inf_off_nat(0, flt, sex, yr)) * exp(sel_inf_off(0, flt, sex, yr));
            Type val = 1.0 / (1.0 + exp(-slope * (x_val - inf)));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        // NonParametricIntegrable (13) under Time_varying_sel = "RandomWalk": the
        // Ianelli curve on the base coefficients plus the running sum of increments
        // (years after the fleet's start year), centred for output only, so with no
        // increments every year is NonParametric's curve. No cap; the increments are
        // scored in ceattle.cpp. Under "Off" or "IID" the form falls to the shared
        // non-parametric branch below.
        case SEL_CASE_NP_INTEGRABLE_WALK: {
          for(int bin = 0; bin < n_sel_bins; bin++){
            Type inc = (yr > flt_sel_start_yr(flt)) ? sel_coff_dev(flt, sex, bin, yr) : Type(0.0);
            np_unc(sex, bin, yr) = (yr > 0 ? np_unc(sex, bin, yr - 1) : Type(0.0)) + inc;
            non_par_sel(flt, sex, bin, yr) = sel_coff(flt, sex, bin) + np_unc(sex, bin, yr);
          }
          { vector<Type> base(n_sel_bins);
            for(int bin = 0; bin < n_sel_bins; bin++) base(bin) = non_par_sel(flt, sex, bin, yr);
            avg_sel(flt, sex, yr) = log_mean_exp(base); }
          for(int bin = n_sel_bins; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) = non_par_sel(flt, sex, n_sel_bins - 1, yr);
          }
          { vector<Type> ls(nbins);
            for(int bin = 0; bin < nbins; bin++) ls(bin) = non_par_sel(flt, sex, bin, yr);
            avgsel_tmp = log_mean_exp(ls); }
          for(int bin = 0; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) -= avgsel_tmp;
            log_non_par_sel(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            non_par_sel(flt, sex, bin, yr) = exp(non_par_sel(flt, sex, bin, yr));
            if (is_length_based) sel_at_length(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            else                 sel_at_age(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
          }
          break;
        }

        case 9: { // NonParametricPM (RTMB "rpm"): random walk on the per-year-
                  // renormalized log-selectivity, then a flat age-cap.
          // sel_coff = base coffs (year styr, ages 0..n_sel_bins-1); sel_coff_dev =
          // RAW per-year increments placed at the year they apply. Ages >= n_sel_bins
          // plateau at the last coff. The UNCAPPED centered curve (np_unc) is carried
          // forward for the walk; the realized curve is that capped flat at
          // flt_sel_cap_bin and re-centered (mean(exp)=1).
          // For years at or before the fleet's selectivity start year, the curve is
          // built directly from the base coefficients (a fresh mean-centering of that
          // single vector each year) rather than by carrying the running random walk.
          // This follows the AMAK convention of setting a survey's base selectivity
          // once at its start year; beyond the start year the curve carries forward as
          // a random walk. For a fleet that starts at styr (start_yr = 0) this is the
          // ordinary base-year build.
          bool from_base = (yr <= flt_sel_start_yr(flt));
          for(int bin = 0; bin < nbins; bin++){
            Type prev = from_base ? sel_coff(flt, sex, (bin < n_sel_bins ? bin : n_sel_bins - 1))
                                  : np_unc(sex, (bin < n_sel_bins ? bin : n_sel_bins - 1), yr - 1);
            Type inc  = (bin < n_sel_bins) ? sel_coff_dev(flt, sex, bin, yr) : Type(0.0);
            np_unc(sex, bin, yr) = prev + inc;
          }
          for(int bin = n_sel_bins; bin < nbins; bin++) np_unc(sex, bin, yr) = np_unc(sex, n_sel_bins - 1, yr);
          // Bins below bin_first_selected (e.g. the acoustic-survey age-1) are held at
          // 0 before each year's mean-centering, following the AMAK convention that a
          // fixed first-bin log-selectivity is 0 prior to normalization. This keeps the
          // excluded bin out of the shared normalization of the selected ages and out
          // of the deviate/curvature penalties.
          for(int bin = 0; bin < bin_first_selected(flt); bin++) np_unc(sex, bin, yr) = 0.0;
          { vector<Type> u(nbins);
            for(int bin = 0; bin < nbins; bin++) u(bin) = np_unc(sex, bin, yr);
            Type m = log_mean_exp(u);
            for(int bin = 0; bin < nbins; bin++) np_unc(sex, bin, yr) -= m; }
          { vector<Type> cl(nbins);
            for(int bin = 0; bin < nbins; bin++) cl(bin) = np_unc(sex, bin, yr);
            int cap = flt_sel_cap_bin(flt);
            if(cap >= 0) for(int bin = cap + 1; bin < nbins; bin++) cl(bin) = cl(cap);
            Type m2 = log_mean_exp(cl);
            for(int bin = 0; bin < nbins; bin++){
              log_non_par_sel(flt, sex, bin, yr) = cl(bin) - m2;
              non_par_sel(flt, sex, bin, yr) = exp(cl(bin) - m2);
              if (is_length_based) sel_at_length(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
              else                 sel_at_age(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            }
          }
          break;
        }

        case 13: // NonParametricIntegrable under "Off"/"IID"; penalties on the base (ceattle.cpp)
        case 2: // Non-parametric (Ianelli style)
          for(int bin = 0; bin < n_sel_bins; bin++) {
            non_par_sel(flt, sex, bin, yr) = sel_coff(flt, sex, bin) + sel_coff_dev(flt, sex, bin, yr);
          }
          // AMAK avgsel: the level of the estimated coefficients, which the
          // mean-centering below removes from the curve. The penalty on it is
          // the only thing holding that level (ceattle.cpp, JNLL_SEL_NONPARAM).
          { vector<Type> base(n_sel_bins);
            for(int bin = 0; bin < n_sel_bins; bin++) base(bin) = non_par_sel(flt, sex, bin, yr);
            avg_sel(flt, sex, yr) = log_mean_exp(base); }

          for(int bin = n_sel_bins; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) = non_par_sel(flt, sex, n_sel_bins - 1, yr);
          }

          { vector<Type> ls(nbins);
            for(int bin = 0; bin < nbins; bin++) ls(bin) = non_par_sel(flt, sex, bin, yr);
            avgsel_tmp = log_mean_exp(ls); }

          for(int bin = 0; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) -= avgsel_tmp;
            log_non_par_sel(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            non_par_sel(flt, sex, bin, yr) = exp(non_par_sel(flt, sex, bin, yr));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            else sel_at_age(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
          }
          break;

        case 3: // Double Logistic
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? (lengths(sp, bin) + 0.5 * binwidth) : Type(bin + 1);
            Type slp1 = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr)
                            + sel_slp_off(0, flt, sex, yr)) + sel_slp_off_nat(0, flt, sex, yr);
            Type inf1 = (sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr)
                         + sel_inf_off_nat(0, flt, sex, yr)) * exp(sel_inf_off(0, flt, sex, yr));
            Type slp2 = exp(log_sel_slp(1, flt, sex) + log_sel_slp_dev(1, flt, sex, yr)
                            + sel_slp_off(1, flt, sex, yr)) + sel_slp_off_nat(1, flt, sex, yr);
            Type inf2 = (sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr)
                         + sel_inf_off_nat(1, flt, sex, yr)) * exp(sel_inf_off(1, flt, sex, yr));
            Type val = (1.0 / (1.0 + exp(-slp1 * (x_val - inf1)))) * (1.0 - (1.0 / (1.0 + exp(-slp2 * (x_val - inf2)))));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        case 4: // Descending Logistic
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? (lengths(sp, bin) + 0.5 * binwidth) : Type(bin + 1);
            Type slp2 = exp(log_sel_slp(1, flt, sex) + log_sel_slp_dev(1, flt, sex, yr)
                            + sel_slp_off(1, flt, sex, yr)) + sel_slp_off_nat(1, flt, sex, yr);
            Type inf2 = (sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr)
                         + sel_inf_off_nat(1, flt, sex, yr)) * exp(sel_inf_off(1, flt, sex, yr));
            Type val = (1.0 - (1.0 / (1.0 + exp(-slp2 * (x_val - inf2)))));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        case 5: // Hake/Taylor Non-parametric
          {
            Type cum_sum = 0;
            for (int bin = bin_first_selected(flt); bin < nbins; bin++) {
              if (bin < n_sel_bins) {
                cum_sum += sel_coff(flt, sex, bin) + sel_coff_dev(flt, sex, bin, yr);
              }
              if (is_length_based) sel_at_length(flt, sex, bin, yr) = cum_sum;
              else sel_at_age(flt, sex, bin, yr) = cum_sum;
            }

            // Normalized per sex here, not in normalize_and_project_selectivity(),
            // so both sexes reach 1: Sel_norm_scope is inert on Hake, which cannot
            // carry a sex difference in level. inst/dev/TODO-selectivity.md
            max_sel = -1e10;
            if (sel_norm_bin1(flt) >= 0 && sel_norm_bin2(flt) < 0) {
              max_sel = is_length_based ? sel_at_length(flt, sex, sel_norm_bin1(flt), yr) : sel_at_age(flt, sex, sel_norm_bin1(flt), yr);
            } else {
              for(int bin=0; bin<nbins; bin++) {
                Type val = is_length_based ? sel_at_length(flt, sex, bin, yr) : sel_at_age(flt, sex, bin, yr);
                max_sel = max2(max_sel, val);   // fold, don't branch (see above)
              }
            }
            for (int bin = 0; bin < nbins; bin++) {
              if (is_length_based) sel_at_length(flt, sex, bin, yr) = exp(sel_at_length(flt, sex, bin, yr) - max_sel);
              else sel_at_age(flt, sex, bin, yr) = exp(sel_at_age(flt, sex, bin, yr) - max_sel);
            }
          }
          break;

        case 8: { // Double Normal (4-param, with right-tail floor; see Doxygen header above)
          // Parameters:
          //   sel_inf(0)     = peak bin (mode of selectivity)
          //   log_sel_slp(0) = log(sigma_ascending)   - ascending limb SD
          //   log_sel_slp(1) = log(sigma_descending)  - descending limb SD
          //   sel_inf(1)     = logit(right_floor) — right-tail floor, analogous to SS3 P6 (end_logit).
          //                    right_floor -> 0: fully dome-shaped; right_floor -> 1: logistic (ascending only).
          // DoubleNormal reuses the logistic slots: sel_inf(0)=peak,
          // log_sel_slp(0/1)=log sigma asc/desc, sel_inf(1)=logit right-floor.
          // The linkage offsets ride in the same position as the deviates, so
          // slp/inf params act on sigma/peak/floor here without special-casing.
          Type peak        = (sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr)
                              + sel_inf_off_nat(0, flt, sex, yr)) * exp(sel_inf_off(0, flt, sex, yr));
          Type sigma_asc   = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr)
                              + sel_slp_off(0, flt, sex, yr)) + sel_slp_off_nat(0, flt, sex, yr);
          Type sigma_desc  = exp(log_sel_slp(1, flt, sex) + log_sel_slp_dev(1, flt, sex, yr)
                              + sel_slp_off(1, flt, sex, yr)) + sel_slp_off_nat(1, flt, sex, yr);
          // right_floor is logit-parameterised (sel_inf(1) = logit(right_floor)), so the
          // linkage offsets act on the logit scale -- additive (nat) and multiplicative
          // (exp) inside the logistic, matching `peak` and the other inflection slots. This
          // keeps right_floor in [0, 1] for any offset (applying exp() to the probability
          // could push it > 1). Bit-identical at zero offset.
          Type right_floor = 1.0 / (1.0 + exp(-((sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr)
                              + sel_inf_off_nat(1, flt, sex, yr)) * exp(sel_inf_off(1, flt, sex, yr)))));
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val      = is_length_based ? (lengths(sp, bin) + 0.5 * binwidth) : Type(bin + 1);
            // Smooth logistic blend: ~0 left of peak, ~1 right of peak.
            Type w          = 1.0 / (1.0 + exp(-20.0 * (x_val - peak)));
            Type asc_gauss  = exp(-0.5 * pow((x_val - peak) / sigma_asc,  2.0));
            Type desc_gauss = exp(-0.5 * pow((x_val - peak) / sigma_desc, 2.0));
            // Right tail approaches right_floor (mirrors SS3 P6 / end_logit behaviour)
            Type val_desc   = right_floor + (1.0 - right_floor) * desc_gauss;
            Type val        = (1.0 - w) * asc_gauss + w * val_desc;
            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else                 sel_at_age(flt, sex, bin, yr) = val;
          }
          break;
        }

        case 6: // 2D AR1-age x year
        case 7: // 3D AR1 conditional variance
          for (int bin = 0; bin < nbins; bin++) {

            // Safely cap the bin index to prevent out-of-bounds access
            int active_bin = (bin < n_sel_bins) ? bin : (n_sel_bins - 1);
            Type val = 1.0 / (1.0 + exp(-(sel_coff(flt, sex, active_bin) + sel_coff_dev(flt, sex, active_bin, yr))));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        case 11: { // LogisticPM (ADMB AMAK "pm" bottom-trawl survey form)
          // Standard 2-parameter logistic over all bins, but with MULTIPLICATIVE
          // time-varying deviations on slope and inflection (matching AMAK's
          //   sel = 1/(1 + exp(-slp*exp(slp_dev) * (age - a50*exp(a50_dev)))) ),
          // plus a FREE first-bin (age-1) log-selectivity independent of the
          // logistic (AMAK sel_age_one_bts*exp(sel_age_one_bts_dev)). The first-bin
          // base and its deviates are stored in the unused descending-limb slots
          // sel_inf(1)/sel_inf_dev(1). No internal normalization (set Sel_norm_bin1
          // = NA): AMAK does not renormalize the BTS curve and age-1 may exceed 1.
          // NOTE: AMAK evaluates the logistic at age_vector(j) = j + 0.5 (mid-age),
          // so the age-based x is (bin + 1) + 0.5 = bin + 1.5, NOT bin + 1 as in
          // the standard Logistic (case 1). This 0.5 shift cannot be folded into a50
          // because the inflection deviate is multiplicative (a50*exp(dev)).
          // LogisticPM deviates are MULTIPLICATIVE (a50 * exp(dev)), so a log-link
          // linkage offset rides inside the same exp and an identity offset adds
          // to the base -- both leaving the log-link = multiplicative-on-natural
          // meaning that holds for the other forms.
          Type slope = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr)
                           + sel_slp_off(0, flt, sex, yr)) + sel_slp_off_nat(0, flt, sex, yr);
          Type inf   = (sel_inf(0, flt, sex) + sel_inf_off_nat(0, flt, sex, yr))
                         * exp(sel_inf_dev(0, flt, sex, yr) + sel_inf_off(0, flt, sex, yr));
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? (lengths(sp, bin) + 0.5 * binwidth) : Type(bin + 1.5);
            Type val = 1.0 / (1.0 + exp(-slope * (x_val - inf)));
            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else                 sel_at_age(flt, sex, bin, yr) = val;
          }
          // Free first-bin (age-1) log-selectivity override
          Type log_s1 = (sel_inf(1, flt, sex) + sel_inf_off_nat(1, flt, sex, yr))
                          * exp(sel_inf_dev(1, flt, sex, yr) + sel_inf_off(1, flt, sex, yr));
          if (is_length_based) sel_at_length(flt, sex, 0, yr) = exp(log_s1);
          else                 sel_at_age(flt, sex, 0, yr) = exp(log_s1);
          break;
        }
        }

      } // End sex
    } // End year

    // --- 2b. PER-SEX APICAL HEIGHT ---
    // One sex's whole curve times exp(log_sel_apical + log-link offset), plus
    // the identity-link offset, the convention the log-scale slopes use. Applied
    // after every form's own centring or normalization (NonParametric,
    // NonParametricPM and Hake rescale inside the switch) and before the shared
    // normalizer, which must pool its reference across sexes (AcrossSexes) or be
    // off for the offset to survive; .check_sel_apical_rows() refuses the rest.
    if (sel_apical_on == 1) {
      for (int yr = 0; yr < nyrs_hind; yr++) {
        for (int sex = 0; sex < nsex(sp); sex++) {
          Type scale = exp(log_sel_apical(flt, sex) + sel_apical_off(flt, sex, yr))
                       + sel_apical_off_nat(flt, sex, yr);
          for (int bin = 0; bin < nbins; bin++) {
            if (is_length_based) sel_at_length(flt, sex, bin, yr) *= scale;
            else                 sel_at_age(flt, sex, bin, yr) *= scale;
          }
        }
      }
    }

    // --- 3. NORMALIZATION & PROJECTION ---
    normalize_and_project_selectivity(
      flt, nyrs_hind, nyrs, flt_spp, flt_sel_type, flt_sel_dim, bin_first_selected, nages, nlengths, nsex,
      sel_norm_bin1, sel_norm_bin2, sel_norm_scope,
      is_length_based ? sel_at_length : sel_at_age
    );

    // --- 4. CONVERT LENGTH TO AGE (If necessary) ---
    if (is_length_based) {
      for (int yr = 0; yr < nyrs; yr++) {
        for (int sex = 0; sex < nsex(sp); sex++) {
          int wtind = nspp * 2 + flt;
          convert_length_selectivity(
            flt, sp, sex, yr, wtind, nages, nlengths,
            sel_at_length, growth_matrix, sel_at_age
          );
        }
      }
    }
  } // End fleet
}

#endif
