#ifndef SELECTIVITY_HPP
#define SELECTIVITY_HPP

/**
 * @brief Zeroes, normalizes and projects ONE fleet's selectivity (age or length based).
 *
 * @details Normalization never applies to Hake (5), LogisticPM (11) or DoubleNormalSS3
 * (15). Under Sel_norm_scope = AcrossSexes the less-selected sex peaks below 1.
 *
 * @param sel_norm_bin1 Normalization bin, 0-based; -99 = the maximum, -999 = off.
 * @param sel_norm_bin2 Upper bin for the mean over a range; negative = no range.
 * @param sel_norm_scope 0 = WithinSex (a reference per sex), 1 = AcrossSexes (pooled).
 * @param selectivity [fleet, sex, age or length bin, year], modified in place.
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
    const vector<int>&  nlengths_pop,
    const vector<int>&  nsex,
    const vector<int>&  sel_norm_bin1,
    const vector<int>&  sel_norm_bin2,
    const vector<int>&  sel_norm_scope,
    array<Type> &selectivity
) {
  int sp = flt_spp(flt);
  int sel_type = flt_sel_type(flt);
  int nbins = (flt_sel_dim(flt) == 0) ? nages(sp) : nlengths_pop(sp);

  // 1. Bins below the first selected bin are not selected
  for (int yr = 0; yr < nyrs_hind; yr++) {
    for (int sex = 0; sex < nsex(sp); sex++) {
      for (int bin = 0; bin < nbins && bin < bin_first_selected(flt); bin++) {
        selectivity(flt, sex, bin, yr) = 0;
      }
    }
  }

  // 2. Normalization. Hake scales itself per sex, LogisticPM reads Sel_norm_bin as a
  // penalty range, and DoubleNormalSS3 peaks at 1. Codes 10, 12 and 14 are unused.
  bool normalize = (sel_norm_bin1(flt) > -500) && (sel_type != 5) && (sel_type != 12)
                   && (sel_type != 11) && (sel_type != 15);
  if (normalize) {
    bool at_bin     = (sel_norm_bin1(flt) >= 0);
    bool over_range = at_bin && (sel_norm_bin2(flt) >= 0);
    bool across     = (sel_norm_scope(flt) == 1);
    vector<Type> ref(nsex(sp));

    for (int yr = 0; yr < nyrs_hind; yr++) {
      for (int sex = 0; sex < nsex(sp); sex++) {
        if (over_range) {
          int n_range = sel_norm_bin2(flt) - sel_norm_bin1(flt) + 1;
          ref(sex) = 0;
          for (int bin = sel_norm_bin1(flt); bin <= sel_norm_bin2(flt); bin++) {
            ref(sex) += selectivity(flt, sex, bin, yr) / n_range;
          }
        } else if (at_bin) {
          // One bin, floored so the division stays finite
          ref(sex) = max2(Type(0.001), selectivity(flt, sex, sel_norm_bin1(flt), yr));
        } else {
          // max2(), not an `if`: TMB records the model once, so an `if` would freeze
          // which bin is the maximum.
          ref(sex) = 0;
          for (int bin = 0; bin < nbins; bin++) {
            ref(sex) = max2(ref(sex), selectivity(flt, sex, bin, yr));
          }
        }
      }

      if (across) {
        Type pooled = ref(0);
        for (int sex = 1; sex < nsex(sp); sex++) pooled = max2(pooled, ref(sex));
        for (int sex = 0; sex < nsex(sp); sex++) ref(sex) = pooled;
      }

      for (int bin = 0; bin < nbins; bin++) {
        for (int sex = 0; sex < nsex(sp); sex++) {
          selectivity(flt, sex, bin, yr) /= ref(sex);
        }
      }
    }
  }

  // 3. Projection years repeat the last hindcast year
  for (int yr = nyrs_hind; yr < nyrs; yr++) {
    for (int sex = 0; sex < nsex(sp); sex++) {
      for (int bin = 0; bin < nbins; bin++) {
        selectivity(flt, sex, bin, yr) = selectivity(flt, sex, bin, nyrs_hind - 1);
      }
    }
  }
}


/**
 * @brief Converts one fleet/sex/year of length-based selectivity to ages: the sum over
 * population length bins of selectivity-at-length x Prob(length | age).
 */
template <class Type>
void convert_length_selectivity(
    const int& flt,
    const int& sp,
    const int& sex,
    const int& yr,
    const int& wtind,
    const vector<int>& nages,
    const vector<int>& nlengths_pop,
    array<Type>& sel_length,
    array<Type>& growth_matrix,
    array<Type>& sel_at_age
) {
  for (int age = 0; age < nages(sp); age++) {
    sel_at_age(flt, sex, age, yr) = 0.0;
    for (int lp = 0; lp < nlengths_pop(sp); lp++) {
      sel_at_age(flt, sex, age, yr) += growth_matrix(wtind, sex, age, lp, yr) * sel_length(flt, sex, lp, yr);
    }
  }
}


/**
 * @brief Calculates selectivity for every fleet, sex, bin (age or length) and year.
 *
 * @details An age-based curve is evaluated at the bin number bin + 1, without minage,
 * so an inflection or peak is on the bin scale; a length-based one at bin midpoints (cm).
 *
 * Forms (fleet_control$Selectivity; R/0-switches.R is the authority):
 * - 0 Fixed: read from emp_sel_obs; age-based.
 * - 1 Logistic.
 * - 2 NonParametric (Ianelli et al. 2018, AMAK): one log coefficient per bin, mean-centred.
 * - 3 DoubleLogistic (Dorn and Methot 1990).
 * - 4 DescendingLogistic.
 * - 5 Hake (Taylor et al. 2014): exp of the cumulative sum of coefficients.
 * - 6 2DAR1: inverse logit of coefficients with AR1 deviates by bin and year.
 * - 7 3DAR1 (Cheng et al. 2024): as 6, with a cohort term.
 * - 8 DoubleNormal: Gaussian dome with a right-tail floor.
 * - 9 NonParametricPM (AMAK "pm"): random walk on the log curve, capped at Sel_cap_bin.
 * - 11 LogisticPM (AMAK "pm" bottom-trawl survey): multiplicative deviates, free first bin.
 * - 13 NonParametricIntegrable: form 2 with deviates random_sel can integrate.
 * - 15 DoubleNormalSS3: Stock Synthesis size pattern 24 (age pattern 20).
 *
 * Which linkage offsets a form accepts is checked in R, by .check_sel_linkage_support().
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
    const vector<int>&  nlengths_pop,
    matrix<Type>& lengths_pop,
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
    // Linkage offsets: log link, and identity link (_nat). The coefficient ones are unread.
    array<Type>& sel_slp_off,      array<Type>& sel_slp_off_nat,
    array<Type>& sel_inf_off,      array<Type>& sel_inf_off_nat,
    array<Type>& sel_coff_off,     array<Type>& sel_coff_off_nat,
    const int& sel_apical_on,
    array<Type>& log_sel_apical,   // [n_flt, max_sex]
    array<Type>& sel_apical_off,   array<Type>& sel_apical_off_nat,  // [n_flt, max_sex, nyrs]
    array<Type>& sel_dn6,          // [6, n_flt, max_sex]
    array<Type>& sel_dn6_off,      array<Type>& sel_dn6_off_nat,     // [6, n_flt, max_sex, nyrs]
    matrix<int>& sel_dn6_ends      // [n_flt, 2]: SS3 initial/final selectivity used (1) or -999 (0)
) {
  sel_at_age.setZero();
  sel_at_length.setZero();

  // --- 1. FIXED SELECTIVITY (code 0) ---
  for (int i = 0; i < emp_sel_obs.rows(); i++) {
    int flt    = emp_sel_ctl(i, 0) - 1;
    int sp     = emp_sel_ctl(i, 1) - 1;
    int sex_in = emp_sel_ctl(i, 2);       // 0 = every sex
    int yr     = emp_sel_ctl(i, 3) - styr;
    if (flt_sel_type(flt) != 0 || yr >= nyrs_hind) continue;

    int sex_first = (sex_in == 0) ? 0 : sex_in - 1;
    int sex_last  = (sex_in == 0 && nsex(sp) == 2) ? 1 : sex_first;
    for (int sex = sex_first; sex <= sex_last; sex++) {
      for (int age = 0; age < nages(sp); age++) {
        if (!isNA(emp_sel_obs(i, age))) sel_at_age(flt, sex, age, yr) = emp_sel_obs(i, age);
      }
    }
  }

  // --- 2. ESTIMATED SELECTIVITY ---
  for (int flt = 0; flt < n_flt; flt++) {
    int sp = flt_spp(flt);
    int sel_type = flt_sel_type(flt);
    if (sel_type == 0) continue;

    // Form 13 is built differently under Time_varying_sel = "RandomWalk" (4). The
    // negative label is internal, not a Selectivity code.
    const int SEL_CASE_NP_INTEGRABLE_WALK = -13;
    int sel_case = (sel_type == 13 && flt_varying_sel(flt) == 4)
                     ? SEL_CASE_NP_INTEGRABLE_WALK : sel_type;

    // A length-based curve runs over the POPULATION length bins. `sel` is that array
    // itself, not a copy.
    bool is_length_based = flt_sel_dim(flt) == 1;
    array<Type>& sel = is_length_based ? sel_at_length : sel_at_age;
    int nbins = is_length_based ? nlengths_pop(sp) : nages(sp);
    int n_sel_bins = flt_n_sel_bins(flt);

    // x-axis: bin + 1 by age; by length, each bin's midpoint from its own width
    // (SS3's len_bins_m), since SS3 grids often widen the tail bins.
    vector<Type> x(nbins);
    for (int b = 0; b < nbins; b++) {
      if (!is_length_based) {
        x(b) = Type(b + 1);
      } else if (nbins < 2) {
        x(b) = lengths_pop(sp, b);
      } else {
        Type w = (b + 1 < nbins) ? (lengths_pop(sp, b + 1) - lengths_pop(sp, b))
                                 : (lengths_pop(sp, b) - lengths_pop(sp, b - 1));
        x(b) = lengths_pop(sp, b) + 0.5 * w;
      }
    }
    // SS3's binwidth2 (SS_readdata_330.tpl:1644), read only by DoubleNormalSS3's peak2
    // (SS_selex.tpl:153).
    Type binwidth = (is_length_based && nbins > 1)
      ? (lengths_pop(sp, nbins / 2) - lengths_pop(sp, nbins / 2 - 1)) : Type(1.0);

    // Random-walk state for forms 9 and 13, carried across years.
    array<Type> log_walk(nsex(sp), nbins, nyrs_hind); log_walk.setZero();

    for (int yr = 0; yr < nyrs_hind; yr++) {
      for (int sex = 0; sex < nsex(sp); sex++) {

        // Local functions: slope and inflection of logistic slot k (0 = ascending, 1 =
        // descending), base + deviate + linkage offsets (log link multiplies, identity adds).
        auto slope = [&](int k) -> Type {
          return exp(log_sel_slp(k, flt, sex) + log_sel_slp_dev(k, flt, sex, yr)
                     + sel_slp_off(k, flt, sex, yr)) + sel_slp_off_nat(k, flt, sex, yr);
        };
        auto inflection = [&](int k) -> Type {
          return (sel_inf(k, flt, sex) + sel_inf_dev(k, flt, sex, yr)
                  + sel_inf_off_nat(k, flt, sex, yr)) * exp(sel_inf_off(k, flt, sex, yr));
        };

        switch (sel_case) {

        case 1: { // Logistic
          Type slp = slope(0);
          Type inf = inflection(0);
          for (int bin = 0; bin < nbins; bin++) {
            sel(flt, sex, bin, yr) = invlogit(slp * (x(bin) - inf));
          }
          break;
        }

        case 2:  // NonParametric
        case 13: // NonParametricIntegrable, "Off" or "IID"
        case SEL_CASE_NP_INTEGRABLE_WALK: { // NonParametricIntegrable, "RandomWalk"
          // Log coefficients: base + this year's deviate, or + the running sum of
          // increments after the start year under the random walk.
          vector<Type> coef(n_sel_bins);
          for (int bin = 0; bin < n_sel_bins; bin++) {
            if (sel_case == SEL_CASE_NP_INTEGRABLE_WALK) {
              Type inc = (yr > flt_sel_start_yr(flt)) ? sel_coff_dev(flt, sex, bin, yr) : Type(0.0);
              log_walk(sex, bin, yr) = (yr > 0 ? log_walk(sex, bin, yr - 1) : Type(0.0)) + inc;
              coef(bin) = sel_coff(flt, sex, bin) + log_walk(sex, bin, yr);
            } else {
              coef(bin) = sel_coff(flt, sex, bin) + sel_coff_dev(flt, sex, bin, yr);
            }
          }
          // AMAK avgsel: the log mean that centring removes. Form 2 penalizes it
          // (ceattle.cpp, JNLL_SEL_NONPARAM).
          avg_sel(flt, sex, yr) = log_mean_exp(coef);

          // Plateau past N_sel_bins, then centre so mean selectivity is 1.
          vector<Type> log_curve(nbins);
          for (int bin = 0; bin < nbins; bin++) {
            int b = (bin < n_sel_bins) ? bin : (n_sel_bins - 1);
            log_curve(bin) = coef(b);
          }
          Type centre = log_mean_exp(log_curve);
          for (int bin = 0; bin < nbins; bin++) {
            log_non_par_sel(flt, sex, bin, yr) = log_curve(bin) - centre;
            non_par_sel(flt, sex, bin, yr) = exp(log_non_par_sel(flt, sex, bin, yr));
            sel(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
          }
          break;
        }

        case 3: { // DoubleLogistic
          // Descending limb as a negative-slope logistic, not 1 minus a logistic, so the
          // oldest ages keep a positive, accurate value on a steep slope.
          Type slp1 = slope(0);
          Type inf1 = inflection(0);
          Type slp2 = slope(1);
          Type inf2 = inflection(1);
          for (int bin = 0; bin < nbins; bin++) {
            sel(flt, sex, bin, yr) = invlogit(slp1 * (x(bin) - inf1)) * invlogit(-slp2 * (x(bin) - inf2));
          }
          break;
        }

        case 4: { // DescendingLogistic
          Type slp2 = slope(1);
          Type inf2 = inflection(1);
          for (int bin = 0; bin < nbins; bin++) {
            sel(flt, sex, bin, yr) = invlogit(-slp2 * (x(bin) - inf2));
          }
          break;
        }

        case 5: { // Hake: log selectivity is the cumulative sum from bin_first_selected
          Type cum_sum = 0;
          for (int bin = bin_first_selected(flt); bin < nbins; bin++) {
            if (bin < n_sel_bins) cum_sum += sel_coff(flt, sex, bin) + sel_coff_dev(flt, sex, bin, yr);
            sel(flt, sex, bin, yr) = cum_sum;
          }
          // Scaled to 1 per sex here, so Sel_norm_scope is inert on Hake (issue #236).
          Type log_ref = -1e10;
          if (sel_norm_bin1(flt) >= 0 && sel_norm_bin2(flt) < 0) {
            log_ref = sel(flt, sex, sel_norm_bin1(flt), yr);
          } else {
            for (int bin = 0; bin < nbins; bin++) log_ref = max2(log_ref, sel(flt, sex, bin, yr));
          }
          for (int bin = 0; bin < nbins; bin++) {
            sel(flt, sex, bin, yr) = exp(sel(flt, sex, bin, yr) - log_ref);
          }
          break;
        }

        case 6: // 2DAR1
        case 7: // 3DAR1
          for (int bin = 0; bin < nbins; bin++) {
            int b = (bin < n_sel_bins) ? bin : (n_sel_bins - 1);   // plateau past N_sel_bins
            sel(flt, sex, bin, yr) = invlogit(sel_coff(flt, sex, b) + sel_coff_dev(flt, sex, b, yr));
          }
          break;

        case 8: { // DoubleNormal
          // Slots: sel_inf(0) peak, log_sel_slp(0)/(1) log SDs of the two limbs, sel_inf(1)
          // logit of the right-tail floor (SS3's P6; a floor of 0 is a full dome, 1 a logistic).
          Type peak        = inflection(0);
          Type sigma_asc   = slope(0);
          Type sigma_desc  = slope(1);
          Type right_floor = invlogit(inflection(1));
          for (int bin = 0; bin < nbins; bin++) {
            Type d          = x(bin) - peak;
            Type w          = invlogit(20.0 * d);   // smooth switch between the limbs
            Type asc_gauss  = exp(-0.5 * pow(d / sigma_asc,  2.0));
            Type desc_gauss = exp(-0.5 * pow(d / sigma_desc, 2.0));
            sel(flt, sex, bin, yr) = (1.0 - w) * asc_gauss
                                     + w * (right_floor + (1.0 - right_floor) * desc_gauss);
          }
          break;
        }

        case 9: { // NonParametricPM
          // Random walk on the log curve; it restarts from the base sel_coff at or before
          // the start year, as AMAK sets a survey's base selectivity once.
          bool from_base = (yr <= flt_sel_start_yr(flt));
          for (int bin = 0; bin < nbins; bin++) {
            if (bin < n_sel_bins) {
              Type prev = from_base ? sel_coff(flt, sex, bin) : log_walk(sex, bin, yr - 1);
              log_walk(sex, bin, yr) = prev + sel_coff_dev(flt, sex, bin, yr);
            } else {
              log_walk(sex, bin, yr) = log_walk(sex, n_sel_bins - 1, yr);
            }
          }
          // Unselected bins are held at 0 rather than walking (AMAK's fixed first-bin value).
          for (int bin = 0; bin < bin_first_selected(flt); bin++) log_walk(sex, bin, yr) = 0.0;

          vector<Type> walk(nbins);
          for (int bin = 0; bin < nbins; bin++) walk(bin) = log_walk(sex, bin, yr);
          Type centre = log_mean_exp(walk);
          for (int bin = 0; bin < nbins; bin++) log_walk(sex, bin, yr) -= centre;

          // The curve: the centred walk held flat above the cap bin, centred again.
          int cap = flt_sel_cap_bin(flt);
          vector<Type> capped(nbins);
          for (int bin = 0; bin < nbins; bin++) {
            capped(bin) = (cap >= 0 && bin > cap) ? log_walk(sex, cap, yr) : log_walk(sex, bin, yr);
          }
          Type recentre = log_mean_exp(capped);
          for (int bin = 0; bin < nbins; bin++) {
            log_non_par_sel(flt, sex, bin, yr) = capped(bin) - recentre;
            non_par_sel(flt, sex, bin, yr) = exp(log_non_par_sel(flt, sex, bin, yr));
            sel(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
          }
          break;
        }

        case 11: { // LogisticPM
          // AMAK evaluates at mid-age (j + 0.5), so x is bin + 1.5; the shift cannot go
          // into a50, whose deviate is multiplicative (so not inflection(0)). Not normalized.
          Type slp = slope(0);
          Type inf = (sel_inf(0, flt, sex) + sel_inf_off_nat(0, flt, sex, yr))
                     * exp(sel_inf_dev(0, flt, sex, yr) + sel_inf_off(0, flt, sex, yr));
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? x(bin) : Type(bin + 1.5);
            sel(flt, sex, bin, yr) = invlogit(slp * (x_val - inf));
          }
          // Free first-bin log selectivity, held in the unused descending slot sel_inf(1).
          Type log_s1 = (sel_inf(1, flt, sex) + sel_inf_off_nat(1, flt, sex, yr))
                        * exp(sel_inf_dev(1, flt, sex, yr) + sel_inf_off(1, flt, sex, yr));
          sel(flt, sex, 0, yr) = exp(log_s1);
          break;
        }

        case 15: { // DoubleNormalSS3 (SS_selex.tpl)
          // SS3's six parameters: P1 peak, P2 top width (logit), P3/P4 limb widths (log),
          // P5/P6 initial and final selectivity (logit). Each is (base + identity offset) x
          // exp(log offset): an SS3 block, then annual deviates, in SS3's own order.
          vector<Type> P(6);
          for (int k = 0; k < 6; k++) {
            P(k) = (sel_dn6(k, flt, sex) + sel_dn6_off_nat(k, flt, sex, yr)) * exp(sel_dn6_off(k, flt, sex, yr));
          }
          Type peak      = P(0);
          Type upselex   = exp(P(2));
          Type downselex = exp(P(3));
          // startbin: first population bin at or above the first data bin
          // (SS_readdata_330.tpl:1642-1646, used at SS_selex.tpl:142).
          int startbin = 0;
          if (is_length_based) {
            while (startbin < nbins - 1 && lengths_pop(sp, startbin) < lengths(sp, 0)) startbin++;
          }
          Type x_first = x(startbin);
          Type x_last  = x(nbins - 1);
          Type peak2   = peak + binwidth + (0.99 * x_last - peak - binwidth) / (1.0 + exp(-P(1)));
          Type point1  = invlogit(P(4));
          Type point2  = invlogit(P(5));
          Type t1min   = exp(-square(x_first - peak) / upselex);
          Type t2min   = exp(-square(x_last - peak2) / downselex);
          for (int bin = startbin; bin < nbins; bin++) {
            Type t1    = x(bin) - peak;
            Type t2    = x(bin) - peak2;
            Type join1 = invlogit(20.0 * t1 / (1.0 + CppAD::abs(t1)));
            Type join2 = invlogit(20.0 * t2 / (1.0 + CppAD::abs(t2)));
            Type asc   = (sel_dn6_ends(flt, 0) == 1)
                         ? point1 + (1.0 - point1) * (exp(-square(t1) / upselex) - t1min) / (1.0 - t1min)
                         : exp(-square(t1) / upselex);
            Type dsc   = (sel_dn6_ends(flt, 1) == 1)
                         ? 1.0 + (point2 - 1.0) * (exp(-square(t2) / downselex) - 1.0) / (t2min - 1.0)
                         : exp(-square(t2) / downselex);
            sel(flt, sex, bin, yr) = asc * (1.0 - join1) + join1 * ((1.0 - join2) + dsc * join2);
          }
          // Quadratic ramp below startbin (SS_selex.tpl:197-201).
          if (is_length_based && startbin > 0) {
            Type sel_start = sel(flt, sex, startbin, yr);
            for (int bin = 0; bin < startbin; bin++) {
              sel(flt, sex, bin, yr) = square(x(bin) / x(startbin)) * sel_start;
            }
          }
          break;
        }
        }

      } // End sex
    } // End year

    // --- 2b. PER-SEX APICAL HEIGHT (only when a linkage names `apical`) ---
    // Skipped otherwise, so other models are unchanged. The normalizer must pool across
    // sexes or be off for it to survive; .check_sel_apical_rows() refuses the rest.
    if (sel_apical_on == 1) {
      for (int yr = 0; yr < nyrs_hind; yr++) {
        for (int sex = 0; sex < nsex(sp); sex++) {
          Type scale = exp(log_sel_apical(flt, sex) + sel_apical_off(flt, sex, yr))
                       + sel_apical_off_nat(flt, sex, yr);
          for (int bin = 0; bin < nbins; bin++) sel(flt, sex, bin, yr) *= scale;
        }
      }
    }

    // --- 3. ZERO, NORMALIZE AND PROJECT ---
    normalize_and_project_selectivity(
      flt, nyrs_hind, nyrs, flt_spp, flt_sel_type, flt_sel_dim, bin_first_selected, nages,
      nlengths_pop, nsex, sel_norm_bin1, sel_norm_bin2, sel_norm_scope, sel);

    // --- 4. LENGTH TO AGE ---
    if (is_length_based) {
      int wtind = nspp * 2 + flt;   // growth_matrix: two entries per species, then one per fleet
      for (int yr = 0; yr < nyrs; yr++) {
        for (int sex = 0; sex < nsex(sp); sex++) {
          convert_length_selectivity(flt, sp, sex, yr, wtind, nages, nlengths_pop,
                                     sel_at_length, growth_matrix, sel_at_age);
        }
      }
    }
  } // End fleet
}

#endif
