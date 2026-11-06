#ifndef PREDATION_HPP
#define PREDATION_HPP


/**
 * @brief Transforms predator-prey preference parameters into vulnerability.
 *
 * Multinomial logistic transform: each predator's prey preferences and its vulnerability
 * to "other food" sum to 1. Parametric suitability only (suitMode > 0). Adopted from
 * Trijoulet et al. 2020
 * (https://github.com/vtrijoulet/Multisp_model_JAE/blob/master/MS_SSM.cpp).
 *
 * @param vulnerability       Bounded predator-prey vulnerability [nspp, nspp], updated.
 * @param vulnerability_other Remaining vulnerability for "other food" [nspp], updated.
 * @param nspp                Number of species.
 * @param suitMode            Suitability/preference mode per species.
 * @param log_phi             Estimated log-preference parameters [nspp, nspp].
 */
template <class Type>
void calculate_vulnerability(matrix<Type> &vulnerability,
                             vector<Type> &vulnerability_other,
                             int nspp,
                             const vector<int> &suitMode,
                             matrix<Type> &log_phi) {

  vulnerability.setZero();
  vulnerability_other.setZero();

  for(int rsp = 0; rsp < nspp; rsp++) {
    if(suitMode(rsp) > 0) {

      Type sum_phi = 0;

      for(int ksp = 0; ksp < nspp; ksp++) {
        sum_phi += exp(log_phi(rsp, ksp));
      }

      for(int ksp = 0; ksp < nspp; ksp++) {
        vulnerability(rsp, ksp) = exp(log_phi(rsp, ksp)) / (Type(1.0) + sum_phi);
      }

      vulnerability_other(rsp) = Type(1.0) - vulnerability.row(rsp).sum();
    }
  }
}

/**
 * @brief Parametric suitability (Gamma, Lognormal, Normal) from size ratios.
 *
 * Evaluates the predator-to-prey size ratio against a density scaled to a peak of 1.
 * `suitMode` sets both the density and the size metric:
 * - Gamma: suitMode 1 (length-based) or 2 (weight-based)
 * - Lognormal: suitMode 3 (length-based) or 4 (weight-based)
 * - Normal: suitMode 5 (length-based) or 6 (weight-based)
 * The length-based modes 1, 3 and 5 are NOT YET AVAILABLE: they need an active growth
 * model, and data_check() blocks them until that validation is added.
 *
 * @param suitability         Size-based suitability, updated
 *                            [nspp*max_sex, nspp*max_sex, max_age, max_age, nyrs].
 * @param suit_other          Baseline vulnerability to "other food", updated;
 *                            [nspp, max_sex, max_age, nyrs] -- species NOT collapsed,
 *                            unlike suitability.
 * @param nspp                Number of species.
 * @param nyrs                Years, hindcast plus projection.
 * @param nsex                Sexes modelled per species: 1 = combined, 2 = both.
 * @param nages               Age bins per species.
 * @param suitMode            Parametric suitability mode per predator species.
 * @param length_hat          Lengths-at-age [nspp*2+n_flt, max_sex, max_age, nyrs].
 * @param weight_hat          Weights-at-age [nspp*2+n_flt, max_sex, max_age, nyrs].
 * @param vulnerability       Bounded predator-prey vulnerability [nspp, nspp].
 * @param vulnerability_other Vulnerability preference for "other food" [nspp].
 * @param gam_a               Shape/mean parameter for predator size-selectivity [nspp].
 * @param gam_b               Scale/sd parameter for predator size-selectivity [nspp].
 */
template <class Type>
void calculate_parametric_suitability(array<Type> &suitability,
                                      array<Type> &suit_other,
                                      int nspp, int nyrs,
                                      const vector<int> &nsex,
                                      const vector<int> &nages,
                                      const vector<int> &suitMode,
                                      array<Type> &length_hat,
                                      array<Type> &weight_hat,
                                      matrix<Type> &vulnerability,
                                      const vector<Type> &vulnerability_other,
                                      const vector<Type> &gam_a,
                                      const vector<Type> &gam_b) {

  for(int rsp = 0; rsp < nspp; rsp++) {
    int smode = suitMode(rsp);
    int wt_idx_rsp = 2 * rsp;

    // Parametric suitability modes (1 through 6)
    if(smode >= 1 && smode <= 6) {
      for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
        int r_idx = rsp + (nspp * r_sex);

        for(int r_age = 0; r_age < nages(rsp); r_age++) {
          for(int ksp = 0; ksp < nspp; ksp++) {
            int wt_idx_ksp = 2 * ksp;

            for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
              int k_idx = ksp + (nspp * k_sex);

              for(int k_age = 0; k_age < nages(ksp); k_age++) {
                for(int yr = 0; yr < nyrs; yr++) {

                  suit_other(rsp, r_sex, r_age, yr) = vulnerability_other(rsp);
                  Type size_ratio = 0.0;

                  // Row 2 * species (wt_idx_*) holds population length (cm) and weight. Modes 1-4
                  // take the log size ratio, 5-6 the raw one; a ratio <= 0 leaves suitability 0.
                  if(smode == 1 || smode == 3) {
                    size_ratio = log(length_hat(wt_idx_rsp, r_sex, r_age, yr) / length_hat(wt_idx_ksp, k_sex, k_age, yr));
                  } else if (smode == 2 || smode == 4) {
                    size_ratio = log(weight_hat(wt_idx_rsp, r_sex, r_age, yr) / weight_hat(wt_idx_ksp, k_sex, k_age, yr));
                  } else if (smode == 5) {
                    size_ratio = length_hat(wt_idx_rsp, r_sex, r_age, yr) / length_hat(wt_idx_ksp, k_sex, k_age, yr);
                  } else if (smode == 6) {
                    size_ratio = weight_hat(wt_idx_rsp, r_sex, r_age, yr) / weight_hat(wt_idx_ksp, k_sex, k_age, yr);
                  }

                  if(size_ratio > 0) {
                    if(smode == 1 || smode == 2) {
                      // Divided by the density at the mode, (gam_a-1)*gam_b.
                      suitability(r_idx, k_idx, r_age, k_age, yr) = vulnerability(rsp, ksp) * dgamma(size_ratio, gam_a(rsp), gam_b(rsp)) /
                        dgamma((gam_a(rsp)-1) * gam_b(rsp), gam_a(rsp), gam_b(rsp));
                    } else if (smode >= 3 && smode <= 6) {
                      // Divided by the density at the mean, gam_a.
                      suitability(r_idx, k_idx, r_age, k_age, yr) = vulnerability(rsp, ksp) * dnorm(size_ratio, gam_a(rsp), gam_b(rsp)) /
                        dnorm(gam_a(rsp), gam_a(rsp), gam_b(rsp));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


/**
 * @brief MSVPA/empirical suitability, from stomach content and prey biomass.
 *
 * For predators with suitMode = 0. stom_div_bio is the observed diet proportion divided
 * by prey biomass (numbers in thousands x weight in kg = mt; the prey term is squared
 * under msmMode = 2, Type III). That ratio is normalized over the modeled prey and "other
 * food", averaged over suit_styr..suit_endyr, and held constant across all years.
 *
 * A prey age whose average numbers approach zero with a positive diet proportion sends
 * stom_div_bio, and so the prey term, to infinity: that one cell takes the whole
 * normalized share, "other food" vanishes, available food collapses with it, so M2
 * diverges. The other-food term is floored at a fixed fraction of the prey term, bounding
 * suit_other below at floor/(1 + floor), leaving sum(suitability) + suit_other = 1 exact.
 * .warn_suitability_floor() reports when it binds.
 *
 * @param stom_div_bio         Ratio of diet proportion to biomass, updated;
 *                             [nspp*max_sex, nspp*max_sex, max_age, max_age, nyrs].
 * @param suitability          Final calculated suitability, updated.
 * @param suit_other           Remaining suitability for "other food", updated.
 * @param nspp                 Number of species.
 * @param nyrs                 Years, hindcast plus projection.
 * @param suit_styr            Start year index of the suitability window, per predator.
 * @param suit_endyr           End year index of the suitability window, per predator.
 * @param nyrs_suit            Years in the suitability window, per predator.
 * @param msmMode              Predation mode; 2 (Type III) squares the prey numbers.
 * @param nsex                 Sexes modelled per species: 1 = combined, 2 = both.
 * @param nages                Age bins per species.
 * @param suitMode             Suitability mode per predator species.
 * @param diet_prop            Observed diet proportions from organize_diet_obs().
 * @param avgN_at_age          Average numbers-at-age (thousands of fish).
 * @param weight_hat           Estimated weights-at-age (kg).
 * @param other_food_diet_prop Residual diet proportion for "other food", per mt of it.
 */
template <class Type>
void calculate_msvpa_suitability(array<Type> &stom_div_bio,
                                 array<Type> &suitability,
                                 array<Type> &suit_other,
                                 int nspp, int nyrs,
                                 const vector<int> &suit_styr, const vector<int> &suit_endyr, const vector<int> &nyrs_suit,
                                 int msmMode,
                                 const vector<int> &nsex,
                                 const vector<int> &nages,
                                 const vector<int> &suitMode,
                                 array<Type> &diet_prop, // const
                                 array<Type> &avgN_at_age,// const
                                 array<Type> &weight_hat,// const
                                 array<Type> &other_food_diet_prop) {// const

  // Working sums per predator species/sex/age/year: stom_div_bio over prey,
  // and diet proportion over prey.
  int max_sex = imax(nsex);
  int max_age = imax(nages);
  array<Type> suma_suit(nspp, max_sex, max_age, nyrs);
  suma_suit.setZero();

  array<Type> diet_prop_sum(nspp, max_sex, max_age, nyrs);
  diet_prop_sum.setZero();

  stom_div_bio.setZero();


  for(int rsp = 0; rsp < nspp; rsp++) {

    if(suitMode(rsp) == 0) { // MSVPA based suitability

      // Pass 1: diet proportion divided by prey biomass, and its per-predator-age sums.
      for(int yr = 0; yr < nyrs; yr++) {
        for(int ksp = 0; ksp < nspp; ksp++) {
          int wt_idx_ksp = 2 * ksp;
          for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
            int k_idx = ksp + (nspp * k_sex);

            for(int k_age = 0; k_age < nages(ksp); k_age++) {

              // prey_N is thousands of fish and prey_wt kg, so the product is mt.
              // Type III MSVPA (msmMode 2) squares prey_N: prey_denom is not a biomass.
              Type prey_wt = weight_hat(wt_idx_ksp, k_sex, k_age, yr);
              Type prey_N = avgN_at_age(ksp, k_sex, k_age, yr);
              if (msmMode == 2) { prey_N *= prey_N; }

              Type prey_denom = prey_wt * prey_N;

              if (prey_denom > Type(0.0)) {
                Type inv_prey_denom = Type(1.0) / prey_denom;

                for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
                  int r_idx = rsp + (nspp * r_sex);
                  for(int r_age = 0; r_age < nages(rsp); r_age++) {

                    suitability(r_idx, k_idx, r_age, k_age, yr) = 0.0;
                    Type d_prop = diet_prop(r_idx, k_idx, r_age, k_age, yr);

                    if (d_prop > Type(0.0)) {
                      Type sdb = d_prop * inv_prey_denom;
                      stom_div_bio(r_idx, k_idx, r_age, k_age, yr) = sdb;
                      suma_suit(rsp, r_sex, r_age, yr) += sdb;
                      diet_prop_sum(rsp, r_sex, r_age, yr) += d_prop;
                    }
                  }
                }
              }
            }
          }
        }
      }

      // Pass 2: normalize over the modeled prey plus "other food", then average the
      // standardized ratio over the suitability window.
      for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
        int r_idx = rsp + (nspp * r_sex);
        for(int r_age = 0; r_age < nages(rsp); r_age++) {

          vector<Type> inv_pred_denom(nyrs);
          for(int yr = suit_styr(rsp); yr <= suit_endyr(rsp); yr++) {
            Type prey_term  = suma_suit(rsp, r_sex, r_age, yr);
            // Other food is held at >= 1e-6 x the modelled-prey term, a unitless
            // ratio, so suit_other cannot reach 0 when a prey age's biomass does.
            Type other_term = max2(other_food_diet_prop(rsp, r_sex, r_age, yr),
                                   Type(1.0e-6) * prey_term);
            Type denom = prey_term + other_term;
            // `diet_prop_sum` gates whether this predator age has diet data; `denom` is
            // what we divide by and can reach 0 or go negative on its own, so floor it.
            if (diet_prop_sum(rsp, r_sex, r_age, yr) > Type(0.0)) {
              inv_pred_denom(yr) = Type(1.0) / max2(denom, Type(1.0e-10));
            } else {
              inv_pred_denom(yr) = Type(0.0);
            }
          }

          Type total_suit_prey = Type(0.0); // Total suitability assigned to modeled prey

          for(int ksp = 0; ksp < nspp; ksp++) {
            for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
              int k_idx = ksp + (nspp * k_sex);
              for(int k_age = 0; k_age < nages(ksp); k_age++) {

                Type suit_sum_yrs = Type(0.0);

                for(int yr = suit_styr(rsp); yr <= suit_endyr(rsp); yr++) {
                  suit_sum_yrs += stom_div_bio(r_idx, k_idx, r_age, k_age, yr) * inv_pred_denom(yr);
                }

                Type final_suit = suit_sum_yrs / Type(nyrs_suit(rsp));
                total_suit_prey += final_suit;

                for(int yr = 0; yr < nyrs; yr++) {
                  suitability(r_idx, k_idx, r_age, k_age, yr) = final_suit;
                }
              }
            }
          }

          // "Other food" takes whatever suitability the modeled prey leave.
          Type final_suit_other = Type(1.0) - total_suit_prey;
          for(int yr = 0; yr < nyrs; yr++) {
            suit_other(rsp, r_sex, r_age, yr) = final_suit_other;
          }
        }
      }
    }
  }
}


/**
 * @brief Available food, predation mortality (M2), and consumed biomass (MSVPA).
 *
 * Available food is suitability-weighted prey biomass plus "other food": prey numbers
 * (thousands of fish) x weight (kg) = mt, and other_food is mt. Under msmMode = 2
 * (Type III) the prey term is squared, so the denominator is not a biomass there. M2
 * divides each predator's consumption (avgN x ration) by that denominator, so a LARGER
 * available-food denominator gives LESS predation mortality. The dB0 (F = 0) and dBF
 * (F = Ftarget) arrays repeat the calculation for the dynamic reference points. Prey
 * biomass and the predator scalars are each computed once outside the deepest loop.
 *
 * @param avail_food          Available food (mt) [nspp, max_sex, max_age, nyrs], updated.
 * @param avail_food_dB0      Available food at F=0, updated.
 * @param avail_food_dBF      Available food at F=Ftarget, updated.
 * @param M2_at_age           Predation mortality (yr^-1), updated;
 *                            [nspp, max_sex, max_age, nyrs].
 * @param M2_at_age_dB0       Predation mortality at F=0, updated.
 * @param M2_at_age_dBF       Predation mortality at F=Ftarget, updated.
 * @param M2_prop             A predator age's contribution to prey M2, not a share of it.
 * @param B_eaten             Prey biomass eaten by a specific predator age (mt), updated.
 * @param B_eaten_as_prey     Prey biomass eaten by all predators (mt), updated.
 * @param diet_prop_hat       Predicted diet proportions, updated.
 * @param other_diet_prop_hat Predicted proportion of "other food" in diet, updated.
 * @param nspp                Number of species.
 * @param nyrs                Years, hindcast plus projection.
 * @param nsex                Sexes modelled per species: 1 = combined, 2 = both.
 * @param nages               Age bins per species.
 * @param msmMode             Predation mode (1 = Type II MSVPA, Holsman et al. 2015;
 *                            2 = Type III MSVPA, which squares the prey term).
 * @param avgN_at_age         Average numbers-at-age (thousands of fish).
 * @param avgN_at_age_dB0     Average numbers-at-age at F=0.
 * @param avgN_at_age_dBF     Average numbers-at-age at F=Ftarget.
 * @param suitability         Calculated suitability.
 * @param suit_other          Suitability for "other food".
 * @param weight_hat          Estimated weights-at-age (kg).
 * @param ration              Individual consumption (kg/predator/yr); consumption_at_age.
 * @param other_food          "Other food" biomass by predator (mt).
 */
template <class Type>
void calculate_msvpa_predation(array<Type> &avail_food,
                               array<Type> &avail_food_dB0,
                               array<Type> &avail_food_dBF,
                               array<Type> &M2_at_age,
                               array<Type> &M2_at_age_dB0,
                               array<Type> &M2_at_age_dBF,
                               array<Type> &M2_prop,
                               array<Type> &B_eaten,
                               array<Type> &B_eaten_as_prey,
                               array<Type> &diet_prop_hat,
                               array<Type> &other_diet_prop_hat,
                               int nspp, int nyrs,
                               const vector<int> &nsex,
                               const vector<int> &nages,
                               int msmMode,
                               array<Type> &avgN_at_age, // const
                               array<Type> &avgN_at_age_dB0, // const
                               array<Type> &avgN_at_age_dBF, // const
                               array<Type> &suitability, // const
                               array<Type> &suit_other, // const
                               array<Type> &weight_hat, // const
                               array<Type> &ration,
                               const vector<Type> &other_food) {

  // Predation is re-entered on every multispecies iteration (niter). M2, available food
  // and B_eaten_as_prey accumulate with +=; the rest are written only where avail > 0.
  avail_food.setZero();
  avail_food_dB0.setZero();
  avail_food_dBF.setZero();
  M2_at_age.setZero();
  M2_at_age_dB0.setZero();
  M2_at_age_dBF.setZero();
  M2_prop.setZero();
  B_eaten.setZero();
  B_eaten_as_prey.setZero();
  diet_prop_hat.setZero();
  other_diet_prop_hat.setZero();

  // 1. Available food: suitability-weighted prey biomass (mt), then "other food".
  for(int yr = 0; yr < nyrs; yr++) {
    for(int ksp = 0; ksp < nspp; ksp++) {
      int wt_idx_ksp = 2 * ksp;
      for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
        int k_idx = ksp + (nspp * k_sex);
        for(int k_age = 0; k_age < nages(ksp); k_age++) {

          // msmMode 2 (Type III) squares the prey numbers; "other food" is never squared.
          Type prey_N = avgN_at_age(ksp, k_sex, k_age, yr);
          Type prey_N_adj = (msmMode == 1) ? prey_N : (prey_N * prey_N);

          Type prey_N_dB0 = avgN_at_age_dB0(ksp, k_sex, k_age, yr);
          Type prey_N_adj_dB0 = (msmMode == 1) ? prey_N_dB0 : (prey_N_dB0 * prey_N_dB0);

          Type prey_N_dBF = avgN_at_age_dBF(ksp, k_sex, k_age, yr);
          Type prey_N_adj_dBF = (msmMode == 1) ? prey_N_dBF : (prey_N_dBF * prey_N_dBF);

          Type prey_wt = weight_hat(wt_idx_ksp, k_sex, k_age, yr);

          Type bio_avail     = prey_N_adj * prey_wt;
          Type bio_avail_dB0 = prey_N_adj_dB0 * prey_wt; // dynamic B0, F = 0
          Type bio_avail_dBF = prey_N_adj_dBF * prey_wt; // dynamic BF, F = Ftarget

          for(int rsp = 0; rsp < nspp; rsp++) {
            for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
              int r_idx = rsp + (nspp * r_sex);
              for(int r_age = 0; r_age < nages(rsp); r_age++) {
                Type suit = suitability(r_idx, k_idx, r_age, k_age, yr);
                avail_food(rsp, r_sex, r_age, yr)     += suit * bio_avail;
                avail_food_dB0(rsp, r_sex, r_age, yr) += suit * bio_avail_dB0;
                avail_food_dBF(rsp, r_sex, r_age, yr) += suit * bio_avail_dBF;
              }
            }
          }
        }
      }
    }

    // "Other food" (mt) is prey outside the model, added once per predator age and year.
    for(int rsp = 0; rsp < nspp; rsp++) {
      for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
        for(int r_age = 0; r_age < nages(rsp); r_age++) {
          Type other_f = other_food(rsp) * suit_other(rsp, r_sex, r_age, yr);
          avail_food(rsp, r_sex, r_age, yr)     += other_f;
          avail_food_dB0(rsp, r_sex, r_age, yr) += other_f;
          avail_food_dBF(rsp, r_sex, r_age, yr) += other_f;
        }
      }
    }
  }

  // 2. Predation mortality (M2) and the prey biomass eaten.
  for(int yr = 0; yr < nyrs; yr++) {
    for(int rsp = 0; rsp < nspp; rsp++) {
      for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
        int r_idx = rsp + (nspp * r_sex);
        for(int r_age = 0; r_age < nages(rsp); r_age++) {

          Type avail = avail_food(rsp, r_sex, r_age, yr);

          if(avail > Type(0.0)) {
            // avgN (thousands of fish) x ration (kg/predator/yr) is mt consumed, over
            // available food in mt.
            Type inv_avail = Type(1.0) / avail;
            Type pred_rat = ration(rsp, r_sex, r_age, yr);
            Type pred_mort_scalar = avgN_at_age(rsp, r_sex, r_age, yr) * pred_rat * inv_avail;

            Type pred_mort_scalar_dB0 = Type(0.0);
            Type avail_dB0 = avail_food_dB0(rsp, r_sex, r_age, yr);
            if(avail_dB0 > Type(0.0)) {
              pred_mort_scalar_dB0 = (avgN_at_age_dB0(rsp, r_sex, r_age, yr) * pred_rat) / avail_dB0;
            }

            Type pred_mort_scalar_dBF = Type(0.0);
            Type avail_dBF = avail_food_dBF(rsp, r_sex, r_age, yr);
            if(avail_dBF > Type(0.0)) {
              pred_mort_scalar_dBF = (avgN_at_age_dBF(rsp, r_sex, r_age, yr) * pred_rat) / avail_dBF;
            }

            other_diet_prop_hat(rsp, r_sex, r_age, yr) = other_food(rsp) * suit_other(rsp, r_sex, r_age, yr) * inv_avail;

            for(int ksp = 0; ksp < nspp; ksp++) {
              int wt_idx_ksp = 2 * ksp;
              for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
                int k_idx = ksp + (nspp * k_sex);
                for(int k_age = 0; k_age < nages(ksp); k_age++) {

                  Type prey_N = avgN_at_age(ksp, k_sex, k_age, yr);
                  Type prey_wt = weight_hat(wt_idx_ksp, k_sex, k_age, yr);
                  Type suit = suitability(r_idx, k_idx, r_age, k_age, yr);

                  Type M2_frac = pred_mort_scalar * suit;
                  M2_at_age(ksp, k_sex, k_age, yr) += M2_frac;
                  M2_prop(r_idx, k_idx, r_age, k_age, yr) = M2_frac;

                  Type eaten = prey_N * prey_wt * M2_frac;
                  B_eaten(r_idx, k_idx, r_age, k_age, yr) = eaten;
                  B_eaten_as_prey(ksp, k_sex, k_age, yr) += eaten;

                  diet_prop_hat(r_idx, k_idx, r_age, k_age, yr) = prey_N * prey_wt * suit * inv_avail;

                  M2_at_age_dB0(ksp, k_sex, k_age, yr) += pred_mort_scalar_dB0 * suit;
                  M2_at_age_dBF(ksp, k_sex, k_age, yr) += pred_mort_scalar_dBF * suit;
                }
              }
            }
          }
        }
      }
    }
  }
}



#endif
