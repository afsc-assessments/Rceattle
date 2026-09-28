#ifndef PREDATION_HPP
#define PREDATION_HPP


/**
 * @brief Transforms predator-prey preference parameters into vulnerability.
 *
 * Uses a multinomial logistic transformation to ensure that the sum of
 * predator-prey preferences and the vulnerability to "other food" sums to 1.
 * This is only applied to species using parametric suitability (suitMode > 0).
 * Adopted from Trijoulet et al. 2020 (https://github.com/vtrijoulet/Multisp_model_JAE/blob/master/MS_SSM.cpp).
 *
 * @param vulnerability       Matrix updated with bounded predator-prey vulnerability [nspp, nspp].
 * @param vulnerability_other Vector updated with the remaining vulnerability for "other food" [nspp].
 * @param nspp                Total number of species in the model.
 * @param suitMode            Vector specifying the suitability/preference mode for each species.
 * @param log_phi             Matrix of estimated log-preference parameters [nspp, nspp].
 */
template <class Type>
void calculate_vulnerability(matrix<Type> &vulnerability,
                             vector<Type> &vulnerability_other,
                             int nspp,
                             const vector<int> &suitMode,
                             matrix<Type> &log_phi) {

  // Reset arrays to ensure clean calculations
  vulnerability.setZero();
  vulnerability_other.setZero();

  for(int rsp = 0; rsp < nspp; rsp++) {
    if(suitMode(rsp) > 0) {

      Type sum_phi = 0; // Sum of predator-prey preference coefficients for multinomial transformation

      // 1. Calculate the sum of exponentiated preferences
      for(int ksp = 0; ksp < nspp; ksp++) {
        sum_phi += exp(log_phi(rsp, ksp));
      }

      // 2. Apply multinomial logistic transformation
      for(int ksp = 0; ksp < nspp; ksp++) {
        vulnerability(rsp, ksp) = exp(log_phi(rsp, ksp)) / (Type(1.0) + sum_phi);
      }

      // 3. Assign the remainder to "other food"
      vulnerability_other(rsp) = Type(1.0) - vulnerability.row(rsp).sum();
    }
  }
}

/**
 * @brief Calculates parametric suitability (Gamma, Lognormal, and Normal forms) based on predator/prey size ratios.
 *
 * This function merges length-based and weight-based suitability calculations into a single,
 * highly optimized pass. It calculates the relative size ratio of predators to prey and evaluates
 * that ratio against a specified probability distribution function to determine suitability.
 * The `suitMode` switch dictates both the distribution and the biological metric used:
 * - suitMode = 1: Gamma distribution (Length-based)       [NOT YET AVAILABLE -- blocked in data_check()]
 * - suitMode = 2: Gamma distribution (Weight-based)
 * - suitMode = 3: Lognormal distribution (Length-based)   [NOT YET AVAILABLE -- blocked in data_check()]
 * - suitMode = 4: Lognormal distribution (Weight-based)
 * - suitMode = 5: Normal distribution (Length-based)     [NOT YET AVAILABLE -- blocked in data_check()]
 * - suitMode = 6: Normal distribution (Weight-based)
 * Note: length-based modes (1, 3, 5) require a growth model to be active and are blocked
 * at the R level until a growth-model validation pass is added to data_check().
 *
 * @param suitability         Array updated with size-based suitability [nspp*max_sex, nspp*max_sex, max_age, max_age, nyrs].
 * @param suit_other          Array updated with the baseline vulnerability to "other food" [nspp, max_sex, max_age, nyrs].
 * @param nspp                Total number of species in the model.
 * @param nyrs                Total number of years evaluated.
 * @param nsex                Vector containing the number of sexes modeled for each species.
 * @param nages               Vector containing the number of age classes for each species.
 * @param suitMode            Vector specifying the parametric suitability mode for each predator species.
 * @param length_hat          Array of estimated lengths-at-age [nspp*2+n_flt, max_sex, max_age, nyrs].
 * @param weight_hat          Array of estimated weights-at-age [nspp*2+n_flt, max_sex, max_age, nyrs].
 * @param vulnerability       Matrix of bounded predator-prey vulnerability preferences [nspp, nspp].
 * @param vulnerability_other Vector of vulnerability preference for "other food" [nspp].
 * @param gam_a               Vector of distribution shape/mean parameters for predator size-selectivity [nspp].
 * @param gam_b               Vector of distribution scale/sd parameters for predator size-selectivity [nspp].
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

                  // 1. Calculate the size ratio metric based on suitMode
                  // - suitMode 1/2 = Gamma (Length/Weight)
                  // - suitMode 3/4 = Lognormal (Length/Weight)
                  // - suitMode 5/6 = Normal (Length/Weight)
                  if(smode == 1 || smode == 3) {
                    size_ratio = log(length_hat(rsp, r_sex, r_age, yr) / length_hat(ksp, k_sex, k_age, yr));
                  } else if (smode == 2 || smode == 4) {
                    size_ratio = log(weight_hat(wt_idx_rsp, r_sex, r_age, yr) / weight_hat(wt_idx_ksp, k_sex, k_age, yr)); // Log ratio of weights
                  } else if (smode == 5) {
                    size_ratio = length_hat(rsp, r_sex, r_age, yr) / length_hat(ksp, k_sex, k_age, yr);
                  } else if (smode == 6) {
                    size_ratio = weight_hat(wt_idx_rsp, r_sex, r_age, yr) / weight_hat(wt_idx_ksp, k_sex, k_age, yr);
                  }

                  // 2. Apply the distribution density function
                  if(size_ratio > 0) {
                    if(smode == 1 || smode == 2) {
                      // Gamma distribution
                      suitability(r_idx, k_idx, r_age, k_age, yr) = vulnerability(rsp, ksp) * dgamma(size_ratio, gam_a(rsp), gam_b(rsp)) /
                        dgamma((gam_a(rsp)-1) * gam_b(rsp), gam_a(rsp), gam_b(rsp));  // Scale to 0,1 by dividing by max
                    } else if (smode >= 3 && smode <= 6) {
                      // Normal/Lognormal distribution
                      suitability(r_idx, k_idx, r_age, k_age, yr) = vulnerability(rsp, ksp) * dnorm(size_ratio, gam_a(rsp), gam_b(rsp)) /
                        dnorm(gam_a(rsp), gam_a(rsp), gam_b(rsp));  // Scale to 0,1 by dividing by max
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
 * @brief Calculates MSVPA/empirical-based suitability based on stomach content and prey biomass.
 *
 * This function calculates suitability for predator species where suitMode == 0.
 * It first computes the ratio of stomach proportion to available prey biomass (stom_div_bio).
 * It then normalizes this ratio across all included prey and "other food" to determine suitability.
 * The suitability is averaged over a specified reference period (suit_styr to suit_endyr)
 * and held constant across all years. Two working arrays are allocated inside:
 * suma_suit, the sum of stom_div_bio over prey, and diet_prop_sum, the sum of
 * diet proportions, each [nspp, max_sex, max_age, nyrs].
 *
 * @param stom_div_bio         Array updated with the ratio of diet proportion to biomass [nspp*max_sex, nspp*max_sex, max_age, max_age, nyrs].
 * @param suitability          Array updated with the final calculated suitability [nspp*max_sex, nspp*max_sex, max_age, max_age, nyrs].
 * @param suit_other           Array updated with the remaining suitability for "other food" [nspp, max_sex, max_age, nyrs].
 * @param nspp                 Total number of species in the model.
 * @param nyrs                 Total number of years evaluated.
 * @param suit_styr            Start year index (per predator species) for the suitability averaging period.
 * @param suit_endyr           End year index (per predator species) for the suitability averaging period.
 * @param nyrs_suit            Number of years in the suitability averaging period, per predator species.
 * @param msmMode              Integer specifying the multi-species predation mode (e.g., 2 = Type III MSVPA).
 * @param nsex                 Vector containing the number of sexes modeled for each species.
 * @param nages                Vector containing the number of age classes for each species.
 * @param suitMode             Vector specifying the suitability mode for each predator species.
 * @param diet_prop            Array of observed or estimated diet proportions.
 * @param avgN_at_age          Array of average numbers-at-age.
 * @param weight_hat           Array of estimated weights-at-age.
 * @param other_food_diet_prop Array of estimated diet proportion belonging to "other food".
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

  // Initialize temporary TMB arrays
  int max_sex = imax(nsex);
  int max_age = imax(nages);
  array<Type> suma_suit(nspp, max_sex, max_age, nyrs);
  suma_suit.setZero();

  array<Type> diet_prop_sum(nspp, max_sex, max_age, nyrs);
  diet_prop_sum.setZero();

  stom_div_bio.setZero();


  for(int rsp = 0; rsp < nspp; rsp++) {

    if(suitMode(rsp) == 0) { // MSVPA based suitability

      // ======================================================================
      // PASS 1: Calculate stom_div_bio, suma_suit, and diet_prop_sum
      // ======================================================================
      for(int yr = 0; yr < nyrs; yr++) {
        for(int ksp = 0; ksp < nspp; ksp++) {
          int wt_idx_ksp = 2 * ksp;
          for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
            int k_idx = ksp + (nspp * k_sex);

            for(int k_age = 0; k_age < nages(ksp); k_age++) {

              // OPTIMIZATION 1: Calculate prey biomass denominator once per prey
              Type prey_wt = weight_hat(wt_idx_ksp, k_sex, k_age, yr);
              Type prey_N = avgN_at_age(ksp, k_sex, k_age, yr);
              if (msmMode == 2) { prey_N *= prey_N; } // Adjust for Type 3 MSVPA

              Type prey_denom = prey_wt * prey_N;

              // Only process if prey is available
              if (prey_denom > Type(0.0)) {
                Type inv_prey_denom = Type(1.0) / prey_denom; // Convert to fast multiplication

                for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
                  int r_idx = rsp + (nspp * r_sex);
                  for(int r_age = 0; r_age < nages(rsp); r_age++) {

                    suitability(r_idx, k_idx, r_age, k_age, yr) = 0.0;
                    Type d_prop = diet_prop(r_idx, k_idx, r_age, k_age, yr);

                    // Only process if predator actually ate this prey
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

      // ======================================================================
      // PASS 2: Normalize and calculate final Suitability
      // ======================================================================
      for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
        int r_idx = rsp + (nspp * r_sex);
        for(int r_age = 0; r_age < nages(rsp); r_age++) {

          // OPTIMIZATION 2: Pre-calculate the normalization denominator for this specific predator
          vector<Type> inv_pred_denom(nyrs);
          for(int yr = suit_styr(rsp); yr <= suit_endyr(rsp); yr++) {
            Type denom = suma_suit(rsp, r_sex, r_age, yr) + other_food_diet_prop(rsp, r_sex, r_age, yr);
            // `diet_prop_sum` gates whether this predator/age has diet data;
            // `denom` is what we divide by and can reach 0 or go negative
            // independently, so floor it.
            if (diet_prop_sum(rsp, r_sex, r_age, yr) > Type(0.0)) {
              inv_pred_denom(yr) = Type(1.0) / max2(denom, Type(1.0e-10));
            } else {
              inv_pred_denom(yr) = Type(0.0);
            }
          }

          Type total_suit_prey = Type(0.0); // Track total suitability assigned to modeled prey

          for(int ksp = 0; ksp < nspp; ksp++) {
            for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
              int k_idx = ksp + (nspp * k_sex);
              for(int k_age = 0; k_age < nages(ksp); k_age++) {

                Type suit_sum_yrs = Type(0.0);

                // Accumulate standardized suitability over the reference years
                for(int yr = suit_styr(rsp); yr <= suit_endyr(rsp); yr++) {
                  suit_sum_yrs += stom_div_bio(r_idx, k_idx, r_age, k_age, yr) * inv_pred_denom(yr);
                  // stom_div_bio(r_idx, k_idx, r_age, k_age, yr) / (suma_suit(rsp, r_sex, r_age, yr ) + other_food_diet_prop(rsp, r_sex, r_age, yr));
                }

                Type final_suit = suit_sum_yrs / Type(nyrs_suit(rsp));
                total_suit_prey += final_suit;

                // Fill all years
                for(int yr = 0; yr < nyrs; yr++) {
                  suitability(r_idx, k_idx, r_age, k_age, yr) = final_suit;
                }
              }
            }
          }

          // OPTIMIZATION 3: Assign "other food" suitability efficiently for all years
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
 * @brief Calculates available food, predation mortality (M2), and consumed biomass (MSVPA).
 *
 * This function calculates the core MSVPA predation dynamics in two passes.
 * Pass 1 calculates the total available food for each predator. Pass 2 allocates
 * predation mortality (M2) and consumption across prey species based on suitability
 * and available food. Dynamic reference points (dB0 and dBF) are calculated simultaneously.
 * * Highly optimized for TMB:
 * 1. Prey loops are inverted in Step 1 to prevent redundant prey biomass calculations.
 * 2. Predator-specific division scalars are hoisted out of the deepest loops in Step 2.
 * 3. Year loops are moved to the outermost level to maximize column-major cache locality.
 *
 * @param avail_food          Array updated with total available food [nspp, max_sex, max_age, nyrs].
 * @param avail_food_dB0      Array updated with available food at F=0.
 * @param avail_food_dBF      Array updated with available food at F=Ftarget.
 * @param M2_at_age           Array updated with calculated predation mortality [nspp, max_sex, max_age, nyrs].
 * @param M2_at_age_dB0       Array updated with predation mortality at F=0.
 * @param M2_at_age_dBF       Array updated with predation mortality at F=Ftarget.
 * @param M2_prop             Array updated with relative predation mortality.
 * @param B_eaten             Array updated with biomass of prey eaten by specific predator.
 * @param B_eaten_as_prey     Array updated with total biomass of prey eaten by all predators.
 * @param diet_prop_hat       Array updated with predicted diet proportions.
 * @param other_diet_prop_hat Array updated with predicted proportion of "other food" in diet.
 * @param nspp                Total number of species in the model.
 * @param nyrs                Total number of years evaluated.
 * @param nsex                Vector containing the number of sexes modeled for each species.
 * @param nages               Vector containing the number of age classes for each species.
 * @param msmMode             Integer specifying the multi-species predation mode (1 = Type II, 2 = Type III).
 * @param avgN_at_age         Array of average numbers-at-age.
 * @param avgN_at_age_dB0     Array of average numbers-at-age at F=0.
 * @param avgN_at_age_dBF     Array of average numbers-at-age at F=Ftarget.
 * @param suitability         Array of calculated suitability [nspp*max_sex, nspp*max_sex, max_age, max_age, nyrs].
 * @param suit_other          Array of suitability for "other food".
 * @param weight_hat          Array of estimated weights-at-age.
 * @param ration              Array of calculated individual consumption/ration (kg/yr).
 * @param other_food          Vector of base biomass for "other food" by predator.
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

  // ======================================================================
  // 0. SAFETY INITIALIZATION (Prevents ghost values during MSM iterations)
  // ======================================================================
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

  // ======================================================================
  // 1. Calculate Available Food (Prey loop on the OUTSIDE for efficiency)
  // ======================================================================
  for(int yr = 0; yr < nyrs; yr++) {
    for(int ksp = 0; ksp < nspp; ksp++) {
      int wt_idx_ksp = 2 * ksp;
      for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
        int k_idx = ksp + (nspp * k_sex);
        for(int k_age = 0; k_age < nages(ksp); k_age++) {

          // Calculate prey biomass ONCE per year/prey combination
          Type prey_N = avgN_at_age(ksp, k_sex, k_age, yr);
          Type prey_N_adj = (msmMode == 1) ? prey_N : (prey_N * prey_N);

          Type prey_N_dB0 = avgN_at_age_dB0(ksp, k_sex, k_age, yr);
          Type prey_N_adj_dB0 = (msmMode == 1) ? prey_N_dB0 : (prey_N_dB0 * prey_N_dB0);

          Type prey_N_dBF = avgN_at_age_dBF(ksp, k_sex, k_age, yr);
          Type prey_N_adj_dBF = (msmMode == 1) ? prey_N_dBF : (prey_N_dBF * prey_N_dBF);

          Type prey_wt = weight_hat(wt_idx_ksp, k_sex, k_age, yr);

          Type bio_avail     = prey_N_adj * prey_wt;
          Type bio_avail_dB0 = prey_N_adj_dB0 * prey_wt; // F=0 cross-calc
          Type bio_avail_dBF = prey_N_adj_dBF * prey_wt; // F=target cross-calc

          // Apply this pre-calculated biomass to all predators
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

    // Add "Other food" once per predator/year
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

  // ======================================================================
  // 2. Calculate Predation Mortality (M2) and Consumption
  // ======================================================================
  for(int yr = 0; yr < nyrs; yr++) {
    for(int rsp = 0; rsp < nspp; rsp++) {
      for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
        int r_idx = rsp + (nspp * r_sex);
        for(int r_age = 0; r_age < nages(rsp); r_age++) {

          Type avail = avail_food(rsp, r_sex, r_age, yr);

          if(avail > Type(0.0)) {
            // OPTIMIZATION: Hoist predator math and division scalars out of the prey loop
            Type inv_avail = Type(1.0) / avail;
            Type pred_rat = ration(rsp, r_sex, r_age, yr);
            Type pred_mort_scalar = avgN_at_age(rsp, r_sex, r_age, yr) * pred_rat * inv_avail;

            // Hoist dynamic reference point scalars
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

            // Assign other food proportion
            other_diet_prop_hat(rsp, r_sex, r_age, yr) = other_food(rsp) * suit_other(rsp, r_sex, r_age, yr) * inv_avail;

            for(int ksp = 0; ksp < nspp; ksp++) {
              int wt_idx_ksp = 2 * ksp;
              for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
                int k_idx = ksp + (nspp * k_sex);
                for(int k_age = 0; k_age < nages(ksp); k_age++) {

                  Type prey_N = avgN_at_age(ksp, k_sex, k_age, yr);
                  Type prey_wt = weight_hat(wt_idx_ksp, k_sex, k_age, yr);
                  Type suit = suitability(r_idx, k_idx, r_age, k_age, yr);

                  // Predation fractions
                  Type M2_frac = pred_mort_scalar * suit;
                  M2_at_age(ksp, k_sex, k_age, yr) += M2_frac;
                  M2_prop(r_idx, k_idx, r_age, k_age, yr) = M2_frac;

                  Type eaten = prey_N * prey_wt * M2_frac;
                  B_eaten(r_idx, k_idx, r_age, k_age, yr) = eaten;
                  B_eaten_as_prey(ksp, k_sex, k_age, yr) += eaten;

                  diet_prop_hat(r_idx, k_idx, r_age, k_age, yr) = prey_N * prey_wt * suit * inv_avail;

                  // Safely calculate dynamic reference points inside the prey loop
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






// 8.2. KINZEY PREDATION EQUATIONS
/*
 if(msmMode > 2) {

 // 8.2.3. Initialize counters
 Type Pred_ratio = 0.0;          // Predator ratio
 Type Prey_ratio = 0.0;          // Prey ratio
 Type NS_Z = 0.0;                // N(k,yr,a) * survival/Z = 0.0;
 Type Tmort = 0.0;               // Mortality on other
 Type Q_ksum_l = 0.0;            // Diet sum
 Type Term = 0.0;                // Linear adjustment for predation


 // 8.2.4. Calculate equilibrium N predators and prey in styr_pred for each species X age: FIXME: May want to have this be the final year of a projection!
 N_pred_eq.setZero();
 N_prey_eq.setZero();
 for(rsp = 0; rsp < nspp; rsp++) {
 for(ksp = 0; ksp < nspp; ksp++) {
 for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
 for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
 for(r_age = 0; r_age < nages(rsp); r_age++) {
 for(k_age = 0; k_age < nages(ksp); k_age++) {
 N_pred_eq(rsp, r_sex, r_age) += N_at_age(rsp, r_sex, r_age, 0) * suitability(r_idx, k_idx, r_age, k_age, 0); // Denominator of Eq. 17 Kinzey and Punt (2009) 1st year - FIXME: probably should use 2015
 N_prey_eq(ksp, k_sex, k_age) += N_at_age(ksp, k_sex, k_age, 0) * suitability(r_idx, k_idx, r_age, k_age, 0); // Denominator of Eq. 16 Kinzey and Punt (2009) 1st year - FIXME: probably should use 2015
 }
 }
 }
 }
 }
 }


 // 8.2.5. Calculate available prey and predator for each year
 N_pred_yrs.setZero();
 N_prey_yrs.setZero();
 for(yr = 0; yr < nyrs; yr++) {
 for(rsp = 0; rsp < nspp; rsp++) {
 for(ksp = 0; ksp < nspp; ksp++) {
 for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
 for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
 for(r_age = 0; r_age < nages(rsp); r_age++) {
 for(k_age = 0; k_age < nages(ksp); k_age++) {
 N_pred_yrs(rsp, r_sex, r_age, yr) += N_at_age(rsp, r_sex, r_age, yr) * suitability(r_idx, k_idx, r_age, k_age, yr); // Numerator of Eq. 17 Kinzey and Punt (2009) 1st year // FIXME: Use averageN?
 N_prey_yrs(ksp, k_sex, k_age, yr) += N_at_age(ksp, k_sex, k_age, yr) * suitability(r_idx, k_idx, r_age, k_age, yr); // Numerator of Eq. 16 Kinzey and Punt (2009) 1st year
 }
 }
 }
 }
 }
 }
 }


 // 8.2.6. Calculate predator functional response (Table 1 Kinzey and Punt (2009))
 for(rsp = 0; rsp < nspp; rsp++) {                          // Predator loop
 for(ksp = 0; ksp < (nspp + 1); ksp++) {                  // Prey loop
 Term = 1.0e-10 + H_1(rsp, ksp) * (Type(1) + H_1a(rsp) * H_1b(rsp) / (Type(r_age) + H_1b(rsp) + Type(1.0e-10))); // Eq. 15 Kinzey and Punt (2009)
 for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
 for(r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
 for(yr = 0; yr < nyrs; yr++) {                         // Year loop

 // Observed species
 if(ksp < nspp){
 for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
 for(k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop

 // Predator-prey ratios
 Pred_ratio = (N_pred_yrs(rsp, r_sex, r_age, yr) + Type(1.0e-10)) / (N_pred_eq(rsp, r_sex, r_age) + Type(1.0e-10)); // Eq. 17 Kinzey and Punt (2009): Predator biomass relative to equilibrium
 Prey_ratio = (N_prey_yrs(ksp, k_sex, k_age, yr) + Type(1.0e-10)) / (N_prey_eq(ksp, k_sex, k_age) + Type(1.0e-10)); // Eq. 16 Prey Kinzey and Punt (2009): biomass relative to equilibrium
 Pred_r(rsp, r_sex, r_age, yr) = Pred_ratio;
 Prey_r(ksp, k_sex, k_age, yr) = Prey_ratio;

 switch (msmMode) {
 case 3: // Holling Type I (linear)
 pred_resp(r_idx, k_idx, r_age, k_age, yr) = 1.0e-10 + Term;
 break;
 case 4: // Holling Type II
 pred_resp(r_idx, k_idx, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
 ( 1 + H_2(rsp, ksp) * Prey_ratio + 1.0e-10);
 break;
 case 5: // Holling Type III
 pred_resp(r_idx, k_idx, r_age, k_age, yr) = 1.0e-10 +
 Term * (1 + H_2(rsp, ksp)) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) /
 (1 + H_2(rsp, ksp) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) + 1.0e-10 );
 break;
 case 6: // predator interference
 pred_resp(r_idx, k_idx, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
 ( 1 + H_2(rsp, ksp) * Prey_ratio + H_3(rsp, ksp) * (Pred_ratio - 1) + 1.0e-10);
 break;
 case 7: // predator preemption
 pred_resp(r_idx, k_idx, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
 ( (1 + H_2(rsp, ksp) * Prey_ratio) * (1 + H_3(rsp, ksp) * (Pred_ratio - 1)) + Type(1.0e-10));
 break;
 case 8: // Hassell-Varley
 pred_resp(r_idx, k_idx, r_age, k_age, yr) = 1.0e-10 + Term * (2 + H_2(rsp, ksp) + 1.0e-10) /
 (1.0 + H_2(rsp, ksp) * Prey_ratio + pow((Prey_ratio + Type(1.0e-10)), H_4(rsp, ksp)) + 1.0e-10 );
 break;
 case 9: // Ecosim
 pred_resp(r_idx, k_idx, r_age, k_age, yr) = 1.0e-10 + Term /
 (1 + H_3(rsp, ksp) * (Pred_ratio - 1 + 1.0e-10));
 break;
 default:
 error("Invalid 'msmMode'");
 }
 }
 }
 }
 // "other food" is linear
 else{
 pred_resp(r_idx, (nspp * 2), r_age, 0, yr)  = 1.0e-10 + Term;
 } // end of r_ages, k_ages loop
 }   // =========================
 }
 }
 }
 }


 // 8.2.7  Predation mortality: Eq. 6 & 7 Kinzey and Punt (2009)
 M2_at_age.setZero();
 for(yr = 0; yr < nyrs; yr++) {
 for(rsp = 0; rsp  <  nspp; rsp++) {
 for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
 for(ksp = 0; ksp  <  nspp; ksp++) {
 for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
 for(r_age = 0; r_age < nages(rsp); r_age++) {
 for(k_age = 0; k_age < nages(ksp); k_age++) {
 M2_prop(r_idx, k_idx, r_age, k_age, yr) = pred_resp(r_idx, k_idx, r_age, k_age, yr)
 * suitability(r_idx, k_idx, r_age, k_age, yr) * N_at_age(rsp, r_sex, r_age, yr);
 M2_at_age(ksp, k_sex, k_age, yr) += M2_prop(r_idx, k_idx, r_age, k_age, yr);
 }
 }
 }
 }
 }
 }
 }

 // 8.2.9. Numbers and mass eaten (of modeled prey species and "other prey"); Equations 8 and 9 from Kinzey and Punt 2009
 for(yr = 0; yr < nyrs; yr++) {
 for(rsp = 0; rsp  <  nspp; rsp++) {
 for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
 for(r_age = 0; r_age  <  nages(rsp); r_age++) {
 for(ksp = 0; ksp  <  nspp+1; ksp++) {
 // Species included
 if(ksp < nspp) {
 for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
 for(k_age = 0; k_age < nages(ksp); k_age++) {
 wt_idx_ksp = 2 * ksp;
 N_eaten(r_idx, k_idx, r_age, k_age, yr) = M2_prop(r_idx, k_idx, r_age, k_age, yr) * avgN_at_age(ksp, k_sex, k_age, yr); // Eq. 8 Kinzey and Punt (2009)
 B_eaten(r_idx, k_idx, r_age, k_age, yr) = M2_prop(r_idx, k_idx, r_age, k_age, yr) * avgN_at_age(ksp, k_sex, k_age, yr) * prey_wt;
 B_eaten_as_pred(rsp, r_sex, r_age, yr) += B_eaten(r_idx, k_idx, r_age, k_age, yr);
 B_eaten_as_prey(ksp, k_sex, k_age, yr) += B_eaten(r_idx, k_idx, r_age, k_age, yr); // Results by all species: Eq. 11 Kinzey and Punt (2009)
 }
 }
 // Other food
 }else{
 B_eaten_as_pred(rsp, r_sex, r_age, yr)  += other_food(rsp) * (Type(1) - exp(-pred_resp(r_idx, nspp*2, r_age, 0, yr)  * N_at_age(rsp, r_sex, r_age, yr)));      // Eq. 11 Kinzey and Punt (2009)
 }
 }
 }
 }
 }
 }


 // 8.2.11. Predicted diet proportion as weight of prey-at-age in diet of predator-at-age (Eqn 15)
 // NOTE: Only including prey species in the model. Does not include other food
 for(rsp = 0; rsp  <  nspp; rsp++) {
 for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
 for(r_age = 0; r_age  <  nages(rsp); r_age++) {
 for(ksp = 0; ksp  <  nspp; ksp++) {
 for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
 for(k_age = 0; k_age < nages(ksp); k_age++) {
 for(yr = 0; yr < nyrs_hind; yr++) {
 diet_prop_hat(r_idx, k_idx, r_age, k_age, yr) = B_eaten(r_idx, k_idx, r_age, k_age, yr)/B_eaten_as_pred(rsp, r_sex, r_age, yr);
 }
 }
 }
 }
 }
 }
 }
 }


 // 8.2.12. Predict ration Eq. 13 Kinzey and Punt (2009)
 Type numer;
 Type denom;

 ration_hat_ave.setZero();
 ration_hat.setZero();
 for(rsp = 0; rsp < nspp; rsp++) {
 for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
 for(r_age = 0; r_age < nages(rsp); r_age++) {
 numer = 0;
 denom = 0;
 // Calculate year-specific values
 for(yr = 0; yr < nyrs; yr++) {
 if(avgN_at_age(rsp, r_sex, r_age, yr) > 0){
 ration_hat(rsp, r_sex, r_age, yr) = B_eaten_as_pred(rsp, r_sex, r_age, yr)/avgN_at_age(rsp, r_sex, r_age, yr); // NOTE: Divide by 365 to make into daily ration
 }
 }

 // Find average over hindcast// FIXME: suit_years?
 for(yr = 0; yr < nyrs_hind; yr++) {
 numer += B_eaten_as_pred(rsp, r_sex, r_age, yr); // NOTE: Divide by 365 to make into daily ration
 denom += avgN_at_age(rsp, r_sex, r_age, yr);
 }
 ration_hat_ave(rsp, r_sex, r_age) = numer / denom;
 }
 }
 }
 // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
 } // End 8.2. Kinzey predation
 // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
 */

#endif
