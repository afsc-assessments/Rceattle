#ifndef BIOENERGETICS_HPP
#define BIOENERGETICS_HPP

/**
 * @brief Temperature-dependent function of consumption (fT) for bioenergetics.
 *
 * Updates the fT matrix in place. Supports multiple bioenergetics models:
 * 1 = Exponential function (Stewart et al. 1983)
 * 2 = Warm-water species temperature dependence (Kitchell et al. 1977)
 * 3 = Cool/cold-water species temperature dependence (Thornton and Lessem 1978)
 * 4 = Constant (fT = 1.0)
 *
 * @param fT         Temperature functions [nspp, nyrs], updated.
 * @param nspp       Number of species.
 * @param nyrs       Years, hindcast plus projection.
 * @param Ceq        Consumption equation per species.
 * @param Cindex     Environmental index driving each species.
 * @param Qc         Q coefficient in equations 1-2, a temperature in equation 3.
 * @param Tcm        Thermal max.
 * @param Tco        Thermal optimum.
 * @param Tcl        Thermal limit for Thornton/Lessem.
 * @param CK1        Ascending limit for Thornton/Lessem.
 * @param CK4        Descending limit for Thornton/Lessem.
 * @param env_index  Environmental predictors, e.g. bottom temperature [nyrs, n_indices].
 */
template <class Type>
void calculate_temperature_function(matrix<Type> &fT,
                                    int nspp,
                                    int nyrs,
                                    const vector<int> &Ceq,
                                    const vector<int> &Cindex,
                                    const vector<Type> &Qc,
                                    const vector<Type> &Tcm,
                                    const vector<Type> &Tco,
                                    const vector<Type> &Tcl,
                                    const vector<Type> &CK1,
                                    const vector<Type> &CK4,
                                    const matrix<Type> &env_index) {

  // Kitchell and Thornton/Lessem intermediates, hoisted out to keep the AD tape small.
  Type Yc = 0.0;
  Type Zc = 0.0;
  Type Vc = 0.0;
  Type Xc = 0.0;
  Type G2 = 0.0;
  Type L2 = 0.0;
  Type G1 = 0.0;
  Type L1 = 0.0;
  Type Ka = 0.0;
  Type Kb = 0.0;

  for(int sp = 0; sp < nspp; sp++) {

    // Equations 1-3 read one column of env_index; equation 4 reads none. Cindex is
    // 0-based (ceattle.cpp decrements the 1-based input), and an out-of-range column
    // would read adjacent memory (the model is built without bounds checking), so with
    // no column to read consumption takes equation 4's fT = 1.
    bool has_env = (Cindex(sp) >= 0) && (Cindex(sp) < env_index.cols());
    int  eq      = has_env ? Ceq(sp) : 4;

    for(int yr = 0; yr < nyrs; yr++) {

      switch(eq) {
      case 1: // Exponential function from Stewart et al. 1983
        fT(sp, yr) = exp(Qc(sp) * env_index(yr, Cindex(sp)));
        break;

      case 2: // Warm-water species temperature dependence, Kitchell et al. 1977
        Yc = log(Qc(sp)) * (Tcm(sp) - Tco(sp) + 2.0);
        Zc = log(Qc(sp)) * (Tcm(sp) - Tco(sp));
        Vc = (Tcm(sp) - env_index(yr, Cindex(sp))) / (Tcm(sp) - Tco(sp));
        Xc = pow(Zc, 2) * pow((1.0 + pow((1.0 + 40.0 / Yc), 0.5)), 2) / 400.0;
        fT(sp, yr) = pow(Vc, Xc) * exp(Xc * (1.0 - Vc));
        break;

      case 3: // Cool/cold-water species temperature dependence, Thornton and Lessem 1978
        G2 = (1.0 / (Tcl(sp) - Tcm(sp))) * log((0.98 * (1.0 - CK4(sp))) / (CK4(sp) * 0.02));
        L2 = exp(G2 * (Tcl(sp) - env_index(yr, Cindex(sp))));
        Kb = (CK4(sp) * L2) / (1.0 + CK4(sp) * (L2 - 1.0));

        G1 = (1.0 / (Tco(sp) - Qc(sp))) * log((0.98 * (1.0 - CK1(sp))) / (CK1(sp) * 0.02));
        L1 = exp(G1 * (env_index(yr, Cindex(sp)) - Qc(sp)));
        Ka = (CK1(sp) * L1) / (1.0 + CK1(sp) * (L1 - 1.0));

        fT(sp, yr) = Ka * Kb;
        break;

      case 4: // Constant
        fT(sp, yr) = 1.0;
        break;

      default:
        // Ceq outside 1-4, which data_check() refuses. fT = 1 rather than the zero from
        // initialization, which would silently zero consumption and M2 with it.
        fT(sp, yr) = 1.0;
        break;
      }
    }
  }
}


/**
 * @brief Annual ration of a predator at age, in kg per predator per year.
 *
 * Updates `consumption_at_age` in place, from the bioenergetics equations (Ceq 1-3) or
 * from input consumption (Ceq 4).
 *
 * @param consumption_at_age Individual consumption, kg/yr [nspp, max_sex, max_age, nyrs].
 * @param nspp               Number of species.
 * @param nyrs               Years, hindcast plus projection.
 * @param nyrs_hind          Years in the hindcast.
 * @param nsex               Sexes modelled per species: 1 = combined, 2 = both.
 * @param nages              Age bins per species.
 * @param Ceq                Consumption equation per species.
 * @param CA                 Allometric mass function intercept (Cmax = CA * W^(1+CB)).
 * @param CB                 Allometric mass function slope.
 * @param fday               Foraging days per year, per species.
 * @param Pvalue             Proportion of maximum consumption (scales Cmax).
 * @param fT                 Temperature-dependence function [nspp, nyrs].
 * @param ration_data        Relative foraging rates (Ceq 1-3), input consumption (Ceq 4).
 * @param weight_hat         Weight-at-age [n_wt_indices, max_sex, max_age, nyrs].
 */
template <class Type>
void calculate_ration(array<Type> &consumption_at_age,
                      int nspp,
                      int nyrs,
                      int nyrs_hind,
                      const vector<int> &nsex,
                      const vector<int> &nages,
                      const vector<int> &Ceq,
                      const vector<Type> &CA,
                      const vector<Type> &CB,
                      const vector<Type> &fday,
                      const vector<Type> &Pvalue,
                      matrix<Type> &fT, // "const" doesnt work for TMB matrices/arrays
                      array<Type> &ration_data,
                      array<Type> &weight_hat) {

  for(int sp = 0; sp < nspp; sp++) {
    for(int sex = 0; sex < nsex(sp); sex++) {
      for(int age = 0; age < nages(sp); age++) {
        for(int yr = 0; yr < nyrs; yr++) {

          // Projection years reuse the last hindcast year's ration data.
          int yr_ind = (yr < nyrs_hind) ? yr : (nyrs_hind - 1);

          if(Ceq(sp) < 4){
            int wt_idx_pop = 2 * sp;

            // Cmax * fT * fday, grams/predator/yr; weight_hat is kg so x1000 gives grams.
            consumption_at_age(sp, sex, age, yr) = CA(sp) * pow(weight_hat(wt_idx_pop, sex, age, yr) * Type(1000.0), 1.0 + CB(sp)) * fT(sp, yr) * fday(sp);

            // Scale by Pvalue and the relative foraging rate, then grams to kg.
            consumption_at_age(sp, sex, age, yr) = consumption_at_age(sp, sex, age, yr) * Pvalue(sp) * ration_data(sp, sex, age, yr_ind) / 1000.0;
          }

          // Ceq 4 bypasses bioenergetics: input consumption, scaled by Pvalue.
          if(Ceq(sp) == 4){
            consumption_at_age(sp, sex, age, yr) = Pvalue(sp) * ration_data(sp, sex, age, yr_ind);
          }
        }
      }
    }
  }
}


#endif
