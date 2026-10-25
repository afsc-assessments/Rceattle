// =====================================================================
// linkage.hpp
//
// Turns the linkage table (a covariate or time effect on a process parameter)
// into year-by-year offsets that recruitment, M, growth, catchability and
// selectivity add to their parameters. Each offset is beta * X(year, column).
//
// Example: build_M1(linkages = list(M1 = linkage_spec(~ temp))) makes two rows per
// species, an (Intercept) held at 0 and a `temp` slope beta. rceattle_apply_M_linkages()
// fills M_linkage_offset(sp, sex, age, yr) = beta * temp_yr, and ceattle.cpp section 6.1
// reads M1_yr = exp(log M1 + M_offset_yr). The equations for every link and
// the row / slot / group terms are at ceattle.cpp section 2.3b.
//
// Process and parameter codes are in lockstep with R/0-linkage_encode.R
// (LINKAGE_PROCESS_CODES / LINKAGE_PARAM_CODES); change one, change both.
// A stratum id of 0 (species, sex, age bin or fleet) means "every level".
// =====================================================================
#ifndef RCEATTLE_LINKAGE_HPP
#define RCEATTLE_LINKAGE_HPP

// Process codes (R/0-linkage_encode.R::LINKAGE_PROCESS_CODES)
#define RCEATTLE_PROC_RECRUIT 0
#define RCEATTLE_PROC_M       1
#define RCEATTLE_PROC_GROWTH  2
#define RCEATTLE_PROC_Q       3
#define RCEATTLE_PROC_SEL     4
#define RCEATTLE_PROC_COMP    5   // prior only, no offset: Dirichlet-multinomial theta

// Growth offsets: log_K, log_L1, log_Linf, log_m (codes 0-3). The SD codes
// 4 and 5 sit at or above this count, so the growth loop skips them.
#define RCEATTLE_N_GROWTH_PARAMS 4

// Recruitment offsets: log R0, log alpha, log Beta (the stock-recruit `rec_pars` columns),
// and R_init, the initial recruitment level, read at year 0 only.
#define RCEATTLE_N_REC_PARAMS 4
#define RCEATTLE_REC_R0    0
#define RCEATTLE_REC_ALPHA 1
#define RCEATTLE_REC_BETA  2
#define RCEATTLE_REC_R_INIT 3

/**
 * @brief Expand a stratum id into the loop range [lo, hi).
 *
 * Clamped to `n_levels`, so an id past the last level gives an empty range.
 *
 * @param id 1-based stratum id, or 0 for "all levels".
 * @param n_levels Number of levels the caller loops over.
 * @param lo,hi [out] First level, and one past the last level.
 */
inline void rceattle_stratum_range(int id, int n_levels, int& lo, int& hi) {
  lo = (id == 0) ? 0 : (id - 1);
  hi = (id == 0) ? n_levels : id;
  if (hi > n_levels) hi = n_levels;
}


// ---------------------------------------------------------------------
// The five functions below share one pattern: for each linkage row of this
// process whose link matches `link_code`, add beta * linkage_X(yr, col) to the
// offset for every stratum and year the row covers. Shared arguments:
//   linkage_X  covariate design matrix, one row per model year;
//   beta       one effect size per linkage row (length 0 = no linkages);
//   link_code  which rows this call adds:
//     0 = identity: beta * x is added to the parameter's natural-scale value
//         (to the stored logit or log for DoubleNormalSS3 P2-P6, DoubleNormal's
//         logit floor and LogisticPM's log first-bin selectivity);
//     1 = log: that value is multiplied by exp(beta * x);
//     2 = logit: reserved, refused in R;
//     3 = exponential, catchability only: (log q + any log-link offset) is
//         multiplied by exp(beta * x), SS3's environmental link type 1.
// Years beyond linkage_X keep an offset of 0, so env_data need not cover the
// projection. See vignette("environmental-linkages-and-priors").
// ---------------------------------------------------------------------


/**
 * @brief Add growth offsets to `growth_offset[sp, sex, yr, param]`.
 *
 * Growth is not age-specific, so the age-bin id is ignored. The SD codes
 * (sd_L1 = 4, sd_Linf = 5) take intercept rows only and are skipped here.
 *
 * @param growth_offset [in,out] Offsets [nspp, max_sex, nyrs, n_growth_params].
 * @param link_code,linkage_X,beta See the note above.
 */
template<class Type>
void rceattle_apply_growth_linkages(
    array<Type>&          growth_offset,
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   nspp,
    const vector<int>&    nsex,
    int                   nyrs)
{
  (void)linkage_age_bin;               // growth is not age-stratified
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_GROWTH) continue;
    if (linkage_link(i)    != link_code)            continue;

    int param = linkage_param(i);
    if (param < 0 || param >= RCEATTLE_N_GROWTH_PARAMS) continue;

    int xc = linkage_X_col(i);
    Type b = beta(i);

    int sp_lo, sp_hi;
    rceattle_stratum_range(linkage_species(i), nspp, sp_lo, sp_hi);

    for (int sp = sp_lo; sp < sp_hi; ++sp) {
      int sx_lo, sx_hi;
      rceattle_stratum_range(linkage_sex(i), nsex(sp), sx_lo, sx_hi);

      for (int sx = sx_lo; sx < sx_hi; ++sx) {
        for (int yr = 0; yr < yr_hi; ++yr) {
          growth_offset(sp, sx, yr, param) += b * linkage_X(yr, xc);
        }
      }
    }
  }
}


/**
 * @brief Add natural-mortality offsets to `M_offset[sp, sex, age, yr]`.
 *
 * Age-specific, unlike growth. M1 (code 0) is the only parameter.
 *
 * @param M_offset [in,out] Offsets [nspp, max_sex, max_age, nyrs].
 * @param link_code,linkage_X,beta See the note above.
 */
template<class Type>
void rceattle_apply_M_linkages(
    array<Type>&          M_offset,
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   nspp,
    const vector<int>&    nsex,
    const vector<int>&    nages,
    int                   nyrs)
{
  (void)linkage_param;                 // only log_M1 is exposed so far
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_M) continue;
    if (linkage_link(i)    != link_code)       continue;

    int xc = linkage_X_col(i);
    Type b = beta(i);

    int sp_lo, sp_hi;
    rceattle_stratum_range(linkage_species(i), nspp, sp_lo, sp_hi);

    for (int sp = sp_lo; sp < sp_hi; ++sp) {
      int sx_lo, sx_hi, ab_lo, ab_hi;
      rceattle_stratum_range(linkage_sex(i),     nsex(sp),  sx_lo, sx_hi);
      rceattle_stratum_range(linkage_age_bin(i), nages(sp), ab_lo, ab_hi);

      for (int sx = sx_lo; sx < sx_hi; ++sx) {
        for (int ab = ab_lo; ab < ab_hi; ++ab) {
          for (int yr = 0; yr < yr_hi; ++yr) {
            M_offset(sp, sx, ab, yr) += b * linkage_X(yr, xc);
          }
        }
      }
    }
  }
}


/**
 * @brief Add recruitment offsets to `rec_offset[sp, param, yr]`.
 *
 * Recruitment parameters are not sex- or age-specific (`sex_ratio` splits
 * recruits later), so those ids are ignored. R0 = 0, alpha = 1, beta = 2, R_init = 3.
 *
 * @param rec_offset [in,out] Offsets [nspp, n_rec_params, nyrs].
 * @param link_code,linkage_X,beta See the note above.
 */
template<class Type>
void rceattle_apply_recruitment_linkages(
    array<Type>&          rec_offset,
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   nspp,
    int                   nyrs)
{
  (void)linkage_sex;
  (void)linkage_age_bin;
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_RECRUIT) continue;
    if (linkage_link(i)    != link_code)             continue;

    int param = linkage_param(i);
    if (param < 0 || param >= RCEATTLE_N_REC_PARAMS) continue;

    int xc = linkage_X_col(i);
    Type b = beta(i);

    int sp_lo, sp_hi;
    rceattle_stratum_range(linkage_species(i), nspp, sp_lo, sp_hi);

    for (int sp = sp_lo; sp < sp_hi; ++sp) {
      for (int yr = 0; yr < yr_hi; ++yr) {
        rec_offset(sp, param, yr) += b * linkage_X(yr, xc);
      }
    }
  }
}


/**
 * @brief Add catchability offsets to `q_offset[fleet, yr]`.
 *
 * By fleet, which already implies the species, so species, sex and age ids are
 * ignored. q (code 0) is the only parameter.
 *
 * @param q_offset [in,out] Offsets [n_flt, nyrs].
 * @param linkage_fleet Fleet id per row (0 = all fleets).
 * @param link_code,linkage_X,beta See the note above.
 */
template<class Type>
void rceattle_apply_q_linkages(
    matrix<Type>&         q_offset,
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_fleet,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   n_flt,
    int                   nyrs)
{
  (void)linkage_param;                 // only q itself is exposed
  (void)linkage_species;
  (void)linkage_sex;
  (void)linkage_age_bin;
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_Q) continue;
    if (linkage_link(i)    != link_code)       continue;

    int xc = linkage_X_col(i);
    Type b = beta(i);

    int fl_lo, fl_hi;
    rceattle_stratum_range(linkage_fleet(i), n_flt, fl_lo, fl_hi);

    for (int flt = fl_lo; flt < fl_hi; ++flt) {
      for (int yr = 0; yr < yr_hi; ++yr) {
        q_offset(flt, yr) += b * linkage_X(yr, xc);
      }
    }
  }
}


/**
 * @brief Add selectivity offsets, by fleet, sex and year, to one of five tensors.
 *
 * Offsets are per fleet; fleets sharing a `Selectivity_index` still share one
 * parameter block through the map. The param code picks the tensor:
 *   0/1  log_sel_slp, ascending/descending  (`slp_offset`)
 *   2/3  sel_inf, ascending/descending      (`inf_offset`)
 *   4    sel_coff, every bin                (`coff_offset`)
 *   5    log_sel_apical                     (`apical_offset`)
 *   6-11 DoubleNormalSS3 P1-P6, SS3 pattern 24 (`dn6_offset`)
 *
 * @param slp_offset [in,out] Slope offsets [2, n_flt, max_sex, nyrs].
 * @param inf_offset [in,out] Inflection offsets [2, n_flt, max_sex, nyrs].
 * @param coff_offset [in,out] Nonparametric offsets [n_flt, max_sex, n_sel_bins, nyrs].
 * @param apical_offset [in,out] Apical-height offsets [n_flt, max_sex, nyrs].
 * @param dn6_offset [in,out] DoubleNormalSS3 offsets [6, n_flt, max_sex, nyrs].
 * @param link_code,linkage_X,beta See the note above.
 */
template<class Type>
void rceattle_apply_sel_linkages(
    array<Type>&          slp_offset,   // [2, n_flt, max_sex, nyrs]
    array<Type>&          inf_offset,   // [2, n_flt, max_sex, nyrs]
    array<Type>&          coff_offset,  // [n_flt, max_sex, n_sel_bins, nyrs]
    array<Type>&          apical_offset,// [n_flt, max_sex, nyrs]
    array<Type>&          dn6_offset,   // [6, n_flt, max_sex, nyrs]
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_fleet,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   n_flt,
    int                   max_sex,
    int                   n_sel_bins,
    int                   nyrs)
{
  (void)linkage_species;
  (void)linkage_age_bin;
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_SEL) continue;
    if (linkage_link(i)    != link_code)         continue;

    int param = linkage_param(i);
    int xc    = linkage_X_col(i);
    Type b    = beta(i);

    int fl_lo, fl_hi, sx_lo, sx_hi;
    rceattle_stratum_range(linkage_fleet(i), n_flt,   fl_lo, fl_hi);
    rceattle_stratum_range(linkage_sex(i),   max_sex, sx_lo, sx_hi);

    for (int flt = fl_lo; flt < fl_hi; ++flt) {
      for (int sx = sx_lo; sx < sx_hi; ++sx) {
        for (int yr = 0; yr < yr_hi; ++yr) {
          Type v = b * linkage_X(yr, xc);
          if (param == 0 || param == 1) {
            slp_offset(param, flt, sx, yr) += v;
          } else if (param == 2 || param == 3) {
            inf_offset(param - 2, flt, sx, yr) += v;
          } else if (param == 4) {
            for (int bin = 0; bin < n_sel_bins; ++bin) {
              coff_offset(flt, sx, bin, yr) += v;
            }
          } else if (param == 5) {
            apical_offset(flt, sx, yr) += v;
          } else if (param == 6) {
            dn6_offset(0, flt, sx, yr) += v;   // P1 peak
          } else if (param == 7) {
            dn6_offset(1, flt, sx, yr) += v;   // P2 top width (logit)
          } else if (param == 8) {
            dn6_offset(2, flt, sx, yr) += v;   // P3 ascending width (log)
          } else if (param == 9) {
            dn6_offset(3, flt, sx, yr) += v;   // P4 descending width (log)
          } else if (param == 10) {
            dn6_offset(4, flt, sx, yr) += v;   // P5 initial selectivity (logit)
          } else if (param == 11) {
            dn6_offset(5, flt, sx, yr) += v;   // P6 final selectivity (logit)
          }
        }
      }
    }
  }
}


/**
 * @brief The random-effect slots of one group, in slot order, which is time order
 * within a group (.assign_re_registry()).
 *
 * @param group_of_slot Each slot's group (ceattle.cpp's linkage_re_sigma: an index, not an sd).
 * @param grp The group.
 * @return The group's slot indices, ascending; empty if it has none.
 */
inline vector<int> rceattle_group_slots(const vector<int>& group_of_slot, int grp) {
  int len = 0;
  for (int s = 0; s < group_of_slot.size(); ++s) if (group_of_slot(s) == grp) len++;
  vector<int> slots(len);
  int j = 0;
  for (int s = 0; s < group_of_slot.size(); ++s) if (group_of_slot(s) == grp) slots(j++) = s;
  return slots;
}


/**
 * @brief Add the negative log density of one group's deviations to two accumulators.
 *
 * Structure 0 = IID, sd sigma; 1 = random walk, sd sigma on first differences (the
 * first deviation has no term of its own; build_map() holds it at its start, 0 by
 * default, so the base parameter carries the level); 2 = stationary AR1 with MARGINAL
 * sd sigma and correlation rho (the glmmTMB / Rogers QAR1 convention).
 *
 * @param structure 0, 1 or 2, as above.
 * @param re The group's deviations, in time order.
 * @param sigma The group's sd, natural scale (AR1: the marginal sd).
 * @param rho The AR1 correlation in (-1, 1); read only when structure = 2.
 * @param nll,nll_unweighted Cells of jnll_comp and unweighted_jnll_comp to add to.
 */
template<class Type>
void rceattle_re_group_nll(int structure, vector<Type> re, Type sigma, Type rho,
                           Type& nll, Type& nll_unweighted)
{
  if (structure == 1) {                    // random walk
    for (int t = 1; t < re.size(); ++t) {
      Type log_dens = dnorm(re(t) - re(t - 1), Type(0), sigma, true);
      nll            -= log_dens;
      nll_unweighted -= log_dens;
    }
  } else if (structure == 2) {             // AR1
    Type neg_log_dens = density::SCALE(density::AR1(rho), sigma)(re);   // a negative log density
    nll            += neg_log_dens;
    nll_unweighted += neg_log_dens;
  } else {                                 // IID
    Type log_dens = sum(dnorm(re, Type(0), sigma, true));
    nll            -= log_dens;
    nll_unweighted -= log_dens;
  }
}

#endif
