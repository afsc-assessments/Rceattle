// =====================================================================
// linkage.hpp
//
// Generic accumulator for the Rceattle long-format linkage table. The
// R side encodes the pooled `Rceattle_linkage_table` into parallel
// IVECTORs/VECTORs and a dense design matrix `linkage_X`. This header
// reads those, evaluates one (yearly) coefficient contribution per
// linkage row, and accumulates the result into a per-process tensor of
// year-varying offsets that the growth / mortality / recruitment
// modules can add to their linear-predictor parameters.
//
// Process and parameter codes are kept in lockstep with R/0-linkage_encode.R:
//   PROC: recruitment=0, M=1, growth=2, q=3, sel=4
//   PARAM (growth): K=0, L1=1, Linf=2, m=3 (natural-scale names)
//   PARAM (M): M1=0
//   ...
//
// Stratum sentinels: a 0 in `linkage_species`, `linkage_sex`, or
// `linkage_age_bin` means "applies to every level of that stratum".
// Otherwise the value is a 1-based index into 1..nspp / 1..nsex(sp) /
// 1..nages(sp).
// =====================================================================
#ifndef RCEATTLE_LINKAGE_HPP
#define RCEATTLE_LINKAGE_HPP

// Process codes (must match R/0-linkage_encode.R::LINKAGE_PROCESS_CODES)
#define RCEATTLE_PROC_RECRUIT 0
#define RCEATTLE_PROC_M       1
#define RCEATTLE_PROC_GROWTH  2
#define RCEATTLE_PROC_Q       3
#define RCEATTLE_PROC_SEL     4
#define RCEATTLE_PROC_COMP    5   // prior-only: DM composition-weighting overdispersion

// Number of growth parameters tracked in the offset tensor
// (log_K, log_L1, log_Linf, log_m).
#define RCEATTLE_N_GROWTH_PARAMS 4

// Number of recruitment parameters tracked in the offset tensor
// (R0, alpha, beta -- natural-scale names; stored as log_R0 etc. in TMB).
#define RCEATTLE_N_REC_PARAMS 3
#define RCEATTLE_REC_R0    0
#define RCEATTLE_REC_ALPHA 1
#define RCEATTLE_REC_BETA  2

/**
 * @brief Expand a stratum sentinel id into the half-open iteration range [lo, hi).
 *
 * A stratum id of 0 means "applies to every level of this stratum"; otherwise it
 * is a 1-based index. The range is clamped to the levels this species actually
 * has (e.g. a sex-2 row applied to a single-sex stock resolves to an empty range).
 *
 * @param id 1-based stratum id, or 0 for "all levels".
 * @param n_levels Number of levels this stratum has for the current species.
 * @param lo [out] Inclusive lower bound of the iteration range.
 * @param hi [out] Exclusive upper bound of the iteration range.
 */
inline void rceattle_stratum_range(int id, int n_levels, int& lo, int& hi) {
  lo = (id == 0) ? 0 : (id - 1);
  hi = (id == 0) ? n_levels : id;
  if (hi > n_levels) hi = n_levels;
}


// ---------------------------------------------------------------------
// The three accumulators below share one shape:
//
//   for each linkage row of this process, on this link scale:
//     expand the stratum sentinels, then add beta * X(yr, col) into the
//     offset tensor for every (stratum, year) the row applies to.
//
// `link_code` selects which rows a call consumes and therefore which
// tensor it fills: 1 = log (added inside the exp at the consume site),
// 0 = identity (added to the natural-scale value afterwards). The
// consumer combines them as
//   value_yr = exp(log_base + log_offset(yr)) + nat_offset(yr)
// so each process is called twice, once per scale, with the matching
// tensor. Years beyond `linkage_X.rows()` keep a zero offset: env_data
// need not span the projection horizon.
// ---------------------------------------------------------------------


/**
 * @brief Accumulate growth linkage offsets into the growth offset tensor.
 *
 * For each growth linkage row consumed on this link scale, expand the stratum
 * sentinels and add `beta * X(yr, col)` into `growth_offset[sp, sex, yr, param]`
 * for every (species, sex, year) the row applies to. Growth is not age-stratified.
 * See the block comment above for the two-scale (log / identity) consume convention.
 *
 * @param growth_offset [in,out] Offset tensor [nspp, max_nsex, nyrs, n_growth_params].
 * @param link_code Link scale to consume (1 = log, 0 = identity).
 * @param linkage_X Environmental covariate matrix; rows are years.
 * @param beta Per-row effect sizes (0-length = no-op).
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
 * @brief Accumulate natural-mortality (log_M1) linkage offsets.
 *
 * Adds `beta * X(yr, col)` into `M_offset[sp, sex, age_bin, yr]` for every
 * (species, sex, age_bin, year) each consumed M row applies to. Age-stratified
 * (unlike growth); only log_M1 is exposed as a target so far.
 *
 * @param M_offset [in,out] Offset tensor [nspp, max_nsex, max_age, nyrs].
 * @param link_code Link scale to consume (1 = log, 0 = identity).
 * @param linkage_X Environmental covariate matrix; rows are years.
 * @param beta Per-row effect sizes (0-length = no-op).
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
 * @brief Accumulate recruitment linkage offsets (R0 / alpha / beta).
 *
 * Adds `beta * X(yr, col)` into `rec_offset[sp, param, yr]` for every
 * (species, year) each consumed recruitment row applies to. Recruitment happens
 * at age 0 with no sex stratification on the parameter itself (sex_ratio splits
 * recruits downstream), so the sex / age_bin sentinels are accepted and ignored.
 *
 * @param rec_offset [in,out] Offset tensor [nspp, n_rec_params, nyrs].
 * @param link_code Link scale to consume (1 = log, 0 = identity).
 * @param linkage_X Environmental covariate matrix; rows are years.
 * @param beta Per-row effect sizes (0-length = no-op).
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
 * @brief Accumulate catchability (q) linkage offsets.
 *
 * Adds `beta * X(yr, col)` into `q_offset[fleet, yr]` for every (fleet, year)
 * each consumed q row applies to. Catchability is indexed by fleet rather than
 * species/sex/age (a fleet already implies its species through fleet_control),
 * so those sentinels are accepted and ignored; only q itself is exposed.
 *
 * @param q_offset [in,out] Offset matrix [n_flt, nyrs].
 * @param link_code Link scale to consume (1 = log, 0 = identity).
 * @param linkage_fleet Fleet stratum id per row (0 = all fleets).
 * @param linkage_X Environmental covariate matrix; rows are years.
 * @param beta Per-row effect sizes (0-length = no-op).
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
 * @brief Accumulate selectivity linkage offsets into three parameter-family tensors.
 *
 * Fills one of three offset tensors per row, routed by the row's param code, for
 * every (fleet, sex, year) the row applies to. Selectivity is indexed by fleet
 * and sex (mirroring is collapsed by the TMB map, so the parameter arrays are
 * effectively per-fleet). Called once per link scale, like the other processes.
 *
 * Param codes: `0`/`1` -> log_sel_slp[asc/desc] (`slp_offset`);
 *              `2`/`3` -> sel_inf[asc/desc]      (`inf_offset`);
 *              `4`     -> sel_coff, all bins     (`coff_offset`);
 *              `5`     -> log_sel_apical         (`apical_offset`);
 *              `6`-`11` -> sel_dn6[0..5], the six DoubleNormalSS3 (SS3
 *                         pattern 24) parameters (`dn6_offset`).
 *
 * @param slp_offset [in,out] Slope offsets [2, n_flt, max_sex, nyrs].
 * @param inf_offset [in,out] Inflection offsets [2, n_flt, max_sex, nyrs].
 * @param coff_offset [in,out] Nonparametric-coefficient offsets [n_flt, max_sex, n_sel_bins, nyrs].
 * @param apical_offset [in,out] Per-sex apical-height offsets [n_flt, max_sex, nyrs].
 * @param dn6_offset [in,out] DoubleNormalSS3 parameter offsets [6, n_flt, max_sex, nyrs].
 * @param link_code Link scale to consume (1 = log, 0 = identity).
 * @param linkage_X Environmental covariate matrix; rows are years.
 * @param beta Per-row effect sizes (0-length = no-op).
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

#endif
