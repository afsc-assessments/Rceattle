#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

/**
 * @file likelihood.hpp
 * @brief Likelihood math for section 13 of ceattle.cpp: the composition row
 * density and linkage prior families that several terms share, and the
 * composition plus-group fold. The loops over observations, and the SIMULATE
 * draws beside each density, stay in ceattle.cpp.
 */


/**
 * @brief Composition tail accumulation (AFSC ac_yng / ac_old).
 *
 * Bins below `yng` are added into `yng` and bins above `old` into `old`, and only
 * bins yng..old are kept. Joint-sex rows are folded per sex block, so a fold never
 * crosses from females into males. yng = 1 with old >= nbins is a no-op.
 *
 * @param obs,hat Observed and expected proportions; replaced by the folded vectors.
 * @param n_comp Number of bins in the row; replaced by the folded count.
 * @param nbins_blk Bins per sex block (ages or lengths).
 * @param nblk Sex blocks in the row: 2 for joint-sex comps, else 1.
 * @param yng,old First and last bin kept, 1-based ordinals on the fleet's comp dimension.
 */
template<class Type>
void fold_comp_tails(vector<Type>& obs, vector<Type>& hat, int& n_comp,
                     int nbins_blk, int nblk, int yng, int old)
{
  if (yng < 1) yng = 1;
  if (old < 1 || old > nbins_blk) old = nbins_blk;
  if (yng > old) yng = old;   // data_check() validates; keeps 1 <= yng <= old
  if (!(yng > 1 || old < nbins_blk)) return;   // no fold needed

  int nkeep = old - yng + 1;
  vector<Type> obs_acc(nblk * nkeep); obs_acc.setZero();
  vector<Type> hat_acc(nblk * nkeep); hat_acc.setZero();
  for (int b = 0; b < nblk; b++) {
    for (int j = 0; j < nbins_blk; j++) {
      int tgt = j;
      if (j < yng - 1) tgt = yng - 1;
      if (j > old - 1) tgt = old - 1;
      obs_acc(b * nkeep + tgt - (yng - 1)) += obs(b * nbins_blk + j);
      hat_acc(b * nkeep + tgt - (yng - 1)) += hat(b * nbins_blk + j);
    }
  }
  obs = obs_acc;
  hat = hat_acc;
  n_comp = nblk * nkeep;
}


/**
 * @brief Negative log-likelihood of one composition row: an age or length
 * composition, or one length bin of conditional age-at-length (CAAL).
 *
 * Likelihood types:
 *  - -1 AFSC/AMAK multinomial pseudo-likelihood, log L = sum N (o + c) log((p + c) / (o + c)).
 *       Not a density, so OSA residualizes it as a multinomial.
 *  -  0 Multinomial, via the conditional-binomial form in comp_osa.hpp.
 *  -  1 Dirichlet-multinomial with concentration N * p * theta (Thorson et al. 2017).
 *
 * The data weight multiplies the multinomial log-likelihood (types -1 and 0); the
 * Dirichlet-multinomial carries its weight as theta instead. `unweighted_nll` gets the
 * row without either, in ordinary fitting only.
 *
 * @param ll_type Likelihood type, -1, 0 or 1.
 * @param osa_mode 0 = ordinary fitting; 1 or 2 = OSA build (2 adds the CDF terms).
 * @param obs_prop,hat_prop Observed and expected proportions, after any tail fold.
 * @param offset Constant c added to both, to keep log(0) out (comp_offset).
 * @param n_samp Input sample size N.
 * @param weight Data weight on the multinomial log-likelihood.
 * @param dm_theta Dirichlet-multinomial theta, exp() of the log weight.
 * @param obsvec,keep,osa_start OSA observation vector, its keep indicator, and this
 *   row's first position in them (< 0 = row not in obsvec).
 * @param type_name Data column named in the error for an unknown type.
 * @param nll,unweighted_nll Cells of jnll_comp and unweighted_jnll_comp to add to.
 */
template<class Type>
void comp_row_nll(int ll_type, int osa_mode,
                  vector<Type> obs_prop, vector<Type> hat_prop,
                  Type offset, Type n_samp, Type weight, Type dm_theta,
                  vector<Type>& obsvec, data_indicator<vector<Type>, Type>& keep,
                  int osa_start, const char* type_name,
                  Type& nll, Type& unweighted_nll)
{
  int n_comp = obs_prop.size();

  // Offset proportions; observed counts are proportions times N.
  vector<Type> obs_n = obs_prop;
  vector<Type> hat_p = hat_prop;
  obs_n += offset;
  hat_p += offset;
  obs_n *= n_samp;

  // DM concentration. The count total N * (1 + offset * nbins) is the effective
  // sample size ddirmultinom() reads, here and in the AFSC linear-DM form.
  vector<Type> alphas = obs_n.sum() * hat_p * dm_theta;
  vector<Type> unweighted_alphas = obs_n.sum() * hat_p;

  if (osa_mode == 0) {
    switch (ll_type) {
    case -1:
      for (int bin = 0; bin < n_comp; bin++) {
        nll -= weight * n_samp * (obs_prop(bin) + offset) *
          log((hat_prop(bin) + offset) / (obs_prop(bin) + offset));
        unweighted_nll -= n_samp * (obs_prop(bin) + offset) *
          log((hat_prop(bin) + offset) / (obs_prop(bin) + offset));
      }
      break;

    case 0: {
      data_indicator<vector<Type>, Type> keep_ones(obs_n, true);
      nll -= weight * dmultinom_osa(obs_n, hat_p, keep_ones, 1, 1);
      unweighted_nll -= dmultinom_osa(obs_n, hat_p, keep_ones, 1, 1);
      break;
    }

    case 1:
      nll -= ddirmultinom(obs_n, alphas, true);
      unweighted_nll -= ddirmultinom(obs_n, unweighted_alphas, true);
      break;

    default:
      error("Invalid '%s'", type_name);
    }
  } else if (osa_start >= 0) {
    // OSA: unweighted and keep-gated, read from obsvec so oneStepPredict() can
    // residualize each bin. The beta-binomial has no elementary CDF, so the DM
    // gets no CDF term.
    vector<Type> osa_x = obsvec.segment(osa_start, n_comp);
    if (ll_type == 1) {
      nll -= ddirmultinom_osa(osa_x, alphas, keep.segment(osa_start, n_comp), 1, 1);
    } else {
      nll -= dmultinom_osa(osa_x, hat_p, keep.segment(osa_start, n_comp), 1, 1, osa_mode == 2);
    }
  }
}


/**
 * @brief Log density of a linkage prior on the natural scale.
 *
 * Families: 1 = normal(p1, p2); 2 = lognormal, log x ~ normal(p1 - b * p2^2 / 2, p2),
 * so exp(p1) is the mean when bias_adjust_proc b = 1 and the median when 0;
 * 3 = gamma(shape p1, rate p2); 4 = beta(p1, p2). Any other code returns 0.
 *
 * @param family Prior family code, from linkage_spec(priors = ).
 * @param value The parameter, on the scale it is stored.
 * @param value_is_log True if `value` is a log, so the prior reads exp(value).
 * @param p1,p2 The family's two parameters, as listed above.
 * @param bias_adjust_proc 1 to centre the lognormal on its mean, 0 on its median.
 * @return The log prior density (0 for an unknown family).
 */
template<class Type>
Type linkage_log_prior(int family, Type value, bool value_is_log,
                       Type p1, Type p2, Type bias_adjust_proc)
{
  Type x = value_is_log ? exp(value) : value;
  switch (family) {
  case 1:
    return dnorm(x, p1, p2, true);
  case 2: {
    // A stored log is used as is rather than re-logged.
    Type log_x = value_is_log ? value : log(x);
    return dnorm(log_x, p1 - bias_adjust_proc * square(p2) / 2.0, p2, true);
  }
  case 3:
    return dgamma(x, p1, Type(1.0) / p2, true);
  case 4:
    return dbeta(x, p1, p2, true);
  default:
    return Type(0.0);
  }
}

#endif
