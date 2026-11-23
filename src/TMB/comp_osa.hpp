/**
 * @file comp_osa.hpp
 * @brief One-step-ahead (OSA) residual helpers for composition data.
 *
 * TMB's oneStepPredict() needs every observation in a flat `obsvec`, its likelihood
 * term multiplied by its own `keep` element, which multivariate discrete
 * composition data cannot supply directly. Following Trijoulet et al. (2023), their
 * reference code (distr.hpp in github.com/vtrijoulet/OSA_multivariate_dists, from
 * which dmultinom_osa() and ddirmultinom_osa() take their names and CDF terms) and
 * WHAM (src/age_comp_osa.hpp), a composition of A bins is decomposed into A-1
 * univariate conditionals -- binomials for the multinomial, beta-binomials for the
 * Dirichlet-multinomial -- each gated by its keep element. The last bin is fixed by
 * the sum-to-N constraint, contributes nothing, and is reported as NA on the R side.
 * With every keep at 1 the conditional log-densities sum to the joint log-density,
 * so OSA changes how the likelihood is *decomposed*, not its value. The `_osa`
 * suffix and the data_indicator argument keep these clear of TMB's dmultinom() and
 * of the ddirmultinom() in helper_functions.hpp that ordinary fitting uses.
 *
 * The file also holds the conditional CDFs oneStepPredict(method = "cdf") reads --
 * continuous binomial for a composition bin, normal / left-truncated normal for an
 * aggregate series. That residual is Phi^-1(F(x)), formed without a conditional
 * mean, so it cannot leave the support of a bounded observation.
 *
 * Shared contract of the CDF terms below -- the two `_cdf_terms` helpers and
 * dmultinom_osa()'s `do_cdf` path. Each contributes to the log-likelihood
 * `keep.cdf_lower` times log F plus `keep.cdf_upper` times log of its complement,
 * which is 0 whenever the cdf gates are 0: every evaluation except the two
 * oneStepPredict(method = "cdf") makes per observation. The terms are added only when
 * the caller asks (`do_cdf`, from osa_mode == 2), so ordinary fitting and the Gaussian
 * OSA methods are unchanged and pay nothing for them. Every F is mapped onto the open
 * interval by osa_squeeze_cdf() first, because an exactly 0 or 1 CDF puts a -Inf on
 * the tape and a -Inf times a zero gate is a NaN.
 */

// Reorder categories so kept (keep > 0.5) entries come first, preserving the
// original order within the kept and dropped groups. Anders Nielsen's order(),
// via WHAM src/age_comp_osa.hpp.
template<class Type>
vector<int> osa_order(vector<Type> k){
  int n = k.size();
  vector<int> ord(n);
  ord.setZero();
  int at = -1;
  for(int i = 0; i < n; ++i){
    if(k(i) > 0.5){ ord(++at) = i; }
  }
  at = n;
  for(int i = n - 1; i >= 0; --i){
    if(k(i) < 0.5){ ord(--at) = i; }
  }
  return ord;
}

// Map the closed [0, 1] interval onto the open interval. Identical to TMB's
// convenience.hpp squeeze() (eps = machine epsilon), which WHAM's OSA composition
// likelihood uses, so the two agree to the last digit.
template<class Type>
Type osa_squeeze(Type u){
  Type eps = std::numeric_limits<double>::epsilon();
  return (Type(1.0) - eps) * (u - Type(0.5)) + Type(0.5);
}

/**
 * @brief Map a CDF onto the open interval before its logarithm is taken.
 *
 * Four machine epsilons wide where osa_squeeze() uses one: oneStepPredict recovers F
 * from 1/(1 + exp(nlcdf.lower - nlcdf.upper)) in double precision, and at one epsilon
 * that denominator sits on the round-half-to-even tie above 1 -- where it rounds to 1
 * so does F, and the residual is infinite. Measured on BS2017SS: 6 of 4538
 * composition residuals infinite at one epsilon, none at four; the cost is where
 * |residual| saturates -- 8.04 against 8.13 -- a tail no Q-Q plot resolves.
 *
 * It maps rather than clamps, so a caller that forms F by subtracting two CDFs relies
 * on it to absorb the cancellation: osa_trunc_norm_cdf_terms()'s numerator can round
 * to about -2e-16 at an observation near zero, which this shift carries back above 0
 * before the log. Narrowing it would return NaN there.
 */
template<class Type>
Type osa_squeeze_cdf(Type u){
  Type eps = Type(4.0) * std::numeric_limits<double>::epsilon();
  return (Type(1.0) - eps) * (u - Type(0.5)) + Type(0.5);
}

/**
 * @brief Conditional binomial CDF P(X <= x), X ~ Binomial(n, p), at fractional x.
 *
 * Written through the regularized incomplete beta, 1 - I_p(x + 1, n - x), which equals
 * the binomial CDF at whole counts and stays smooth at the fractional ones composition
 * data carry: a bin count is proportion * effective sample size, an integer only by
 * coincidence. Same expression as the Trijoulet et al. (2023) reference implementation
 * (distr.hpp, `dists::pbinom`, which TMB:::install.contrib() installs).
 *
 * @param x Count in this bin.
 * @param n Count still unplaced across this bin and every later one.
 * @param p Probability of landing in this bin rather than a later one.
 */
template<class Type>
Type osa_pbinom(Type x, Type n, Type p){
  // n - x is the count left for the later bins. comp_offset (default 1e-5) is added to
  // every bin, so it is positive unless comp_offset is 0, where pbeta at a zero shape
  // returns NaN rather than the limit -- and a NaN times a zero gate is still a NaN.
  Type nx = CppAD::CondExpLt(n - x, Type(1e-12), Type(1e-12), n - x);
  return Type(1.0) - pbeta(p, x + Type(1.0), nx);
}

/**
 * @brief One-step-ahead CDF terms for a normal observation.
 *
 * Every aggregate series in CEATTLE is normal on the scale `obsvec` stores it in --
 * log for lognormal catch and index, natural scale for a "Normal" index, the whitened
 * L^-1 scale for "MVN" -- so its conditional CDF is `pnorm` of the standardized
 * residual.
 *
 * @param x Observation, on the scale it is stored in.
 * @param mu Predicted value on that same scale.
 * @param sd Standard deviation on that same scale.
 * @param cdf_lower,cdf_upper The observation's `keep.cdf_lower`/`keep.cdf_upper`.
 */
template<class Type>
Type osa_norm_cdf_terms(Type x, Type mu, Type sd, Type cdf_lower, Type cdf_upper){
  Type F = osa_squeeze_cdf(pnorm((x - mu) / sd));
  return cdf_lower * log(F) + cdf_upper * log(Type(1.0) - F);
}

/**
 * @brief One-step-ahead CDF terms for a normal left-truncated at zero.
 *
 * With z = (x - mu)/sd and m = mu/sd, the truncated CDF is
 * F(x) = [Phi(z) - Phi(-m)] / Phi(m) and its complement is 1 - F(x) = Phi(-z) / Phi(m),
 * because Phi(m) + Phi(-m) = 1. Both are written in that form rather than as 1 - F, so
 * the upper tail keeps its precision. This is exact wherever the truncation carries
 * mass, which is what the Gaussian OSA methods cannot see: the truncation enters the
 * density only through log Phi(m), a function of the prediction and not of the
 * observation.
 *
 * @param x Observation (natural scale, > 0).
 * @param mu Predicted index (natural scale).
 * @param sd Observation standard deviation (absolute, natural scale).
 * @param cdf_lower,cdf_upper The observation's `keep.cdf_lower`/`keep.cdf_upper`.
 */
template<class Type>
Type osa_trunc_norm_cdf_terms(Type x, Type mu, Type sd,
                              Type cdf_lower, Type cdf_upper){
  Type m  = mu / sd;
  Type Zm = pnorm(m);
  Type F  = osa_squeeze_cdf((pnorm((x - mu) / sd) - pnorm(-m)) / Zm);
  Type G  = osa_squeeze_cdf(pnorm(-(x - mu) / sd) / Zm);
  return cdf_lower * log(F) + cdf_upper * log(G);
}

/**
 * @brief Keep-aware multinomial log-density via conditional binomials.
 *
 * @param x Observed counts (length A; may be non-integer, e.g. Neff * prop).
 * @param p Predicted proportions (length A; normalized internally).
 * @param keep Data indicator for these A bins (from keep.segment(...)).
 * @param do_osa 1 = conditional-binomial decomposition; 0 = plain multinomial.
 * @param give_log 1 = return the log-density.
 * @param do_cdf 1 = also add the CDF terms (osa_mode == 2); 0 skips the `pbeta`.
 */
template <class Type>
Type dmultinom_osa(vector<Type> x, vector<Type> p,
                   data_indicator<vector<Type>, Type> keep,
                   int do_osa, int give_log, int do_cdf = 0)
{
  Type logres = 0;
  // Normalize predicted proportions, guarding the denominator so a bin vector with no
  // predicted mass (sum(p) == 0, e.g. a growth-transition row exactly zero in the
  // tails) gives p_x == 0 rather than 0/0 = NaN, matching WHAM for vanishing p.
  vector<Type> p_x = p / (sum(p) + Type(1e-300));

  if(do_osa){
    vector<Type> k = keep;
    // The cdf gates are separate data_indicator members, so they are reordered by
    // hand alongside keep; keep.segment() at the call site already carries them
    // (tmb_core.hpp), and gating the wrong bin is silent.
    //
    // osa_order() branches on keep, a PARAMETER under oneStepPredict, so the sort is
    // taped once at the initial values -- 1 for every bin that will ever be
    // residualized or conditioned on -- and later per-observation keep changes do not
    // re-sort. Within a composition that pattern is 1 on every bin but the last, so
    // the order is the identity and this reorder is a no-op on every sequence the
    // package generates. Kept because it is what makes the decomposition right for an
    // arbitrary keep (test-likelihood-osa-cdf.R supplies one directly).
    vector<Type> cl = keep.cdf_lower;
    vector<Type> cu = keep.cdf_upper;
    vector<int> o = osa_order(k);
    x = x(o); p_x = p_x(o); k = k(o);     // kept bins first
    if(do_cdf){ cl = cl(o); cu = cu(o); }
    Type nUnused = asDouble(sum(x));
    Type pUsed = 0;
    Type cdf = 1;
    for(int i = 0; i < x.size(); ++i){
      if(i != (x.size() - 1)){
        // Bin i versus all remaining bins: a binomial written as a 2-category
        // multinomial. Conditional success probability p_i / (1 - pUsed).
        vector<Type> x2(2), p2(2);
        x2(0) = x(i);
        x2(1) = nUnused - x(i);
        p2(0) = osa_squeeze(p_x(i)) / osa_squeeze(Type(1) - pUsed);
        p2(0) = osa_squeeze(p2(0));
        // When one bin holds ~all remaining mass the separately-squeezed numerator and
        // denominator can push the ratio marginally above 1; clamp so p2(1) = 1 - p2(0)
        // stays non-negative: dmultinom of a negative probability is log(negative) = NaN.
        Type eps_clamp = std::numeric_limits<double>::epsilon();
        p2(0) = CppAD::CondExpGt(p2(0), Type(1) - eps_clamp, Type(1) - eps_clamp, p2(0));
        p2(1) = Type(1) - p2(0);
        logres += k(i) * dmultinom(x2, p2, true);
        // Same conditional binomial the density above scores, read as a CDF.
        // Taken before nUnused is decremented, so the size matches the density.
        if(do_cdf) cdf = osa_pbinom(x(i), nUnused, p2(0));
        nUnused -= x(i);
        pUsed   += p_x(i);
        // Keep pUsed inside (0, 1) so 1 - pUsed never reaches 0 or goes negative. WHAM
        // squeezes only 1 - pUsed, so this extra squeeze is redundant protection that
        // contracts pUsed toward 0.5 by one machine epsilon per bin: a composition of A
        // bins sits O(A * eps), ~1e-15 relative, from WHAM's value. Removing it would
        // match WHAM exactly at the cost of re-pinning the OSA references.
        pUsed    = osa_squeeze(pUsed);
      } else {
        logres += k(i) * Type(0);         // last bin: fixed by sum-to-N
        cdf = Type(1);                    // ... so it sits at the top of its own CDF
      }
      if(do_cdf){
        // Onto the open interval before the logs; the width matters (osa_squeeze_cdf()).
        cdf = osa_squeeze_cdf(cdf);
        logres += cl(i) * log(cdf);
        logres += cu(i) * log(Type(1) - cdf);
      }
    }
  } else {
    logres = dmultinom(x, p_x, true);
  }
  return give_log ? logres : exp(logres);
}

/**
 * @brief Keep-aware Dirichlet-multinomial log-density via conditional beta-binomials.
 *
 * No `do_cdf` counterpart: a beta-binomial has no elementary CDF, and the reference
 * route -- summing the pmf over 0..floor(x)
 * (contrib/OSA_multivariate_dists-main/distr.hpp, `pbetabinom`) -- is a step function
 * of a fractional count costing O(x) beta functions per bin, thousands of times the
 * multinomial path at CEATTLE sample sizes. So a Dirichlet-multinomial fleet cannot be
 * residualized with oneStepPredict(method = "cdf"); osa_residuals() sends it to a
 * Gaussian method instead.
 *
 * @param obs Observed counts (length A).
 * @param alpha Dirichlet-multinomial concentration parameters (length A).
 * @param keep Data indicator for these A bins.
 * @param do_osa 1 = conditional beta-binomial decomposition; 0 = plain D-M.
 * @param give_log 1 = return the log-density.
 */
template<class Type>
Type ddirmultinom_osa(vector<Type> obs, vector<Type> alpha,
                      data_indicator<vector<Type>, Type> keep,
                      int do_osa, int give_log)
{
  Type ll = 0.0;
  if(do_osa){
    vector<Type> k = keep;
    vector<int> o = osa_order(k);
    obs = obs(o); alpha = alpha(o); k = k(o);   // kept bins first
    int dim = obs.size();
    Type alp_sum = alpha.sum();
    Type obs_sum = asDouble(obs.sum());
    for(int a = 0; a < dim - 1; ++a){
      // Bin a versus all remaining bins: a beta-binomial written as a
      // 2-category Dirichlet-multinomial (alpha_a vs sum of remaining alphas).
      obs_sum -= obs(a);
      alp_sum -= alpha(a);
      vector<Type> obs_a(2), alphas_a(2);
      obs_a(0)    = obs(a);   obs_a(1)    = obs_sum;
      alphas_a(0) = alpha(a); alphas_a(1) = alp_sum;
      ll += k(a) * ddirmultinom(obs_a, alphas_a, 1);
    }
  } else {
    ll = ddirmultinom(obs, alpha, 1);
  }
  return give_log ? ll : exp(ll);
}
