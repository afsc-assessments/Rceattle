/**
 * @file comp_osa.hpp
 * @brief One-step-ahead (OSA) residual helpers for composition data.
 *
 * TMB's oneStepPredict() requires every observation entering the likelihood to
 * live in a flat `obsvec` with a `keep` indicator, and each likelihood term to
 * be multiplied by its keep element. Multinomial and Dirichlet-multinomial
 * composition data are multivariate and discrete, so oneStepPredict cannot use
 * them directly. Following Trijoulet et al. (2023) and WHAM
 * (src/age_comp_osa.hpp), a composition of A bins is decomposed into A-1
 * univariate conditional distributions -- binomials for the multinomial,
 * beta-binomials for the Dirichlet-multinomial -- each gated by its keep
 * element. The final bin is fixed by the sum-to-N constraint and contributes
 * nothing (its residual is not defined and is reported as NA on the R side).
 *
 * Key property: with keep == 1 the sum of the conditional log-densities equals
 * the joint multinomial / Dirichlet-multinomial log-density. Turning OSA on
 * therefore changes how the likelihood is *decomposed*, not its value.
 *
 * These functions are named with an `_osa` suffix and take a data_indicator so
 * they never collide with TMB's built-in dmultinom() or the plain
 * ddirmultinom() in helper_functions.hpp used during ordinary fitting.
 *
 * The file also holds the conditional CDFs oneStepPredict(method = "cdf") asks
 * for -- the continuous binomial for a composition bin, and the normal /
 * left-truncated normal for the aggregate series -- gated by
 * keep.cdf_lower / keep.cdf_upper. Under that method the residual is
 * Phi^-1(F(x)) and no conditional mean is ever formed, so nothing can leave the
 * support of a bounded observation. They are computed only when the caller asks
 * (`do_cdf`, driven by osa_mode == 2), so ordinary fitting and the Gaussian OSA
 * methods pay nothing for them.
 */

// Reorder categories so kept (keep > 0.5) entries come first, preserving the
// original order within the kept and dropped groups.
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

// Map the closed [0, 1] interval onto the open interval for numerical stability.
// Identical to TMB's convenience.hpp squeeze() (eps = machine epsilon), which is
// what WHAM uses in its OSA composition likelihood, so the two implementations
// agree to the last digit.
template<class Type>
Type osa_squeeze(Type u){
  Type eps = std::numeric_limits<double>::epsilon();
  return (Type(1.0) - eps) * (u - Type(0.5)) + Type(0.5);
}

/**
 * @brief Map a CDF onto the open interval before its logarithm is taken.
 *
 * Four machine epsilons wide, where osa_squeeze() uses one, and the extra width
 * is load-bearing. oneStepPredict recovers F from 1/(1 + exp(nlcdf.lower -
 * nlcdf.upper)) in double precision; at one epsilon that denominator is
 * 1 + 1.11e-16, which is exactly the round-half-to-even tie above 1, so it lands
 * on 1 or not depending on the last bits of an objective summed over thousands
 * of observations -- and a denominator of 1 gives F = 1 and an infinite
 * residual. Measured on BS2017SS: 6 of 4538 composition residuals came back
 * infinite at one epsilon and none at four.
 *
 * The cost is where |residual| saturates: 8.04 here against 8.13 at one epsilon
 * and 8.21 for a CDF read straight off. That ceiling is a property of double
 * precision, not of the model -- an observation that far into the tail of its
 * own conditional is beyond what any Q-Q plot resolves.
 *
 * The width is load-bearing a second way: this maps, it does not clamp, so a
 * caller that forms F by subtracting two CDFs relies on it to absorb the
 * cancellation. osa_trunc_norm_cdf_terms() below is the case -- its numerator
 * can round to about -2e-16 at an observation near zero, which this shift
 * carries back above 0 before the log. Narrowing it would return NaN there.
 */
template<class Type>
Type osa_squeeze_cdf(Type u){
  Type eps = Type(4.0) * std::numeric_limits<double>::epsilon();
  return (Type(1.0) - eps) * (u - Type(0.5)) + Type(0.5);
}

/**
 * @brief Conditional binomial CDF P(X <= x), X ~ Binomial(n, p), at fractional x.
 *
 * Written through the regularized incomplete beta, 1 - I_p(x + 1, n - x). That
 * equals the binomial CDF at whole counts and stays defined and smooth at the
 * fractional ones composition data carry: a bin count is proportion * effective
 * sample size, so it is only an integer by coincidence. Same expression as the
 * Trijoulet et al. (2023) reference implementation TMB ships
 * (contrib/OSA_multivariate_dists-main/distr.hpp, `dists::pbinom`), so the two
 * agree term for term.
 *
 * @param x Count in this bin.
 * @param n Count still unplaced across this bin and every later one.
 * @param p Probability of landing in this bin rather than a later one.
 */
template<class Type>
Type osa_pbinom(Type x, Type n, Type p){
  // n - x is what is left for the later bins. comp_offset adds a positive
  // amount to every bin, so it is strictly positive on every path that reaches
  // here; the floor is defensive, because pbeta at a zero shape returns NaN
  // rather than the limit, and a NaN times a zero cdf gate is still a NaN.
  Type nx = CppAD::CondExpLt(n - x, Type(1e-12), Type(1e-12), n - x);
  return Type(1.0) - pbeta(p, x + Type(1.0), nx);
}

/**
 * @brief One-step-ahead CDF terms for a normal observation.
 *
 * Every aggregate series in CEATTLE is normal on the scale `obsvec` stores it in
 * -- log for lognormal catch and index, natural scale for a "Normal" index, and
 * the whitened L^-1 scale for "MVN" -- so its conditional CDF is `pnorm` of the
 * standardized residual. Returns the term to ADD to the log-likelihood; it is 0
 * whenever the cdf gates are 0, which is every evaluation except the two
 * oneStepPredict(method = "cdf") makes per observation.
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
 * F(x) = [Phi(z) - Phi(-m)] / Phi(m) and its complement is 1 - F(x) =
 * Phi(-z) / Phi(m), because Phi(m) + Phi(-m) = 1. Both are written in that form
 * rather than as 1 - F, so the upper tail keeps its precision.
 *
 * This is exact wherever the truncation carries mass, which is what the Gaussian
 * OSA methods cannot see: the truncation enters the density only through
 * log Phi(m), a function of the prediction and not of the observation.
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
 * With `do_osa`, decomposes the A-bin multinomial into A-1 conditional binomials
 * (each gated by its keep element) so oneStepPredict can consume it; the sum
 * equals the joint multinomial log-density when all keeps are 1.
 *
 * @param x Observed counts (length A; may be non-integer, e.g. Neff * prop).
 * @param p Predicted proportions (length A; normalized internally).
 * @param keep Data indicator for these A bins (from keep.segment(...)).
 * @param do_osa 1 = conditional-binomial decomposition; 0 = plain multinomial.
 * @param give_log 1 = return the log-density.
 * @param do_cdf 1 = also add the keep.cdf_lower / keep.cdf_upper terms that
 *   oneStepPredict(method = "cdf") reads (osa_mode == 2). 0 skips the `pbeta`
 *   entirely, so ordinary fitting and the Gaussian OSA methods are unchanged.
 */
template <class Type>
Type dmultinom_osa(vector<Type> x, vector<Type> p,
                   data_indicator<vector<Type>, Type> keep,
                   int do_osa, int give_log, int do_cdf = 0)
{
  Type logres = 0;
  // Normalize predicted proportions. Guard the denominator so a length/age bin
  // with all-zero predicted mass (sum(p) == 0, e.g. a growth-transition row that
  // is exactly zero in the tails) yields p_x == 0 rather than 0/0 = NaN; the
  // conditional binomials below then read each squeezed-to-eps probability,
  // matching WHAM's behaviour for vanishing predicted proportions.
  vector<Type> p_x = p / (sum(p) + Type(1e-300));

  if(do_osa){
    vector<Type> k = keep;
    // The cdf gates are separate members of the data_indicator, so they have to
    // be reordered by hand alongside the keep vector. keep.segment() at the call
    // site already carries them (tmb_core.hpp documents that its segment IS
    // applied to cdf_lower/cdf_upper); the reorder here is the step that is not
    // done for us, and getting it wrong gates the wrong bin silently.
    //
    // osa_order() branches on keep, and under oneStepPredict keep is a PARAMETER,
    // so that branch is taped once at the initial values -- 1 for every bin that
    // will ever be residualized or conditioned on, 0 for the rest -- and the
    // per-observation keep changes afterwards do not re-sort. Within a
    // composition that initial pattern is 1 on every bin but the last, so the
    // order is the identity and each bin stays where the conditional
    // decomposition needs it. The reorder below is therefore a no-op on every
    // sequence this package generates; it is kept because it is what makes the
    // decomposition right for an arbitrary keep, and it is pinned by a test that
    // supplies one directly (test-likelihood-osa-cdf.R).
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
        // When one bin holds ~all remaining mass the separately-squeezed
        // numerator/denominator can push the ratio marginally above 1; clamp so
        // p2(1) = 1 - p2(0) stays non-negative (avoids dmultinom log(negative) =
        // NaN). A no-op for well-behaved compositions.
        Type eps_clamp = std::numeric_limits<double>::epsilon();
        p2(0) = CppAD::CondExpGt(p2(0), Type(1) - eps_clamp, Type(1) - eps_clamp, p2(0));
        p2(1) = Type(1) - p2(0);
        logres += k(i) * dmultinom(x2, p2, true);
        // Same conditional binomial the density above scores, read as a CDF.
        // Taken before nUnused is decremented, so the size matches the density.
        if(do_cdf) cdf = osa_pbinom(x(i), nUnused, p2(0));
        nUnused -= x(i);
        pUsed   += p_x(i);
        // Keep pUsed inside (0, 1) so 1 - pUsed never reaches 0 or goes
        // negative. WHAM does NOT squeeze pUsed itself; it relies on the
        // squeeze of 1 - pUsed above, which is also applied here. This extra
        // squeeze is therefore redundant protection and contracts pUsed toward
        // 0.5 by one machine epsilon per bin, so a composition of A bins sits
        // O(A * eps) from WHAM's value -- ~1e-15 relative, far below the
        // reference tolerances, but not bit-identical. Removing it would match
        // WHAM exactly at the cost of re-pinning the OSA references.
        pUsed    = osa_squeeze(pUsed);
      } else {
        logres += k(i) * Type(0);         // last bin: fixed by sum-to-N
        cdf = Type(1);                    // ... so it sits at the top of its own CDF
      }
      if(do_cdf){
        // Onto the open interval before taking logs: an exactly 0 or 1 CDF puts
        // a -Inf on the tape, and a -Inf times a zero gate is a NaN. See
        // osa_squeeze_cdf() for why the width matters and what it caps.
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
 * With `do_osa`, decomposes the A-bin Dirichlet-multinomial into A-1 conditional
 * beta-binomials (each gated by its keep element); the sum equals the joint D-M
 * log-density when all keeps are 1.
 *
 * There is no `do_cdf` counterpart here, so a Dirichlet-multinomial fleet cannot
 * be residualized with oneStepPredict(method = "cdf") and osa_residuals() sends
 * it to a Gaussian method instead. A beta-binomial has no elementary CDF, and
 * the reference implementation's route -- summing the pmf over 0..floor(x)
 * (contrib/OSA_multivariate_dists-main/distr.hpp, `pbetabinom`) -- is not open
 * to us on two counts: it is a step function of a fractional count, and it costs
 * O(x) beta functions per bin, which at CEATTLE sample sizes is thousands of
 * times the multinomial path.
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
