
/**
 * @file helper_functions.hpp
 * @brief Small differentiable (AD-safe) utilities used across the CEATTLE
 *   template: numeric guards and transforms, the Dirichlet-multinomial density,
 *   and the separable year x age x cohort GMRF precision matrix (construct_Q).
 */

// Max of an IVECTOR, returned as an int.
template <class Type>
Type imax(const vector<Type> &x)
{
  int res = x[0];
  for(int i=0; i < x.size(); i++){
    int res2 = x[i];
    if(res < res2){
      res = res2;
    }
  }
  return res;
}

// Smooth positive constraint: returns x where x >= eps, otherwise a value in
// (0, eps), and adds a quadratic penalty for the shortfall.
template<class Type>
Type posfun(Type x, Type eps, Type &penalty) {
  Type denom = 2;
  denom -= x/eps;
  Type ans = CppAD::CondExpGe(x, eps, x, eps/denom);
  penalty += CppAD::CondExpGe(x, eps, Type(0), 0.01 * (x - eps) * (x - eps));
  return ans;
}


// Dirichlet-multinomial density (log when do_log = 1), from WHAM
// src/age_comp_osa.hpp (GPL-3).
template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> alpha, int do_log)
{
  int dim = obs.size();
  Type N = obs.sum();
  Type phi=sum(alpha);
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for(int a = 0; a < dim; a++) ll += -lgamma(obs(a) + 1.0) + lgamma(obs(a) + alpha(a)) - lgamma(alpha(a));
  if(do_log == 1) return ll;
  else return exp(ll);
}

// From TMB examples, modified name, constrains between -1 and 1
template <class Type>
Type rho_trans(Type x){
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}


template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


template<class Type>
bool isFinite(Type x){
  return R_finite(asDouble(x));
}


template <class Type> Type square(Type x){return x*x;}


/**
 * @brief Concentrated (analytical) sd of a lognormal observation likelihood.
 *
 * Catch and index observations are fit as `log(obs) ~ N(log(pred) - b * sigma^2 / 2,
 * sigma)`, with `b` the bias-adjustment flag `bias_adjust_obs` (1 = mean-unbiased
 * predictions, Rceattle's default; 0 = median-unbiased). Minimising that density over
 * sigma with `S` = mean squared log residual gives
 *
 *   sigma^2 = 2 * S / (sqrt(1 + b^2 * S) + 1)
 *
 * the Walters and Ludwig (1994) estimator `sigma = sqrt(S)` at `b = 0`, shrunk as the
 * bias term takes up part of the residual. It is the rationalised form of `2 (sqrt(1 +
 * b^2 S) - 1) / b^2`, which needs no branch at `b = 0` and stays accurate for small `b`,
 * where that expression cancels to zero.
 *
 * @param mean_sq_resid  S, the mean of squared `log(obs) - log(pred)`.
 * @param bias_adjust    b, the bias-adjustment flag applied by the density.
 * @return The sd that minimises the lognormal negative log-likelihood.
 */
template <class Type>
Type concentrated_lognormal_sd(Type mean_sq_resid, Type bias_adjust)
{
  return sqrt( Type(2.0) * mean_sq_resid /
               (sqrt(Type(1.0) + square(bias_adjust) * mean_sq_resid) + Type(1.0)) );
}


template <class Type>
vector<Type> first_difference(const vector<Type> &x)
{
  int length = x.size() - 1;

  if(length <= 0){
    std::cerr << "Error -- vector size too small in first_difference()" << std::endl;
    return vector<Type>(0); // no first difference exists below two elements
  }

  vector<Type> tmp(length);

  for(int i = 0; i < length; i++){
    tmp(i) = x(i+1) - x(i);
  }

  return tmp;
}


// Differentiable max of two values.
// https://groups.google.com/g/tmb-users/c/QAVmgj66OC0
template <class Type>
Type max2(Type x, Type y){
  Type ans = 0.5*(CppAD::abs(x-y)+x+y);
  return ans;
}


/**
 * @brief log(mean(exp(x))) over a non-empty vector, computed so it cannot overflow.
 *
 * The non-parametric selectivity forms rescale their log-scale coefficients against
 * log(mean(exp(coefficients))) once a year, holding a curve whose level is otherwise
 * unidentified at mean selectivity 1. Summed directly, exp() overflows to infinity once a
 * coefficient passes about 709, the largest argument it takes in double precision, and
 * the objective is NaN thereafter. Subtracting the largest coefficient first holds every
 * exponent at or below 0; the shift cancels out of both the result and its derivative.
 *
 * @param x non-empty vector of log-scale values.
 * @return log of the mean of exp(x).
 */
template <class Type>
Type log_mean_exp(const vector<Type> &x)
{
  Type mx = x(0);
  for(int i = 1; i < x.size(); i++) mx = max2(mx, x(i));

  Type total = 0;
  for(int i = 0; i < x.size(); i++) total += exp(x(i) - mx);

  return mx + log(total / Type(x.size()));
}


/**
 * @brief Construct a sparse GMRF precision matrix separable along year, age, and cohort.
 *
 * From Matthew Cheng's GMRF_WAA (https://github.com/chengmatt/GMRF_WAA; Cheng et
 * al. 2023, Fisheries Research 266: 106755). Builds Q = (I - B^T) * Omega * (I - B), where B carries the partial correlation to each
 * neighbouring node and Omega is the diagonal precision (1 / variance). `Var_Param`
 * selects the conditional (0) or marginal (1) variance parameterization of the GMRF.
 *
 * The sole caller reads these names the other way round: ceattle.cpp's 3D-AR1 selectivity
 * passes its (bin, year, cohort) correlations in and fills `ay_Index` as (bin, year), so
 * the two inversions cancel. Read the note at that call before renaming anything here.
 *
 * @param n_years Number of years.
 * @param n_ages Number of ages.
 * @param ay_Index [total_n x 2] age/year index locating each node.
 * @param rho_y Partial correlation to the node one back in `ay_Index` column 0.
 * @param rho_a Partial correlation to the node one back in column 1.
 * @param rho_c Partial correlation to the node one back in both -- the diagonal,
 *   which is the cohort when the columns are filled (age, year).
 * @param log_sigma2 Log variance parameter governing the GMRF.
 * @param Var_Param 0 = conditional variance, 1 = marginal variance.
 * @return Sparse precision matrix [total_n x total_n], total_n = n_years * n_ages.
 */
template<class Type>
Eigen::SparseMatrix<Type> construct_Q(int n_years,
                                      int n_ages,
                                      matrix<Type> ay_Index,
                                      Type rho_y,
                                      Type rho_a,
                                      Type rho_c,
                                      Type log_sigma2,
                                      int Var_Param
) {

  int total_n = n_years * n_ages;

  Eigen::SparseMatrix<Type> B(total_n,total_n);
  Eigen::SparseMatrix<Type> I(total_n,total_n);
  I.setIdentity();
  Eigen::SparseMatrix<Type> Omega(total_n,total_n);
  Eigen::SparseMatrix<Type> Q_sparse(total_n, total_n);

  for(int n = 0; n < total_n; n++) {

    Type age = ay_Index(n,0);
    Type year = ay_Index(n,1);

    // B holds the partial correlation from node n to each of its neighbours.
    if(age > 1) {

      // One step back in ay_Index column 0, same column 1.
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year)
          B.coeffRef(n, n1) = rho_y;
      } // n1 loop

    } // end age > 1

    if(year > 1) {

      // One step back in column 1, same column 0.
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age && ay_Index(n1, 1) == year - 1)
          B.coeffRef(n, n1) = rho_a;
      } // n1 loop

    } // if year > 1

    if(year > 1 && age > 1) {

      // One step back in both columns: the diagonal (cohort) neighbour.
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year - 1)
          B.coeffRef(n,n1) = rho_c;
      } // n1 loop

    } // if both year and age > 1

  } // end n loop

  // Omega is diagonal and holds precisions, 1 / variance.
  if(Var_Param == 0) { // Conditional variance

    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/exp(log_sigma2);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop

  } // end if conditional variance

  if(Var_Param == 1) { // Marginal Variance

    matrix<Type> L(total_n, total_n);
    matrix<Type> tmp_I_B = I-B;
    L =  tmp_I_B.inverse();
    vector<Type> d(total_n);   // conditional variance per node

    for(int n = 0; n < total_n; n++) {
      if(n == 0) {
        d(n) = exp(log_sigma2); // marginal variance parameter
      } else{

        Type cumvar = 0;

        for(int n1 = 0; n1 < n; n1++) {
          cumvar += L(n, n1) * d(n1) * L(n, n1);
        } // n1 loop

        // Set d(n) so node n's MARGINAL variance comes out at exp(log_sigma2).
        d(n) = (exp(log_sigma2) - cumvar) / pow(L(n, n), 2);

      } // else loop
    } // n loop

    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/d(i);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop

  } // end if marginal variance

  Eigen::SparseMatrix<Type> B_transpose = B.transpose();

  Q_sparse = (I - B_transpose) * Omega * (I-B);

  return(Q_sparse);

} // end construct_Q function

