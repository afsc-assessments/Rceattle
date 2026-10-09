
/**
 * @file helper_functions.hpp
 * @brief Small differentiable (AD-safe) utilities used across the CEATTLE
 *   template: integer max (imax), positive-constraint penalty (posfun),
 *   Dirichlet-multinomial density (ddirmultinom), correlation transform
 *   (rho_trans), NA / finite tests, square, first_difference, a differentiable
 *   pairwise max (max2), and the separable year x age x cohort GMRF
 *   precision-matrix constructor (construct_Q).
 */

// Function for getting max of an IVECTOR and return an int
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

// Function to maintain > limit
template<class Type>
Type posfun(Type x, Type eps, Type &penalty) {
  Type denom = 2;
  denom -= x/eps;
  Type ans = CppAD::CondExpGe(x, eps, x, eps/denom);   // tune a new candidate x given eps
  penalty += CppAD::CondExpGe(x, eps, Type(0), 0.01 * (x - eps) * (x - eps));  // penalize via x vs eps
  return ans;
}


// Dirichlet multinomial (from WHAM)
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


// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


// Function for detecting Inf
template<class Type>
bool isFinite(Type x){
  return R_finite(asDouble(x));
}


//  Function to getting sqaure
template <class Type> Type square(Type x){return x*x;}


/**
 * @brief Concentrated (analytical) sd of a lognormal observation likelihood.
 *
 * Catch and index observations are fit as
 * `log(obs) ~ N(log(pred) - b * sigma^2 / 2, sigma)`, where `b` is the
 * bias-adjustment flag `bias_adjust_obs` (1 = predictions are mean-unbiased,
 * Rceattle's default; 0 = median-unbiased). Minimising that density over sigma
 * with `S` = mean squared log residual gives
 *
 *   sigma^2 = 2 * S / (sqrt(1 + b^2 * S) + 1)
 *
 * which is the Walters and Ludwig (1994) estimator `sigma = sqrt(S)` when
 * `b = 0`, and shrinks it as the bias term takes up part of the residual. The
 * form above is the rationalised version of `2 (sqrt(1 + b^2 S) - 1) / b^2`;
 * it needs no branch at `b = 0` and stays accurate for small `b`, where that
 * expression cancels to zero.
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


// Function to calculate the first difference
template <class Type>
vector<Type> first_difference(const vector<Type> &x)
{
  int length = x.size() - 1;

  if(length <= 0){
    std::cerr << "Error -- vector size too small in first_difference()" << std::endl;
    return vector<Type>(0); // explicit empty vector (was an ambiguous `return(0)`)
  }

  vector<Type> tmp(length);

  for(int i = 0; i < length; i++){
    tmp(i) = x(i+1) - x(i);
  }

  return tmp;
}


// Function to get max of two values (differentiable)
// https://groups.google.com/g/tmb-users/c/QAVmgj66OC0
template <class Type>
Type max2(Type x, Type y){
  Type ans = 0.5*(CppAD::abs(x-y)+x+y);
  return ans;
}


/**
 * @brief log(mean(exp(x))) over a non-empty vector, computed so it cannot overflow.
 *
 * The non-parametric selectivity forms rescale their log-scale coefficients
 * against log(mean(exp(coefficients))) once per year, which is how a curve whose
 * level is otherwise unidentified is held at mean selectivity 1. Summed directly,
 * exp() overflows to infinity once a coefficient passes about 709 -- the largest
 * argument it takes in double precision -- and the objective is NaN from there on.
 * Subtracting the largest coefficient first holds every exponent at or below 0, so
 * the sum is at most the number of bins whatever the coefficients are.
 *
 * The shift cancels out of both the result and its derivative.
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
 * Builds Q = (I - B^T) * Omega * (I - B), where B carries the partial
 * correlations by year (rho_y), age (rho_a), and cohort (rho_c), and Omega is
 * the diagonal (co)variance. `Var_Param` selects the conditional (0) or marginal
 * (1) variance parameterization of the GMRF.
 *
 * @param n_years Number of years.
 * @param n_ages Number of ages.
 * @param ay_Index [total_n x 2] age/year index locating each node.
 * @param rho_y Partial correlation across years.
 * @param rho_a Partial correlation across ages.
 * @param rho_c Partial correlation across cohorts.
 * @param log_sigma2 Log variance parameter governing the GMRF.
 * @param Var_Param 0 = conditional variance, 1 = marginal variance.
 * @return Sparse precision matrix [total_n x total_n], total_n = n_years * n_ages.
 */
template<class Type>
Eigen::SparseMatrix<Type> construct_Q(int n_years, // Integer of years
                                      int n_ages, // Integer of ages
                                      matrix<Type> ay_Index, // Index matrix to construct
                                      Type rho_y, // Partial correlation by years
                                      Type rho_a, // Partial correlation by ages
                                      Type rho_c, // Partial correlation by cohort
                                      Type log_sigma2, // Variance parameter governing GMRF
                                      int Var_Param // Parameterization of Variance ==0 (Conditional), == 1(Marginal)
) {

  // Dimension to construct matrices
  int total_n = n_years * n_ages;

  // Construct matrices for precision matrix
  Eigen::SparseMatrix<Type> B(total_n,total_n); // B matrix
  Eigen::SparseMatrix<Type> I(total_n,total_n); // Identity matrix
  I.setIdentity(); // Set I to identity matrix
  Eigen::SparseMatrix<Type> Omega(total_n,total_n); // Omega matrix (variances)
  Eigen::SparseMatrix<Type> Q_sparse(total_n, total_n); // Precision matrix

  for(int n = 0; n < total_n; n++) {

    // Define year and age objects
    Type age = ay_Index(n,0);
    Type year = ay_Index(n,1);

    // Constructing B matrix to determine where the correlation pars should go
    if(age > 1) {

      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year)
          B.coeffRef(n, n1) = rho_y;
      } // n1 loop

    } // end age > 1

    if(year > 1) {

      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age && ay_Index(n1, 1) == year - 1)
          B.coeffRef(n, n1) = rho_a;
      } // n1 loop

    } // if year > 1

    if(year > 1 && age > 1) {

      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year - 1)
          B.coeffRef(n,n1) = rho_c; // correlation by cohort
      } // n1 loop

    } // if both year and age > 1

  } // end n loop

  // Fill in Omega matrix here (variances)
  if(Var_Param == 0) { // Conditional variance

    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/exp(log_sigma2);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop

  } // end if conditional variance

  if(Var_Param == 1) { // Marginal Variance

    // Construct container objects
    matrix<Type> L(total_n, total_n); // L Matrix
    matrix<Type> tmp_I_B = I-B; // Temporary Matrix to store I-B
    L =  tmp_I_B.inverse(); // Invert to get L
    vector<Type> d(total_n); // Store variance calculations

    for(int n = 0; n < total_n; n++) {
      if(n == 0) {
        d(n) = exp(log_sigma2); // marginal variance parameter
      } else{

        Type cumvar = 0; // Cumulative Variance Container

        for(int n1 = 0; n1 < n; n1++) {
          cumvar += L(n, n1) * d(n1) * L(n, n1);
        } // n1 loop

        // Calculate diagonal values for omega
        d(n) = (exp(log_sigma2) - cumvar) / pow(L(n, n), 2);

      } // else loop
    } // n loop

    // Now fill in our diagonals for Omega
    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/d(i);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop

  } // end if marginal variance

  // Now, do calculations to construct (Q = (I - t(B)) %*% Omega %*% (I-B))
  Eigen::SparseMatrix<Type> B_transpose = B.transpose(); // transpose B matrix

  // Calculate Precision Matrix
  Q_sparse = (I - B_transpose) * Omega * (I-B);

  return(Q_sparse);

} // end construct_Q function

