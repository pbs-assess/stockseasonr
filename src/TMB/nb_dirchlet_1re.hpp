/// @file nb_dirchlet_1re.hpp

#ifndef nb_dirchlet_1re_hpp
#define nb_dirchlet_1re_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// Negative log-likelihood of the normal distribution.
template<class Type>
Type nb_dirchlet_1re(objective_function<Type>* obj) {
  // DATA
  // Abundance 
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  DATA_IVECTOR(factor1k_i);
  DATA_INTEGER(nk1);
  DATA_MATRIX(X1_pred_ij);
  // vector of higher level aggregates used to generate predictions; length
  // is equal to the number of predictions made
  DATA_IVECTOR(pred_factor2k_h);
  DATA_IVECTOR(pred_factor2k_levels);
  // Composition 
  DATA_MATRIX(y2_ig);   // matrix of observed distribuons for g 
  DATA_MATRIX(X2_ij);       // model matrix for fixed effects
  DATA_IVECTOR(factor2k_i); // vector of random factor levels
  DATA_INTEGER(nk2);    // number of factor levels
  DATA_MATRIX(X2_pred_ij);    // prediction matrix for compositon
  DATA_INTEGER(random_walk);
  
  // PARAMETERS 
  // Abundance 
  PARAMETER_VECTOR(b1_j);
  PARAMETER(log_phi);
  PARAMETER_VECTOR(z1_k);
  PARAMETER(log_sigma_zk1);
  // Composition
  PARAMETER_MATRIX(b2_jg);   // matrix of fixed int. (rows = fixed cov, cols = g)
  PARAMETER_VECTOR(z2_k);    // vector of random int.
  PARAMETER(log_sigma_zk2);  // among random int SD
  
  
  // DEFINE INTERMEDIATES
  // Abundance 
  int n1 = y1_i.size();
  vector<Type> linear_predictor1_i(n1);
  // int n_preds1 = X1_pred_ij.rows();   // number of abund predictions (finer)
  // Composition 
  int n2 = y2_ig.rows();              // number of observations
  int n_cat = y2_ig.cols();           // number of categories
  // int n_preds2 = X2_pred_ij.rows();   // number of comp predictions (coarser)
  matrix<Type> total_eff(n2, n_cat);
  
  
  // LINEAR PREDICTORS
  // Abundance
  linear_predictor1_i = X1_ij * b1_j;
  
  // Composition linear
  matrix<Type> fx_eff = X2_ij * b2_jg;
  
  for (int i = 0; i < n2; ++i) {
    for(int k = 0; k < n_cat; k++) {
      total_eff(i, k) = fx_eff(i, k) + z2_k(factor2k_i(i));
    }
  }
  
  matrix<Type> gamma = exp(total_eff.array()); // add random effect
  vector<Type> n_plus = y2_ig.rowwise().sum(); // row sum of response
  vector<Type> gamma_plus = gamma.rowwise().sum(); // row sum of gamma
  
  
  // LOG-LIKELIHOOD
  Type jll = 0; // initialize joint log-likelihood
  // Composition ll
  for(int i = 0; i <= (n2 - 1); i++){
    jll = jll + lgamma((n_plus(i) + 1));
    jll = jll + lgamma(gamma_plus(i));
    jll = jll - lgamma((n_plus(i) + gamma_plus(i)));
    for(int k = 0; k <= (n_cat - 1); k++){
      jll += lgamma((y2_ig(i, k) + gamma(i, k)));
      jll -= lgamma(gamma(i, k));
      jll -= lgamma((y2_ig(i, k) + 1));
    }
  }
  
  Type jnll;
  jnll = -jll;
  
  // Abundance nll
  Type s1, s2;
  for(int i = 0; i < n1; i++){
    s1 = linear_predictor1_i(i) + z1_k(factor1k_i(i)); //mu
    s2 = 2.0 * (s1 - log_phi); //scale
    jnll -= dnbinom_robust(y1_i(i), s1, s2, true);
  }
  
  // Probability of random composition coefficients
  if (random_walk == 1) {
    for (int k = 0; k < nk1; k++) {
      if (k == 0) {
        jnll -= dnorm(z1_k(k), Type(0.0), exp(log_sigma_zk1), true);  
      }
      if (k > 0) {
        jnll -= dnorm(z1_k(k), z1_k(k - 1), exp(log_sigma_zk1), true);
      }
    }
    for (int k = 0; k < nk2; k++) {
      if (k == 0) {
        jnll -= dnorm(z2_k(k), Type(0.0), exp(log_sigma_zk2), true);  
      }
      if (k > 0) {
        jnll -= dnorm(z2_k(k), z2_k(k - 1), exp(log_sigma_zk2), true);
      }
    }
  } else {
    for (int k = 0; k < nk1; k++){
      jnll -= dnorm(z1_k(k), Type(0.0), exp(log_sigma_zk1), true);
    }
    for (int k = 0; k < nk2; k++) {
      jnll -= dnorm(z2_k(k), Type(0.0), exp(log_sigma_zk2), true);
    }
  }
  
  
  // FIXED EFFECTS PREDICTIONS
  // Abundance 
  vector<Type> log_pred_abund(X1_pred_ij.rows());
  log_pred_abund = X1_pred_ij * b1_j;
  matrix<Type> pred_abund = exp(log_pred_abund.array());
  
  // REPORT(pred_abund);
  ADREPORT(log_pred_abund);
  
  // Calculate predicted abundance based on higher level groupings
  int n_preds = pred_factor2k_h.size();
  int n_pred_levels = pred_factor2k_levels.size();
  vector<Type> agg_pred_abund(n_pred_levels);
  vector<Type> log_agg_pred_abund(n_pred_levels);
  
  for (int h = 0; h < n_preds; h++) {
    for (int m = 0; m < n_pred_levels; m++) {
      if (pred_factor2k_h(h) == pred_factor2k_levels(m)) {
        agg_pred_abund(m) += pred_abund(h);
        log_agg_pred_abund(m) = log(agg_pred_abund(m));
      }
    }
  }
  
  // ADREPORT(agg_pred_abund);
  ADREPORT(log_agg_pred_abund);
  
  
  // Composition 
  matrix<Type> pred_eff(n_pred_levels, n_cat);    //pred effects on log scale
  matrix<Type> pred_gamma(n_pred_levels, n_cat);  //transformed pred effects 
  vector<Type> pred_gamma_plus(n_pred_levels);        
  vector<Type> pred_theta(n_pred_levels); 
  matrix<Type> pred_pi(n_pred_levels, n_cat);      // predicted counts in real 
  vector<Type> pred_n_plus(n_pred_levels); 
  matrix<Type> pred_pi_prop(n_pred_levels, n_cat); // predicted counts as ppn.
  matrix<Type> inv_logit_pred_pi_prop(n_pred_levels, n_cat); // pred. ppn link
  
  pred_eff = X2_pred_ij * b2_jg; 
  pred_gamma = exp(pred_eff.array());
  pred_gamma_plus = pred_gamma.rowwise().sum();
  pred_theta = 1 / (pred_gamma_plus + 1);
  for(int m = 0; m < n_pred_levels; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_pi(m, k) = pred_gamma(m, k) / pred_theta(m);
    }
  }
  pred_n_plus = pred_pi.rowwise().sum();
  for(int m = 0; m < n_pred_levels; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_pi_prop(m, k) = pred_pi(m, k) / pred_n_plus(m);
      inv_logit_pred_pi_prop(m, k) = invlogit(pred_pi_prop(m, k));
    }
  }
  
  // REPORT(pred_pi_prop);
  // ADREPORT(pred_pi_prop);
  ADREPORT(inv_logit_pred_pi_prop);
  
  
  // Combined predictions
  matrix<Type> pred_abund_mg(n_pred_levels, n_cat);
  
  for (int m = 0; m < n_pred_levels; ++m) {
    for (int k = 0; k < n_cat; ++k) {
      pred_abund_mg(m, k) = agg_pred_abund(m) * pred_pi_prop(m, k);
    }
  }
  matrix<Type> log_pred_abund_mg = log(pred_abund_mg.array());
  
  // REPORT(pred_abund_mg);
  // ADREPORT(pred_abund_mg);
  ADREPORT(log_pred_abund_mg);
  
  return jnll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
