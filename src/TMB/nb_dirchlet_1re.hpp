/// @file nb_dirchlet_1re.hpp

#ifndef nb_dirchlet_1re_hpp
#define nb_dirchlet_1re_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// Negative log-likelihood of the normal distribution.
template<class Type>
Type nb_dirchlet_1re(objective_function<Type>* obj) {
  
  // DATA
  // abundance data
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  DATA_IVECTOR(factor1k_i);
  DATA_INTEGER(nk1);
  DATA_MATRIX(X1_pred_ij);
  // vector of higher level regions used to generate predictions; length
  // is equal to the number of predictions made
  DATA_IVECTOR(pred_factor2k_h);
  DATA_IVECTOR(pred_factor2k_levels);
  // composition data
  DATA_MATRIX(y2_ig);   // matrix of observed distributions for g 
  DATA_MATRIX(X2_ij);       // model matrix for fixed effects
  DATA_IVECTOR(factor2k_i); // vector of random factor levels
  DATA_INTEGER(nk2);    // number of factor levels
  DATA_MATRIX(X2_pred_ij);    // prediction matrix for compositon
  // conditionals
  DATA_INTEGER(abundance_component);
  DATA_INTEGER(random_walk);

  // PARAMETERS
  //abundance parameters
  PARAMETER_VECTOR(b1_j);
  PARAMETER(log_phi);
  PARAMETER_VECTOR(z1_k);
  PARAMETER(log_sigma_zk1);
  // composition parameters
  PARAMETER_MATRIX(b2_jg);   // matrix of fixed int. (rows = fixed cov, cols = g)
  PARAMETER_VECTOR(z2_k);    // vector of random int.
  PARAMETER(log_sigma_zk2);  // among random int SD

  Type jnll = 0.0;
  
  // CALCULATIONS
  if (abundance_component) {
    // intermediate storage
    int n1 = y1_i.size();
    vector<Type> linear_predictor1_i(n1);
    
    // linear predictors
    linear_predictor1_i = X1_ij * b1_j;
    
    // fixed-effects nll
    Type s1, s2;
    for(int i = 0; i < n1; i++){
      s1 = linear_predictor1_i(i) + z1_k(factor1k_i(i)); //mu
      s2 = 2.0 * (s1 - log_phi); //scale
      jnll -= dnbinom_robust(y1_i(i), s1, s2, true);
    }
  
    // random-effects nll
    if (random_walk) {
      for (int k = 0; k < nk1; k++) {
        if (k == 0) {
          jnll -= dnorm(z1_k(k), Type(0.0), exp(log_sigma_zk1), true);  
        }
        if (k > 0) {
          jnll -= dnorm(z1_k(k), z1_k(k - 1), exp(log_sigma_zk1), true);
        }
      }
    } else {
      for (int k = 0; k < nk1; k++){
        jnll -= dnorm(z1_k(k), Type(0.0), exp(log_sigma_zk1), true);
      }
    }    
  }


  // COMPOSITION COMPONENT
  // intermediate storage
  int n2 = y2_ig.rows();              // number of observations
  int n_cat = y2_ig.cols();           // number of categories
  matrix<Type> total_eff(n2, n_cat);
  
  // linear predictor
  matrix<Type> fx_eff = X2_ij * b2_jg;
  
  for (int i = 0; i < n2; ++i) {
    for(int k = 0; k < n_cat; k++) {
      total_eff(i, k) = fx_eff(i, k) + z2_k(factor2k_i(i));
    }
  }
  
  matrix<Type> gamma = exp(total_eff.array()); // add random effect
  vector<Type> n_plus = y2_ig.rowwise().sum(); // row sum of response
  vector<Type> gamma_plus = gamma.rowwise().sum(); // row sum of gamma
  
  // fixed-effects nll
  Type jll = 0; // initialize joint log-likelihood
  // positive because the mental math hurt my head...
  for (int i = 0; i <= (n2 - 1); i++){
    jll = jll + lgamma((n_plus(i) + 1));
    jll = jll + lgamma(gamma_plus(i));
    jll = jll - lgamma((n_plus(i) + gamma_plus(i)));
    for(int k = 0; k <= (n_cat - 1); k++){
      jll += lgamma((y2_ig(i, k) + gamma(i, k)));
      jll -= lgamma(gamma(i, k));
      jll -= lgamma((y2_ig(i, k) + 1));
    }
  }
  // now it's negative!
  jnll -= jll;
  
  // random-effects nll
  if (random_walk) {
    for (int k = 0; k < nk2; k++) {
      if (k == 0) {
        jnll -= dnorm(z2_k(k), Type(0.0), exp(log_sigma_zk2), true);  
      }
      if (k > 0) {
        jnll -= dnorm(z2_k(k), z2_k(k - 1), exp(log_sigma_zk2), true);
      }
    }
  } else {
    for (int k = 0; k < nk2; k++) {
      jnll -= dnorm(z2_k(k), Type(0.0), exp(log_sigma_zk2), true);
    }
  }
  
  
  // FIXED-EFFECTS PREDICTIONS
  // Composition
  int n_pred_levels = X2_pred_ij.rows();     
  matrix<Type> pred_eff(n_pred_levels, n_cat);    //pred effects on log scale
  matrix<Type> pred_gamma(n_pred_levels, n_cat);  //transformed pred effects 
  vector<Type> pred_gamma_plus(n_pred_levels);        
  vector<Type> pred_theta(n_pred_levels); 
  matrix<Type> pred_pi(n_pred_levels, n_cat);      // predicted counts in real 
  vector<Type> pred_n_plus(n_pred_levels); 
  matrix<Type> pred_pi_prop(n_pred_levels, n_cat); // predicted counts as ppn.
  matrix<Type> logit_pred_pi_prop(n_pred_levels, n_cat); // pred. ppn link 
  
  pred_eff = X2_pred_ij * b2_jg; 
  pred_gamma = exp(pred_eff.array());
  pred_gamma_plus = pred_gamma.rowwise().sum();
  pred_theta = 1 / (pred_gamma_plus + 1);
  for (int m = 0; m < n_pred_levels; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_pi(m, k) = pred_gamma(m, k) / pred_theta(m);
    }
  }
  pred_n_plus = pred_pi.rowwise().sum();
  for (int m = 0; m < n_pred_levels; m++) {
    for (int k = 0; k < n_cat; k++) {
      pred_pi_prop(m, k) = pred_pi(m, k) / pred_n_plus(m);
      logit_pred_pi_prop(m, k) = logit(pred_pi_prop(m, k));
    }
  }
  
  ADREPORT(logit_pred_pi_prop);
  
  if (abundance_component) {
    // Abundance predictions 
    vector<Type> log_pred_abund(X1_pred_ij.rows());
    log_pred_abund = X1_pred_ij * b1_j;
    matrix<Type> pred_abund = exp(log_pred_abund.array());
    
    // Calculate predicted abundance based on higher level groupings
    int n_preds = pred_factor2k_h.size();
    vector<Type> reg_pred_abund(n_pred_levels);
    vector<Type> log_reg_pred_abund(n_pred_levels);
    
    for (int h = 0; h < n_preds; h++) {
      for (int m = 0; m < n_pred_levels; m++) {
        if (pred_factor2k_h(h) == pred_factor2k_levels(m)) {
          reg_pred_abund(m) += pred_abund(h);
          log_reg_pred_abund(m) = log(reg_pred_abund(m));
        }
      }
    }
  
    // Combined predictions
    matrix<Type> pred_abund_mg(n_pred_levels, n_cat);
    
    for (int m = 0; m < n_pred_levels; ++m) {
      for (int k = 0; k < n_cat; ++k) {
        pred_abund_mg(m, k) = reg_pred_abund(m) * pred_pi_prop(m, k);
      }
    }
    matrix<Type> log_pred_abund_mg = log(pred_abund_mg.array());
    
    // Report
    ADREPORT(log_pred_abund);
    ADREPORT(log_reg_pred_abund);
    ADREPORT(log_pred_abund_mg);
  }
  
  return jnll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
