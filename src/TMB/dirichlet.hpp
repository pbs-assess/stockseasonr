/// @file dirichlet.hpp

#ifndef dirichlet_hpp
#define dirichlet_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type dirichlet(objective_function<Type>* obj) {

  // load namespace with multivariate distributions
  // using namespace density;

  // DATA ----------------------------------------------------------------------

  DATA_MATRIX(Y2_ik);    // response matrix (a n-by-k matrix)
  DATA_MATRIX(X2_ij);    // covariate matrix (a n-by-j matrix)
  // DATA_IVECTOR(rfac2);    // vector of random factor levels
  // DATA_INTEGER(n_rfac2);  // number of random factor levels
  DATA_IMATRIX(re_index2);   // matrix of random intercept levels (a n-by-h matrix) 
  DATA_IVECTOR(ln_sigma_re_index2);  // vector of random intercept deviance estimates (h)
  DATA_IVECTOR(nobs_re2);   // number of random intercepts (h)
  DATA_IVECTOR(rw_index2);   // vector flagging first level of RI for RW

  //predictions
  DATA_MATRIX(pred_X2_ij);    // model matrix for predictions

  DATA_INTEGER(random_walk); // should RIs be random walk
  DATA_INTEGER(has_preds);  // whether or not predictions included
  DATA_INTEGER(n_predX);        // number of predictions (same for both model components)


  // PARAMETERS ----------------------------------------------------------------

  PARAMETER_MATRIX(B2_jk); // parameter matrix
  PARAMETER_VECTOR(re2);  // vector of random intercepts
  PARAMETER_VECTOR(ln_sigma_re2); // among random intercept SD


  // DERIVED QUANTITIES --------------------------------------------------------

  int n2 = Y2_ik.rows();         // number of observations
  int n_re2 = re_index2.cols();      // number of random intercepts
  int n_cat = Y2_ik.cols();         // number of categories
  // int n_predX2 = pred_X2_ij.rows();   // number of covariates to make predictions on

  // Matrix for intermediate objects
  matrix<Type> Mu2_ik(n2, n_cat); // matrix of combined fixed/random eff

  Type jll = 0; // initialize joint log-likelihood


  // LINEAER PREDICTOR ---------------------------------------------------------
  
  matrix<Type> Mu2_fx_ik = X2_ij * B2_jk; // fixed effects

  // Add random intercepts
  vector<Type> eta_re2_i(n2);
  eta_re2_i.setZero();
  vector<Type> mu2_i(n2);
  mu2_i.setZero();
  for (int i = 0; i < n2; ++i) {
    int temp = 0;
    for (int g = 0; g < n_re2; g++) {
      if (g == 0) eta_re2_i(i) += re2(re_index2(i, g));
      if (g > 0) {
        temp += nobs_re2(g - 1);
        eta_re2_i(i) += re2(re_index2(i, g) + temp);
      } 
    }
    for (int k = 0; k < n_cat; k++) {
      Mu2_ik(i, k) = Mu2_fx_ik(i, k) + eta_re2_i(i);
    }
   }

  matrix<Type> Gamma = exp(Mu2_ik.array()); 
  vector<Type> n_plus = Y2_ik.rowwise().sum(); // row sum of response
  vector<Type> Gamma_plus = Gamma.rowwise().sum(); // row sum of gamma
  
  
  // LIKELIHOOD ----------------------------------------------------------------

  for(int i = 0; i <= (n2 - 1); i++){
    jll = jll + lgamma((n_plus(i) + 1));
    jll = jll + lgamma(Gamma_plus(i));
    jll = jll - lgamma((n_plus(i) + Gamma_plus(i)));
    for(int k = 0; k <= (n_cat - 1); k++){
      jll += lgamma((Y2_ik(i, k) + Gamma(i, k)));
      jll -= lgamma(Gamma(i, k));
      jll -= lgamma((Y2_ik(i, k) + 1));
    }
  }

  Type jnll;
  jnll = -jll;

  
  // Probability of random intercepts
  if (random_walk) {
    for (int h = 0; h < re2.size(); h++) {
      if (rw_index2(h) == 0) {
        jnll -= dnorm(re2(h), Type(0.0), exp(ln_sigma_re2(ln_sigma_re_index2(h))), true);  
      }
      if (rw_index2(h) > 0) {
        jnll -= dnorm(re2(h), re2(h - 1), exp(ln_sigma_re2(ln_sigma_re_index2(h))), true);
      }
    } 
  } else {
    for (int h = 0; h < re2.size(); h++) {
      jnll -= dnorm(re2(h), Type(0.0), exp(ln_sigma_re2(ln_sigma_re_index2(h))), true);
    }
  }
  
  
  // PREDICTIONS ---------------------------------------------------------------
  
  if (has_preds) {
    matrix<Type> pred_Mu2(n_predX, n_cat);    //pred FE on log scale
    matrix<Type> pred_Gamma(n_predX, n_cat);  //transformed pred effects 
    vector<Type> pred_Gamma_plus(n_predX);        
    vector<Type> pred_theta(n_predX); 
    matrix<Type> pred_Pi(n_predX, n_cat);      // predicted counts in real 
    vector<Type> pred_n_plus(n_predX); 
    matrix<Type> pred_Pi_prop(n_predX, n_cat); // predicted counts as ppn.
    matrix<Type> logit_pred_Pi_prop(n_predX, n_cat); 
  
    pred_Mu2 = pred_X2_ij * B2_jk; 
  
    pred_Gamma = exp(pred_Mu2.array());
    pred_Gamma_plus = pred_Gamma.rowwise().sum();
    pred_theta = 1 / (pred_Gamma_plus + 1);
    for(int m = 0; m < n_predX; m++) {
      for(int k = 0; k < n_cat; k++) {
        pred_Pi(m, k) = pred_Gamma(m, k) / pred_theta(m);
      }
    }
    pred_n_plus = pred_Pi.rowwise().sum();
    for(int m = 0; m < n_predX; m++) {
      for(int k = 0; k < n_cat; k++) {
        pred_Pi_prop(m, k) = pred_Pi(m, k) / pred_n_plus(m);
        logit_pred_Pi_prop(m, k) = logit(pred_Pi_prop(m, k));
      }
    }
  
    ADREPORT(pred_Mu2);
    ADREPORT(logit_pred_Pi_prop);
  }
  
  // Return negative loglikelihood
  return jnll;
  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
