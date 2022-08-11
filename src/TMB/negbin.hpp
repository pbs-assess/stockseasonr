/// @file negbin.hpp

#ifndef negbin_hpp
#define negbin_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// List of matrices
template <class Type>
struct LOM_t : vector<matrix<Type> > {
  LOM_t(SEXP x){  // x = list passed from R
(*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

template<class Type>
Type negbin(objective_function<Type>* obj) {

  // DATA ----------------------------------------------------------------------
  
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  DATA_IMATRIX(re_index1);
  DATA_IVECTOR(ln_sigma_re_index1);
  DATA_IVECTOR(nobs_re1);
  DATA_STRUCT(Zs, LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix
  DATA_VECTOR(offset_i); // optional offset
  DATA_IVECTOR(rw_index1);   // vector flagging first level of RI for RW
  // for penalized regression splines
  DATA_INTEGER(has_smooths);  // whether or not smooths are included
  DATA_IVECTOR(b_smooth_start);
  
  //predictions
  DATA_MATRIX(pred_X1_ij); // matrix for FE predictions
  DATA_STRUCT(pred_Zs, LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(pred_Xs); // smoother linear effect matrix 
  // DATA_IVECTOR(pred_re_index1);

  DATA_INTEGER(random_walk); // should RIs be random walk
  DATA_INTEGER(has_preds);  // whether or not predictions included
  DATA_INTEGER(n_predX);        // number of predictions (same for both model components)

  // PARAMETERS ----------------------------------------------------------------
  
  PARAMETER_VECTOR(b1_j); // fixed effects parameters
  PARAMETER(ln_phi);     // variance
  PARAMETER_VECTOR(bs);   // smoother linear effects
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included
  // random effects
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER_VECTOR(re1);
  PARAMETER_VECTOR(ln_sigma_re1);
  

  // DERIVED -------------------------------------------------------------------

  int n1 = y1_i.size();
  int n_re1 = re_index1.cols();      // number of random intercepts
  // int n_predX = pred_X1_ij.rows(); // number of predictions   

  Type jnll = 0.0; // initialize joint negative log likelihood


  // LINEAER PREDICTOR ---------------------------------------------------------

  vector<Type> eta_fx_i = X1_ij * b1_j;

  // Smooths
  vector<Type> eta_smooth_i(X1_ij.rows());
  eta_smooth_i.setZero();

  if (has_smooths) {
    for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
      vector<Type> beta_s(Zs(s).cols());
      beta_s.setZero();
      for (int j = 0; j < beta_s.size(); j++) {
        beta_s(j) = b_smooth(b_smooth_start(s) + j);
        jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
      }
      eta_smooth_i += Zs(s) * beta_s;
    }
    eta_smooth_i += Xs * bs;
  }

  // Combine smooths, linear and offset
  vector<Type> eta_i(n1);
  eta_i.setZero();
  for (int i = 0; i < n1; i++) {
    eta_i(i) = eta_fx_i(i) + eta_smooth_i(i) + offset_i(i); 
  }    
  
  
 // Add random intercepts
 vector<Type> eta_re_i(n1);
 eta_re_i.setZero();
 vector<Type> mu_i(n1);
 mu_i.setZero();
 for (int i = 0; i < n1; i++) {
    int temp = 0;
    for (int g = 0; g < n_re1; g++) {
      if (g == 0) eta_re_i(i) += re1(re_index1(i, g));
      if (g > 0) {
        temp += nobs_re1(g - 1);
        eta_re_i(i) += re1(re_index1(i, g) + temp);
      } 
    }
   mu_i(i) = exp(eta_i(i) + eta_re_i(i));
 }


  // LIKELIHOOD ----------------------------------------------------------------

  Type s1, s2;
  for(int i = 0; i < n1; i++){
    s1 = log(mu_i(i)); // log(mu_i)
    s2 = 2. * s1 - ln_phi; // log(var - mu)
    jnll -= dnbinom_robust(y1_i(i), s1, s2, true);
  }
  // std::cout << "HERE-3" << "\n";

  // Report for residuals
  // ADREPORT(s1);
  // ADREPORT(s2);

  // Probability of random coefficients (add counter for random walk)
  // for(int h = 0; h < re1.size(); h++){
  //   // if (h == 0) {
  //       jnll -= dnorm(re1(h), Type(0.0), exp(ln_sigma_re1(ln_sigma_re_index1(h))), true);  
  // }
  // Probability of random intercepts
  if (random_walk) {
    for (int h = 0; h < re1.size(); h++) {
      if (rw_index1(h) == 0) {
        jnll -= dnorm(re1(h), Type(0.0), exp(ln_sigma_re1(ln_sigma_re_index1(h))), true);  
      }
      if (rw_index1(h) > 0) {
        jnll -= dnorm(re1(h), re1(h - 1), exp(ln_sigma_re1(ln_sigma_re_index1(h))), true);
      }
    } 
  } else {
    for (int h = 0; h < re1.size(); h++) {
      jnll -= dnorm(re1(h), Type(0.0), exp(ln_sigma_re1(ln_sigma_re_index1(h))), true);
    }
  }


  // PREDICTIONS ---------------------------------------------------------------

  if (has_preds) {
    vector<Type> pred_mu1 = pred_X1_ij * b1_j;
  
    // smoothers
    vector<Type> pred_smooth_i(n_predX);
    pred_smooth_i.setZero();
    
    if (has_smooths) {
      for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
        vector<Type> beta_s(pred_Zs(s).cols());
        beta_s.setZero();
        for (int j = 0; j < beta_s.size(); j++) {
          beta_s(j) = b_smooth(b_smooth_start(s) + j);
        }
        pred_smooth_i += pred_Zs(s) * beta_s;
      }
      pred_smooth_i += pred_Xs * bs;
    }
  
    // combine fixed and smoothed predictions
    for (int i = 0; i < n_predX; i++) {
      pred_mu1(i) += pred_smooth_i(i);
    }
  
    // // add random intercepts 
    // for (int i = 0; i < n_predX; i++) {
    //   pred_mu1(i) += re1(pred_re_index1(i));
    // }    
  
    REPORT(pred_mu1);
    ADREPORT(pred_mu1);
  }
  
  return jnll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
