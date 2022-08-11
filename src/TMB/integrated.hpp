/// @file integrated.hpp

#ifndef integrated_hpp
#define integrated_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// NOTE: identical to LOM_t in negbin.hpp but had to be renamed
// List of matrices 
template <class Type>
struct LOM_tt : vector<matrix<Type> > {
  LOM_tt(SEXP x){  // x = list passed from R
(*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

template<class Type>
Type integrated(objective_function<Type>* obj) {
  
  // DATA ----------------------------------------------------------------------

  // abundance data
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  DATA_IMATRIX(re_index1);
  DATA_IVECTOR(ln_sigma_re_index1);
  DATA_IVECTOR(nobs_re1);
  DATA_STRUCT(Zs, LOM_tt); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix
  DATA_VECTOR(offset_i); // optional offset
  DATA_IVECTOR(rw_index1);   // vector flagging first level of RI for RW
  // for penalized regression splines
  DATA_INTEGER(has_smooths);  // whether or not smooths are included
  DATA_IVECTOR(b_smooth_start);

  // composition data
  DATA_MATRIX(Y2_ik);    // response matrix (a n-by-k matrix)
  DATA_MATRIX(X2_ij);    // covariate matrix (a n-by-j matrix)
  DATA_IMATRIX(re_index2);   // matrix of random intercept levels (a n-by-h matrix) 
  DATA_IVECTOR(ln_sigma_re_index2);  // vector of random intercept deviance estimates (h)
  DATA_IVECTOR(nobs_re2);   // number of random intercepts (h)
  DATA_IVECTOR(rw_index2);   // vector flagging first level of RI for RW
  
  //abundance predictions
  DATA_MATRIX(pred_X1_ij);      // matrix for FE predictions
  DATA_STRUCT(pred_Zs, LOM_tt);  // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(pred_Xs);         // smoother linear effect matrix 
  //composition predictions
  DATA_MATRIX(pred_X2_ij);      // matrix for composition FE predictions

  // shared 
  DATA_INTEGER(random_walk);    // should RIs be random walk
  DATA_INTEGER(has_preds);      // whether or not predictions included
  DATA_INTEGER(n_predX);        // number of predictions (same for both model components)

  // PARAMETERS ----------------------------------------------------------------

  // abundance parameters
  PARAMETER_VECTOR(b1_j); // fixed effects parameters
  PARAMETER(ln_phi);     // variance
  PARAMETER_VECTOR(bs);   // smoother linear effects
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER_VECTOR(re1);
  PARAMETER_VECTOR(ln_sigma_re1);
  
  // composition parameters
  PARAMETER_MATRIX(B2_jk); // parameter matrix
  PARAMETER_VECTOR(re2);  // vector of random intercepts
  PARAMETER_VECTOR(ln_sigma_re2); // among random intercept SD


  // DERIVED QUANTITIES --------------------------------------------------------

  int n1 = y1_i.size();               // number of catch obs
  int n_re1 = re_index1.cols();       // number of abundance RIs
  int n2 = Y2_ik.rows();              // number of composition obs
  int n_re2 = re_index2.cols();       // number of composition RIs
  int n_cat = Y2_ik.cols();           // number of composition categories

  // Matrix for intermediate objects
  matrix<Type> Mu2_ik(n2, n_cat); // matrix of combined fixed/random eff

  Type jnll = 0; // initialize joint log-likelihood


  // LINEAER PREDICTOR ---------------------------------------------------------
  
  // abundance predictor
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


  // Composition linear predictor
  matrix<Type> Mu2_fx_ik = X2_ij * B2_jk;

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
  
  // abundance likelihood
  Type s1, s2;
  for (int i = 0; i < n1; i++) {
    s1 = log(mu_i(i)); // log(mu_i)
    s2 = 2. * s1 - ln_phi; // log(var - mu)
    jnll -= dnbinom_robust(y1_i(i), s1, s2, true);
  }

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

  // composition likelihood
  Type jll = 0; // initialize joint log-likelihood

  for (int i = 0; i <= (n2 - 1); i++) {
    jll = jll + lgamma((n_plus(i) + 1));
    jll = jll + lgamma(Gamma_plus(i));
    jll = jll - lgamma((n_plus(i) + Gamma_plus(i)));
    for(int k = 0; k <= (n_cat - 1); k++){
      jll += lgamma((Y2_ik(i, k) + Gamma(i, k)));
      jll -= lgamma(Gamma(i, k));
      jll -= lgamma((Y2_ik(i, k) + 1));
    }
  }

  jnll -= jll; 
  
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
    
    // abundance predictions
    vector<Type> pred_mu1 = pred_X1_ij * b1_j;
  
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
  
    // combine fixed and smoothed predictions in link space
    for (int m = 0; m < n_predX; m++) {
      pred_mu1(m) += pred_smooth_i(m);
    }

    // composition predictions
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
  
    // combined predictions
    vector<Type> real_pred_mu1 = exp(pred_mu1); // calculate real values for summing
    matrix<Type> pred_mu1_Pi(n_predX, n_cat);
    for (int m = 0; m < n_predX; m++) {
      for (int k = 0; k < n_cat; k++) {
        pred_mu1_Pi(m, k) = real_pred_mu1(m) * pred_Pi_prop(m, k);
      }
    }
    matrix<Type> log_pred_mu1_Pi = log(pred_mu1_Pi.array());
  

    ADREPORT(pred_mu1);
    ADREPORT(pred_Mu2);
    ADREPORT(logit_pred_Pi_prop);
    ADREPORT(log_pred_mu1_Pi);

  }
  
  return jnll;
  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
