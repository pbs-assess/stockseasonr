#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () 
{
  // load namespace with multivariate distributions
  using namespace density;

  // DATA ----------------------------------------------------------------------

  DATA_MATRIX(Y2_ik);    // response matrix (a n-by-k matrix)
  DATA_MATRIX(X2_ij);    // covariate matrix (a n-by-j matrix)
  DATA_IVECTOR(rfac2);    // vector of random factor levels
  DATA_INTEGER(n_rfac2);  // number of random factor levels
  
  //predictions
  DATA_MATRIX(pred_X2_ij);    // model matrix for predictions
  // DATA_IVECTOR(pred_rfac2); // vector of predicted random intercepts


  // PARAMETERS ----------------------------------------------------------------

  PARAMETER_MATRIX(B2_jk); // parameter matrix
  // PARAMETER_MATRIX(A2_hk);  // matrix of random intercepts (n_rfac2 x n_cat)
  PARAMETER(ln_sigma_A2); // among random intercept SD
  PARAMETER_VECTOR(A2_h);  // vector of random intercepts


  // DERIVED QUANTITIES --------------------------------------------------------

  int n2 = Y2_ik.rows();         // number of observations
  int n_cat = Y2_ik.cols();         // number of categories
  int n_predX2 = pred_X2_ij.rows();   // number of covariates to make predictions on

  // Matrix for intermediate objects
  matrix<Type> Mu2_ik(n2, n_cat); // matrix of combined fixed/random eff

  // Covariance matrix for MVN random intercepts
  matrix<Type> cov_mat(n_cat, n_cat);
  for (int j = 0; j < n_cat; j++) {
    for (int jj = 0; jj < n_cat; jj++) {
      if (j == jj) {
        cov_mat(j, jj) = exp(ln_sigma_A2) * exp(ln_sigma_A2);
      } else {
        cov_mat(j, jj) = 0;
      }
    }
  }

  Type jll = 0; // initialize joint log-likelihood


  // LINEAER PREDICTOR ---------------------------------------------------------
  
  matrix<Type> Mu2_fx_ik = X2_ij * B2_jk; // fixed effects

  for (int i = 0; i < n2; ++i) {
    for(int k = 0; k < n_cat; k++) {
      Mu2_ik(i, k) = Mu2_fx_ik(i, k) + A2_h(rfac2(i));
    }
  }

  matrix<Type> Gamma = exp(Mu2_ik.array()); // add random effect
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

  
  // Probability of multivariate random intercepts
  // for (int h = 0; h < n_rfac2; h++) {
  //   vector<Type> A2_hk_vec = A2_hk.row(h);
  //   MVNORM_t<Type> neg_log_dmvnorm(cov_mat);
  //   jnll += neg_log_dmvnorm(A2_hk_vec);
  // }

  // Probability of random intercepts
  for (int h = 0; h < n_rfac2; h++) {
    jnll -= dnorm(A2_h(h), Type(0.0), exp(ln_sigma_A2), true);
  }

  Type sigma_rfac2 = exp(ln_sigma_A2);
  ADREPORT(sigma_rfac2);
  

  // PREDICTIONS ---------------------------------------------------------------
  
  // matrix<Type> pred_Mu2_fx(n_predX2, n_cat);    //pred fixed effects on log scale
  matrix<Type> pred_Mu2(n_predX2, n_cat);    //pred FE on log scale
  matrix<Type> pred_Gamma(n_predX2, n_cat);  //transformed pred effects 
  vector<Type> pred_Gamma_plus(n_predX2);        
  vector<Type> pred_theta(n_predX2); 
  matrix<Type> pred_Pi(n_predX2, n_cat);      // predicted counts in real 
  vector<Type> pred_n_plus(n_predX2); 
  matrix<Type> pred_Pi_prop(n_predX2, n_cat); // predicted counts as ppn.
  matrix<Type> logit_pred_Pi_prop(n_predX2, n_cat); 

  pred_Mu2 = pred_X2_ij * B2_jk; 

  // for (int m = 0; m < n_predX2; m++) {
  //   for(int k = 0; k < n_cat; k++) {
  //     pred_Mu2(m, k) = pred_Mu2_fx(m, k) + A2_hk(pred_rfac2(m), k);
  //   }
  // }

  pred_Gamma = exp(pred_Mu2.array());
  pred_Gamma_plus = pred_Gamma.rowwise().sum();
  pred_theta = 1 / (pred_Gamma_plus + 1);
  for(int m = 0; m < n_predX2; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_Pi(m, k) = pred_Gamma(m, k) / pred_theta(m);
    }
  }
  pred_n_plus = pred_Pi.rowwise().sum();
  for(int m = 0; m < n_predX2; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_Pi_prop(m, k) = pred_Pi(m, k) / pred_n_plus(m);
      logit_pred_Pi_prop(m, k) = logit(pred_Pi_prop(m, k));
    }
  }

  ADREPORT(pred_Mu2);
  ADREPORT(logit_pred_Pi_prop);
  
  // Return negative loglikelihood
  return jnll;
  
}
