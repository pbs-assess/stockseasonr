// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_stockseasonr_TMBExports
#include <TMB.hpp>
#include "dirichlet.hpp"
#include "integrated.hpp"
#include "nb_dirchlet_1re.hpp"
#include "negbin.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "dirichlet") {
    return dirichlet(this);
  } else if(model == "integrated") {
    return integrated(this);
  } else if(model == "nb_dirchlet_1re") {
    return nb_dirchlet_1re(this);
  } else if(model == "negbin") {
    return negbin(this);
  } else {
    Rf_error("Unknown model.");
  }
  return 0;
}
