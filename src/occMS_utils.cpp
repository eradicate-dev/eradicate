#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

void getDetVec2(int y, vec& detVec, double mp) {
  if(y == 0) {
    detVec(1) = 1/(1 + exp(mp));
  } else {
    detVec(0) = 0;
    detVec(1) = exp(mp)/(1 + exp(mp));
  }
}


SEXP getSingleDetVec(SEXP y_, SEXP mp_, SEXP K_) {
  int y = as<int>(y_);
  int K = as<int>(K_) + 1;
  vec detVec(K);

  double mp = as<double>(mp_);
  for(int k = 0; k != K; ++k) {
    detVec(k) = 1;
  }
  getDetVec2(y, detVec, mp);
  return wrap(detVec);
}


SEXP getDetVecs(SEXP y_arr, SEXP mp_arr, SEXP J_i, SEXP tin, SEXP K_) {
  vec mp_dims = mp_arr.size();
  int nDMP = mp_dims[0], J = mp_dims[1], nY = mp_dims[2], M = mp_dims[3];
  int K = INTEGER_VALUE(K_) + 1;
  int t = INTEGER_VALUE(tin) - 1;
  void (*getDetVecPtr) (int, double*, double*);
  switch (K) {
  case 2:
    getDetVecPtr = getDetVec2;
    break;
    case 4:
      getDetVecPtr = getDetVec4;
      break;
  }
  SEXP detVec;
  PROTECT(detVec = NEW_NUMERIC(K*M));
  double *mp = NUMERIC_POINTER(mp_arr), *detVecPtr = NUMERIC_POINTER(detVec);
  int *J_ip = INTEGER_POINTER(J_i), *y = INTEGER_POINTER(y_arr);
  int detVec_ind = 0, mp_ind, y_ind;
  for (int i = 0; i != M; ++i) {
    for(int k = 0; k != K; ++k) {
      detVecPtr[detVec_ind + k] = 1;  // initialize to 1.
    }
    for(int j = 0; j != J_ip[i]; ++j) {
      y_ind = i + t*M + j*M*nY;
      mp_ind = j*nDMP + t*nDMP*J + i*nDMP*J*nY;
      if(y[y_ind] != 99) {
        getDetVecPtr(y[y_ind], detVecPtr + detVec_ind, mp + mp_ind);
      }
    }
    detVec_ind += K;
  }
  UNPROTECT(1);
  return detVec;
}
