#include "eigenTypedefs.h"

namespace Utils{

  inline void RSDecomp(const mat2& A, mat2& R, mat2& S){
    
    if(A.norm() < 1e-12){
      R = mat2::Identity();
      S = mat2::Zero();
      return;
    }

    auto detA = A.determinant();
    R = A + (detA < 0 ? -1.0f : 1.0f)*
      (mat2() << A(1,1), -A(1,0), -A(0,1), A(0,0)).finished();
    R.col(0) *= 1.0/R.col(0).norm();
    R.col(1) *= 1.0/R.col(1).norm();

    S = R.transpose()*A;
  }


}
