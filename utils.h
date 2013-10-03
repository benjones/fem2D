#include "eigenTypedefs.h"

namespace Utils{

  inline void svd2(const mat2 A, mat2& U, mat2& V, 
		   vec2& sigma){
    
    //check if orthogonal:
    float norm1, norm2;
    norm1 = A.col(0).norm();
    norm2 = A.col(1).norm();
    //std::cout << "A: " << A << std::endl;
    //std::cout << "n1: " << norm1 << " norm2: " << norm2 << std::endl;
    //std::cout << "dotproduct: " << A.col(0).transpose()*A.col(1) << std::endl;
    if( (fabs(norm1 - norm2) < .0001) &&
	(fabs(A.col(0).transpose()*A.col(1)) < .0001)){
      //A is orthogonal
      
      U = A;
      U.col(0) *= (1.0f/norm1);
      U.col(1) *= (1.0f/norm2);
      sigma = vec2(norm1, norm2); //n1 == n2
      V = mat2::Identity();
      if(A.determinant() < 0){
	//std::cout << "flipping sign" << std::endl;
	sigma(1) *= -1;
	U.col(1) *= -1;
      } 
      
      /*std::cout << "orthogonal matrix: " << A << std::endl
		<< "U : " << U << std::endl
		<< "V: " << V << std::endl
		<< "Sigma: " << sigma << std::endl
		<< "U UT " << U*U.transpose() << std::endl
		<< "V VT " << V*V.transpose() << std::endl
		<< "U*sigma*V: " << U*sigma.asDiagonal()*V << std::endl;*/
      return;
    }


    mat2 uProd = A*A.transpose();
    mat2 vProd = A.transpose()*A;
    //uRotation
    float theta = 0.5f*atan2(2*uProd(1,0), uProd(0,0) - uProd(1,1));
    // v rotation
    float phi = 0.5f*atan2(2*vProd(1,0), vProd(0,0) - vProd(1,1));
    //cossinf optimization would be nice
    float ct = cosf(theta);
    float st = sinf(theta);
    float cp = cosf(phi);
    float sp = sinf(phi);
    
    U << ct , -st , st , ct;
    V << cp , -sp , sp , cp;
    
    //compute singular values:
    float sqrtTerm = sqrtf(pow(uProd(0,0) - uProd(1,1),2) + 4*uProd(0,1)*uProd(0,1));
    sigma(0) = sqrtf(0.5f*(uProd(0,0) + uProd(1,1) + sqrtTerm));
    sigma(1) = sqrtf(0.5f*(uProd(0,0) + uProd(1,1) - sqrtTerm));

    //compute signs of the SVs
    mat2 signMatrix = U.transpose()*A*V;
    bool n0, n1; //negative sv 0, negative sv 1
    n0 = (signMatrix(0,0) < 0 ) ? true : false;
    n1 = (signMatrix(1,1) < 0 ) ? true : false;

    if(n0){
      if(n1){
	//just negate the matrix (lets say V)
	V *= -1;
      } else {
	//negate V, and sigma2
	V *= -1;
	sigma(1) *= -1;	  
      }
    } else if (n1){
      sigma(1) *= -1; //smaller singular value needs to be negative
    } 
    
    //looks like everythign is working correctly
    /*if( ! (U*sigma.asDiagonal()*V.transpose()).isApprox(A) ){
      std::cout << "SVD not close: " << U*sigma.asDiagonal()*V.transpose() - A << std::endl;
      }*/

    //assert((U*sigma.asDiagonal()*V.transpose() - A).norm() < .001);
    //std::cout << "U*UT" << U*U.transpose() << "V*VT " << V*V.transpose() << std::endl;
  }



  inline void RSDecomp(const mat2& A, mat2& R, mat2& S){

    if(A.norm() < 1e-12){
      R = mat2::Identity();
      S = mat2::Zero();
      return;
    }
    
    auto detA = A.determinant();
    
    if(true ||detA <= 0.001){
      //use the SVD
      mat2 U, V;
      vec2 sigma;
      svd2(A, U, V, sigma);
      vec2 c(1, (U*V.transpose()).determinant());
      R = U*c.asDiagonal()*V.transpose();
      
    } else {
      
      R = A + 
	(mat2() << A(1,1), -A(1,0), -A(0,1), A(0,0)).finished();
      R.col(0) *= 1.0/R.col(0).norm();
      R.col(1) *= 1.0/R.col(1).norm();
      
      //assert(fabs(R.determinant() - 1) < 0.0001);
      if(R.determinant() < 0.9 ||  !(R*R.transpose()).isApprox(mat2::Identity())){
	std::cout << "bad polor of " << A << std::endl
		  << "R: " << R << std::endl
		  << "Rdet: " << R.determinant() << std::endl;
	exit(1);
      }
      
    }
    
    S = R.transpose()*A;
  }



    

}
