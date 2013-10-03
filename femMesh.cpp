#include "femMesh.h"
#include "range.cpp"
#include "utils.h"
#include "planeObstacle.h"
#include "cghs.h"
#include <SDL_opengl.h>


void FemMesh::initialize(const std::vector<double>& vertices,
			 const std::vector<unsigned>& triangles,
			 MaterialParameters _materialParameters){
  materialParameters = _materialParameters;

  nodes.resize(vertices.size()/2);
  for(auto i : range(nodes.size())){
    nodes[i].worldPos(0) = vertices[2*i];
    nodes[i].worldPos(1) = vertices[2*i+1];

    nodes[i].referencePos = nodes[i].worldPos;
    nodes[i].velocity.setZero();
    nodes[i].mass = 0;
  }

  elements.resize(triangles.size()/3);
  for(auto i : range(elements.size())){
    elements[i].nodes[0] = triangles[3*i];
    elements[i].nodes[1] = triangles[3*i + 1];
    elements[i].nodes[2] = triangles[3*i + 2];
  }

  //compute mass and beta matrices for each element

  for(auto &e : elements){

    auto& u0 = nodes[e.nodes[0]].referencePos;
    auto& u1 = nodes[e.nodes[1]].referencePos;
    auto& u2 = nodes[e.nodes[2]].referencePos;
    
    e.beta.col(0) = u1 - u0;
    e.beta.col(1) = u2 - u0;

    //area (length) weighted face normals
    //since I'm not normalizing, the vectors are already "area weighted"
    e.faceNormals[0] = /*(u1 - u0).norm()* */vec2(u1.y() - u0.y(), u0.x() - u1.x());
    e.faceNormals[1] = /*(u2 - u1).norm()* */vec2(u2.y() - u1.y(), u1.x() - u2.x());
    e.faceNormals[2] = /*(u0 - u2).norm()* */vec2(u0.y() - u2.y(), u2.x() - u0.x());

    
    double area = 0.5*(e.beta.col(0)(0)*e.beta.col(1)(1) - 
		       e.beta.col(0)(1)*e.beta.col(1)(0));
    
    for(auto i : range(3)){
      nodes[e.nodes[i]].mass += area*materialParameters.density;
    }
    
    e.beta = e.beta.inverse().eval();
    
    //compute the un-rotated stiffness matrix stuff
    //outer product matrices we'll reuse below
    mat2 outers[6] = { e.faceNormals[1]*e.faceNormals[1].transpose(), //0,0
		       e.faceNormals[1]*e.faceNormals[2].transpose(), //0, 1
		       e.faceNormals[1]*e.faceNormals[0].transpose(), //0, 2
		       e.faceNormals[2]*e.faceNormals[2].transpose(), //1, 1
		       e.faceNormals[2]*e.faceNormals[0].transpose(), //1, 2
		       e.faceNormals[0]*e.faceNormals[0].transpose()}; //2, 2
		       
    double dots[6] = {e.faceNormals[1].dot(e.faceNormals[1]), //0,0
		      e.faceNormals[1].dot(e.faceNormals[2]), //0, 1
		      e.faceNormals[1].dot(e.faceNormals[0]), //0, 2
		      e.faceNormals[2].dot(e.faceNormals[2]), //1, 1
		      e.faceNormals[2].dot(e.faceNormals[0]), //1, 2
		      e.faceNormals[0].dot(e.faceNormals[0])}; //2, 2
		      
    for(auto i : range(6)){
      e.stiffnessMatrices[i] = (1/area)*
	materialParameters.lambda*outers[i] +
	materialParameters.mu*(dots[i]*mat2::Identity() + 
			       outers[i].transpose());
      
    }
    //add the stiffnessmatrices*rest_pos to each node's force offsets
    //force offsets
    e.forceOffset[0] = e.stiffnessMatrices[0]*u0 +
      e.stiffnessMatrices[1]*u1 + e.stiffnessMatrices[2]*u2;
    e.forceOffset[1] = e.stiffnessMatrices[1].transpose()*u0 +
      e.stiffnessMatrices[3]*u1 + e.stiffnessMatrices[4]*u2;
    e.forceOffset[2] = e.stiffnessMatrices[2].transpose()*u0 +
      e.stiffnessMatrices[4].transpose()*u1 + e.stiffnessMatrices[5]*u2;
  }
  
  for(auto & n : nodes){
    n.mass /= 3;
  }

  stiffnessMatrix.computeAdjacency(*this);


  
}



void FemMesh::renderOpenGL(){

  
  glColor3f(0,1,0);
  glBegin(GL_LINES);
  for(auto& e : elements){
    for(auto i : range(3)){
      
      glVertex2d(nodes[e.nodes[i]].worldPos(0),
		 nodes[e.nodes[i]].worldPos(1));
      glVertex2d(nodes[e.nodes[(i+1)%3]].worldPos(0),
		 nodes[e.nodes[(i+1)%3]].worldPos(1));
    }
  }
  glEnd();

}

void FemMesh::setForcesToGravity(){
  //clear forces:
  for(auto& n : nodes){
    if(true){ //gravity
      n.force = vec2(0, -9.81)* n.mass;
    } else {
      n.force.setZero();
    }
  }
}

void FemMesh::computeDeformationGradient(bool isExplicit){

  //compute elastic/damping forces
  for(auto & e : elements){

    auto& x0 = nodes[e.nodes[0]].worldPos;
    auto& x1 = nodes[e.nodes[1]].worldPos;
    auto& x2 = nodes[e.nodes[2]].worldPos;
    
    

    mat2 xMatrix;
    xMatrix.col(0) = x1 - x0;
    xMatrix.col(1) = x2 - x0;

    mat2 vMatrix;

    auto deformationGradient = xMatrix*e.beta;
    
    mat2 S;
    Utils::RSDecomp(deformationGradient, e.elementRotation, S);
    
    if(isExplicit){
    
      auto& v0 = nodes[e.nodes[0]].velocity;
      auto& v1 = nodes[e.nodes[1]].velocity;
      auto& v2 = nodes[e.nodes[2]].velocity;
      
      vMatrix.col(0) = v1 - v0;
      vMatrix.col(1) = v2 - v0;
      
      auto velocityGradient = vMatrix*e.beta;
      
      
      
      auto fTwiddle = e.elementRotation.transpose()*deformationGradient;
      auto vTwiddle = e.elementRotation.transpose()*velocityGradient;
      
      auto strain = 0.5*(fTwiddle + fTwiddle.transpose()) - mat2::Identity();
      auto strainRate = 0.5*(vTwiddle +
			     vTwiddle.transpose());
      
      auto stress = materialParameters.lambda*strain.trace()*mat2::Identity() +
	2*materialParameters.mu*strain;
      
      auto viscousStress = 
	materialParameters.dampLambda*strainRate.trace()*mat2::Identity() +
	2*materialParameters.dampMu*strainRate;
      
      auto forceProduct = e.elementRotation*(stress + viscousStress);
      for(auto i : range(3)){
	nodes[e.nodes[i]].force += forceProduct*e.faceNormals[(i+1)%3];
      }
    }
  }
  
}


void FemMesh::integrate(double dt){
  this->dt = dt;
   
  setForcesToGravity();
  
  computeDeformationGradient(true);
  
  //compute deformation gradient

  //forward symplectic euler
  //update velocities
  for(auto& n : nodes){
    n.velocity += n.force*(dt/n.mass);
    n.worldPos += n.velocity*dt;
  }
  
  constrainNodes();
}

void FemMesh::constrainNodes(){
  //fix 2 vertices of the triangle for testing
  for(auto i : {2}){
    nodes[i].worldPos = nodes[i].referencePos;
    nodes[i].velocity.setZero();
  }
}

void FemMesh::collidePlanes(const std::vector<PlaneObstacle>& planes){

  for(const auto& p : planes){
    for(auto& n : nodes){
      
      if(n.worldPos.dot(p.normal) < p.offset){

	//BAD BAD BAD
	n.worldPos.y() =0;//*= -1;
	n.velocity.y() *= -1;
      }
    }
  }

}


void mult(const FemMesh& mesh, double* x, double* y){
  
  auto mFactor = (1 + mesh.dt*
		  mesh.materialParameters.raleighAlpha);
  auto kFactor = mesh.dt*(mesh.dt + 
			  mesh.materialParameters.raleighBeta);
  
  for(auto i : range(mesh.nodes.size())){
    auto vi = vec2(x[2*i], x[2*i +1]);
    vec2 entry = mesh.nodes[i].mass*mFactor*vi;
    
    for(auto j : range(mesh.stiffnessMatrix.offsets[i],
		       mesh.stiffnessMatrix.offsets[i+1])){
      auto vj = vec2(x[2*mesh.stiffnessMatrix.columns[j]], 
		     x[2*mesh.stiffnessMatrix.columns[j]+1]);
      entry += mesh.stiffnessMatrix.data[j]*kFactor*vj;
    }
    y[2*i]   = entry(0);
    y[2*i+1] = entry(1);
  }


}


void FemMesh::integrateBackwardsEuler(double dt){
  this->dt = dt;

  setForcesToGravity();
  computeDeformationGradient(false);

  stiffnessMatrix.assemble(*this);
  
  //fill in rhs = dt*fExt + M*v - dt K (x - u)
  rhsSolve.resize(2*nodes.size());
  velocityGuess.resize(2*nodes.size());
  for(auto i : range(nodes.size())){

    rhsSolve[2*i] = dt*nodes[i].force(0) +
      nodes[i].mass*nodes[i].velocity(0);
    
    rhsSolve[2*i+1] = dt*nodes[i].force(1) +
      nodes[i].mass*nodes[i].velocity(1);

    velocityGuess[2*i    ] = nodes[i].velocity(0);
    velocityGuess[2*i + 1] = nodes[i].velocity(1);
  }

  //probably slow, and should get mixed in with some other step
  for(auto& e : elements){
    for(auto i : range(3)){
      vec2 rotatedOffset = e.elementRotation*e.forceOffset[i];
      rhsSolve[2*e.nodes[i]   ] += dt*rotatedOffset(0);
      rhsSolve[2*e.nodes[i] +1] += dt*rotatedOffset(1);
    }
  }

  //loop over the stiffness matrx for the last bit
  vec2 entry;

  for(auto row : range(nodes.size())){
    for(auto j : range(stiffnessMatrix.offsets[row],
		       stiffnessMatrix.offsets[row+1])){
      entry = -dt* ( stiffnessMatrix.data[j]*
		    (nodes[stiffnessMatrix.columns[j]].worldPos));// -
		    //		     nodes[stiffnessMatrix.columns[j]].referencePos));
      rhsSolve[2*row   ] += entry(0);
      rhsSolve[2*row +1] += entry(1);
    }
  }


  
  //do the solve
  cghsIterationLimit(nodes.size()*2, *this, &(rhsSolve[0]), &(velocityGuess[0]), 1e-10, 30, false);
  
  //do the update
  for(auto i : range(nodes.size())){
    nodes[i].velocity(0) = velocityGuess[2*i];
    nodes[i].velocity(1) = velocityGuess[2*i +1];
    nodes[i].worldPos += dt*nodes[i].velocity;
  }


  constrainNodes();

}



