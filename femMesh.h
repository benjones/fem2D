#pragma once

#include <vector>
#include "eigenTypedefs.h"
#include "stiffnessMatrix.h"

class PlaneObstacle;

class FemNode {
public:
  vec2 referencePos, worldPos, velocity;
  double mass;
  vec2 force;

  vec2 tensileForces, compressiveForces;

  mat2 tensileTensor, compressiveTensor;

  std::vector<size_t> triangles;

};

class FemElement {
public:
  unsigned nodes[3];
  
  vec2 faceNormals[3];

  mat2 beta;
  mat2 bMatrix;
  
  mat2 stiffnessMatrices[6];
  //00, 01, 02, 11, 12, 22, symmetric, so only store these

  mat2 stress, elementRotation;

  double area;

  vec2 forceOffset[3];
};


struct MaterialParameters{
  double density, lambda, mu, dampLambda, dampMu;
  double raleighAlpha, raleighBeta;
  double toughness;
  //C = alpha*massMatrix + beta*stiffnessMatrix
};

class FemMesh {

public:
  FemMesh() {}

  void initialize (const std::vector<double>& vertices,
		   const std::vector<unsigned>& triangles,
		   MaterialParameters _materialParameters);

  void integrate(double dt);

  void integrateBackwardsEuler(double dt);

  void collidePlanes(const std::vector<PlaneObstacle>& planes);

  void renderOpenGL();


  void setForcesToGravity();
  void computeDeformationGradient(bool isExplicit);
  void constrainNodes();
  void doFracture();


  std::vector<FemNode> nodes;
  std::vector<FemElement> elements;
  MaterialParameters materialParameters;

  StiffnessMatrix stiffnessMatrix;

  std::vector<double> rhsSolve, velocityGuess;
  
  std::vector<size_t> constrainedNodes;

  double dt;

};

//compute y = Ax, with A = M*(1+dt*alpha) + K(dt^2 + dt*beta)
void mult(const FemMesh& mesh, double* x, double* y);
