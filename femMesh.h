#pragma once

#include <vector>
#include "eigenTypedefs.h"

class PlaneObstacle;

class FemNode {
public:
  vec2 referencePos, worldPos, velocity;
  double mass;
  vec2 force;

};

class FemElement {
public:
  unsigned nodes[3];
  
  vec2 faceNormals[3];

  mat2 beta;

};


struct MaterialParameters{
  double density, lambda, mu, dampLambda, dampMu;
};

class FemMesh {

public:
  FemMesh() {}

  void initialize (const std::vector<double>& vertices,
		   const std::vector<unsigned>& triangles,
		   MaterialParameters _materialParameters);

  void integrate(double dt);
  void collidePlanes(const std::vector<PlaneObstacle>& planes);

  void renderOpenGL();


  std::vector<FemNode> nodes;
  std::vector<FemElement> elements;
  MaterialParameters materialParameters;

};
