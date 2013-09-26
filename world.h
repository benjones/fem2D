#pragma once
#include <vector>
#include "femMesh.h"
#include "planeObstacle.h"

class World{
public:
  std::vector<FemMesh> femObjects;
  std::vector<PlaneObstacle> planeObstacles;

  double dt;


  World(std::string filename);
  
  void integrate();
  void renderOpenGL();

};
