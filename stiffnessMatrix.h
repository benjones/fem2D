#pragma once

//sort of a misnomer, it's the stiffnes matrix + mass and damping matrices
#include <vector>
#include "eigenTypedefs.h"

class FemMesh;

class StiffnessMatrix{

public:
  
  StiffnessMatrix() {}
  

  void computeAdjacency(const FemMesh& mesh);

  void assemble(const FemMesh& mesh);

  void checkConsistency();

  size_t findOffset(size_t row, size_t column);

  std::vector<mat2> data;
  std::vector<size_t> columns, offsets;

  

};

