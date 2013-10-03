
#include "stiffnessMatrix.h"
#include "femMesh.h"
#include "range.h"

void StiffnessMatrix::computeAdjacency(const FemMesh& mesh){

  columns.clear();
  offsets.resize(mesh.nodes.size()+1);

  std::vector<std::vector<size_t>> neighbors(mesh.nodes.size());
  
  for(auto& e : mesh.elements){
    for(auto i : range(3)){
      neighbors[e.nodes[i]].push_back(e.nodes[i]);
      for(auto j : range(3)){
	if(e.nodes[i] < e.nodes[j]){
	  neighbors[e.nodes[i]].push_back(e.nodes[j]);
	  neighbors[e.nodes[j]].push_back(e.nodes[i]);
	}
      }
    }
  }

  auto columnsSize = 0;
  for(auto i : range(mesh.nodes.size())){
    columnsSize += neighbors[i].size();
  }
  
  columns.resize(columnsSize);

  offsets[0] = 0;
  for(auto i : range(mesh.nodes.size())){
    offsets[i+1] = offsets[i] + neighbors[i].size();
    std::partial_sort_copy(neighbors[i].begin(), neighbors[i].end(),
			   columns.begin() + offsets[i],
			   columns.begin() + offsets[i+1]);
  }  
}


size_t StiffnessMatrix::findOffset(size_t row, size_t column){
  
  auto it = std::lower_bound(columns.begin() + offsets[row],
			     columns.begin() + offsets[row+1],
			     column);

  assert(*it == column);
  return std::distance(columns.begin(), it);
  

}


void StiffnessMatrix::assemble(const FemMesh& mesh){

  data.assign(columns.size(), mat2::Zero());
  for(auto& e : mesh.elements){
    auto& R = e.elementRotation;
    auto& RT = R.transpose();

    //0, 0
    mat2 entry = R*e.stiffnessMatrices[0]*RT;
    data[findOffset(e.nodes[0], e.nodes[0])] += entry;
    //0, 1 and 1, 0
    entry = R*e.stiffnessMatrices[1]*RT;
    data[findOffset(e.nodes[0], e.nodes[1])] += entry;
    data[findOffset(e.nodes[1], e.nodes[0])] += entry.transpose();
    //0, 2 and 2, 0
    entry = R*e.stiffnessMatrices[2]*RT;
    data[findOffset(e.nodes[0], e.nodes[2])] += entry;
    data[findOffset(e.nodes[2], e.nodes[0])] += entry.transpose();
    //1, 1
    entry = R*e.stiffnessMatrices[3]*RT;
    data[findOffset(e.nodes[1], e.nodes[1])] += entry;
    //1, 2 and 2,1
    entry = R*e.stiffnessMatrices[4]*RT;
    data[findOffset(e.nodes[1], e.nodes[2])] += entry;
    data[findOffset(e.nodes[2], e.nodes[1])] += entry.transpose();
    //2,2
    entry = R*e.stiffnessMatrices[5]*RT;
    data[findOffset(e.nodes[2], e.nodes[2])] += entry;


  }

  //checkConsistency();
}


void StiffnessMatrix::checkConsistency(){
  
  for(auto i : range(offsets.size() -1)){
    mat2 sum = mat2::Zero();
    for(auto j : range(offsets[i], offsets[i+1])){
      sum += data[j];
    }
    if(!sum.isZero(0.0001)){
      std::cout << "i : " << i << std::endl;
      std::cout << "sum: " << sum << std::endl;
      exit(1);
    }
  }  
  std::cout << "consistent" << std::endl;
}



