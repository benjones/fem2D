#include "readNodeAndEle.h"

#include <fstream>
#include <cassert>

#include "range.h"

bool readNodeFile(std::string filename,
		  std::vector<double>& vertices);

bool readEleFile(std::string filename,
		 std::vector<unsigned>& triangles);

bool readNodeAndEle(std::string filename,
		    std::vector<double>& vertices,
		    std::vector<unsigned>& triangles){
  


  return readNodeFile(filename + ".node", vertices) &&
    readEleFile(filename + ".ele", triangles);

  
  
}


bool readNodeFile(std::string filename, std::vector<double>& vertices){

  std::ifstream ins(filename.c_str());

  unsigned numVertices, dimension, numAttributes, border;

  ins >> numVertices >> dimension >> numAttributes >> border;
  assert(dimension == 2);
  
  vertices.resize(2*numVertices);
  unsigned trash;
  for(auto i : range(numVertices)){
    ins >> trash >> vertices[2*i] >> vertices[2*i +1];
    if(border) ins >> trash; //get border attribute too
  }

  //for(auto& v : vertices){
  //  std::cout << v << " " << std::endl;
  //}
  return true;
}


bool readEleFile(std::string filename, std::vector<unsigned>& triangles){

  std::ifstream ins(filename.c_str());

  unsigned numTriangles, nodesPerTriangle, numAttributes;
  
  ins >> numTriangles >> nodesPerTriangle >> numAttributes;
  assert(nodesPerTriangle == 3);

  triangles.resize(3*numTriangles);
  //assume 1- indexed
  unsigned triNumber;
  for(auto i : range(numTriangles)){
    ins >> triNumber >> triangles[3*i] >> triangles[3*i +1] >> triangles[3*i +2];
    --triangles[3*i];
    --triangles[3*i + 1];
    --triangles[3*i + 2];
  }

  //for(auto t: triangles){
  //    std::cout << t << std::endl;
  //}
  return true;
}
