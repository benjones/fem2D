
#include "readObj.h"
#include <fstream>
#include <sstream>

bool readObj(const char* filename,
	     std::vector<double>& vertices,
	     std::vector<unsigned>& triangles,
	     std::vector<double>& BBmin,
	     std::vector<double>& BBmax){

  
  vertices.clear();
  triangles.clear();
  BBmin.resize(3);  
  BBmax.resize(3);

  std::ifstream in(filename);

  char line[255];

  BBmin[0] = BBmin[1] = BBmin[2] = 1000000;
  BBmax[0] = BBmax[1] = BBmax[2] = -1000000;

  
  std::string word;

  while(in.good()){
    in.getline(line, 254);
    if(line[0] == 'v'){
      double x, y, z;
      std::istringstream iss(line);
      iss >> word >> x>> y >> z;
      if(x < BBmin[0])
	BBmin[0] = x;
      if(y < BBmin[1])
	BBmin[1] = y;
      if(z < BBmin[2])
	BBmin[2] = z;
      if(x > BBmax[0])
	BBmax[0] = x;
      if(y > BBmax[1])
	BBmax[1] = y;
      if(z > BBmax[2])
	BBmax[2] = z;


      vertices.push_back(x);
      vertices.push_back(y);
      vertices.push_back(z);
    } else if(line[0] == 'f'){
      double i, j, k;
      std::istringstream iss(line);
      iss >> word >> i >> j >> k;
      triangles.push_back(i -1);
      triangles.push_back(j -1);
      triangles.push_back(k -1);
      
    }

  }
  in.close();
  
  return true;
}
