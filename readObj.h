

#include <vector>

bool readObj(const char* filename,
	     std::vector<double>& vertices,
	     std::vector<unsigned>& triangles,
	     std::vector<double>& BBmin,
	     std::vector<double>& BBmax);
