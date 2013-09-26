#pragma once

#include <vector>
#include <string>

bool readNodeAndEle(std::string filename,
		    std::vector<double>& vertices,
		    std::vector<unsigned>& triangles);
