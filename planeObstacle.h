#pragma once

#include "eigenTypedefs.h"
#include <ostream>

class PlaneObstacle{
public:
  vec2 normal;
  double offset;
  
  PlaneObstacle(vec2 _normal, double _offset) : normal(_normal), offset(_offset) {}

  void renderOpenGL();

};


inline std::ostream& operator<<(std::ostream& outs, const PlaneObstacle& obs){
  outs << "x: " << obs.normal.x() << " y: " << obs.normal.y() << " offset: " << obs.offset;
  return outs;
}
