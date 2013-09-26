#include "planeObstacle.h"
#include <SDL_opengl.h>

void PlaneObstacle::renderOpenGL(){


  vec2 tangent(normal.y(), -normal.x());

  vec2 pointOnLine = normal.x() > 0.001 ? vec2(offset/normal.x(), 0) : vec2(0,1);



  glColor3f(0,0.3, 0.8);
  glBegin(GL_LINES);
  glVertex2f(pointOnLine.x() - 100*tangent.x(), pointOnLine.x() - 100*tangent.y());
  glVertex2f(pointOnLine.x() + 100*tangent.x(), pointOnLine.x() + 100*tangent.y());
  glEnd();
}
