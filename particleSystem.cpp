
#include "readObj.h"
#include "range.h"
#include <iostream>
#include <set>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <thread>
#include <cassert>

#include <SDL.h>
#include <SDL_opengl.h>
//#include <OpenGL/gl.h>
#include <OpenGL/glu.h>

struct Edge{

  unsigned e1, e2; //indeces of the vertices on the edge
  unsigned t1, t2; //indeces of the points on the triangles incedent on the edge

  Edge(unsigned _e1, unsigned _e2, unsigned _t1): e1(_e1), e2(_e2), t1(_t1), t2(100000){}

  double edgeLength, dihedralLength;

};

template<typename P>
typename std::iterator_traits<P>::value_type  l2Norm(P p1, P p2){

  auto dx = *p1 - *p2;
  auto dy = *(p1 + 1) - *(p2 +1);
  auto dz = *(p1 + 2) - *(p2 + 2);
  
  auto ret =  sqrt(dx*dx + dy*dy + dz*dz);
  //  std::cout << "ret: " << ret << std::endl;
  return ret;
}


bool operator <(const Edge& lhs, const Edge& rhs){
  return lhs.e1 < rhs.e1 || (lhs.e1 == rhs.e1 && lhs.e2 < rhs.e2);
}






SDL_Window* surface;


void simulationLoop(std::vector<double>& vertices, std::vector<Edge>& edges, double dt);

int main(int argc, char**argv){


  std::vector<double> vertices, BBmin, BBmax;
  std::vector<unsigned> triangles;

  if(!readObj(argv[1], vertices, triangles, BBmin, BBmax)){
    std::cout << "couldn't read OBJ" << std::endl;
    exit(1);
  }

  std::cout << "bbmin: " << BBmin[0] << " " << BBmin[1] << " " << BBmin[2] << std::endl;
  std::cout << "bbmax: " << BBmax[0] << " " << BBmax[1] << " " << BBmax[2] << std::endl;

  for(auto i : range(vertices.size()/3)){
    vertices[3*i + 1] += 5;
  }

  std::cout << "vertices size: " << vertices.size() << std::endl;

  
  std::set<Edge> edgesSet;

  for(auto i : range(triangles.size()/3)){

    auto t = 3*i;

    for( auto j : range(3)){
      auto e1 = triangles[t + j];
      auto e2 = triangles[t + ((j +1)%3)];

      if( e2 < e1){
	std::swap(e1, e2);
      }

      auto e = Edge(e1, e2, triangles[t + ((j + 2)%3)]);
      auto it = edgesSet.find(e);

      if(it == edgesSet.end()){
	edgesSet.insert(it, e);
      } else {
	e.t2 = it->t1;
	edgesSet.erase(it);
	edgesSet.insert(e);
	
      }

    }
  }





  std::vector<Edge> edges(edgesSet.begin(), edgesSet.end());

  //  for(const auto& e : edges){
  // std::cout << "e1: " << e.e1 << " e2: " << e.e2 << " t1 " << e.t1 << " t2 : " << e.t2 << " e length " << e.edgeLength << " d length : " << e.dihedralLength << std::endl;
  // }

  std::cout << "about to set edge lengths" << std::endl;

  for(auto& e : edges){
    e.edgeLength = l2Norm(vertices.begin() + 3*e.e1,
			  vertices.begin() + 3*e.e2);
    
    if(e.t1 < vertices.size() && e.t2 < vertices.size()){
      e.dihedralLength = l2Norm(vertices.begin() + 3*e.t1,
				vertices.begin() + 3*e.t2);
    }
  }

  //for(auto i : range (vertices.size()/3)){
  //  std::cout << "vertex: " << i << " " << vertices[i*3] << " " << vertices[i*3 + 1] << " " << vertices[i*3 + 2] << std::endl;
  // }

  //for(const auto& e : edges){
  // std::cout << "e1: " << e.e1 << " e2: " << e.e2 << " t1 " << e.t1 << " t2 : " << e.t2 << " e length " << e.edgeLength << " d length : " << e.dihedralLength << std::endl;
  //}


  std::cout << "computed edge rest lengths" << std::endl;

  if(SDL_Init(SDL_INIT_EVERYTHING) < 0) {
    std::cout << "couldn't init SDL" << std::endl;
    exit(1);
  }
  
  if((surface = SDL_CreateWindow("3.2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,                
				 800, 600, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN)) == NULL){
    std::cout << "couldn't create SDL surface" << std::endl;
    exit(1);
  }

  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  auto* context = SDL_GL_CreateContext(surface);




  double dt = atof(argv[2]);
  
  std::cout << "about to loop" << std::endl;

  simulationLoop(vertices, edges, dt);
  

  return 0;
}


void drawFrame(std::vector<double>& vertices, std::vector<Edge>& edges){
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(20, 800/600.0, 0.1, 50); 

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 1, 10,
	    0, 1, 0,
	    0, 1, 0);

  glClearColor(0,0,0, 0);
  glClear(GL_COLOR_BUFFER_BIT);
  

  glColor3d(0.7, 0.7, 0.7);
  glBegin(GL_QUADS);
  glVertex3d(-100, 0, -100);
  glVertex3d(-100, 0, 100);
  glVertex3d(100, 0, 100);
  glVertex3d(100, 0, -100);
  glEnd();


  glPointSize(2);
  glColor3d(1,0,1);
  glBegin(GL_POINTS);
  for(auto i : range(vertices.size()/3)){
    glVertex3d(vertices[3*i],
	       vertices[3*i +1],
	       vertices[3*i +2]);
    
  }
  glEnd();

  glColor3d(0,1,0);
  glBegin(GL_LINES);
  for(auto& e : edges){
    
    
    glVertex3d(vertices[3*e.e1],
	       vertices[3*e.e1+1],
	       vertices[3*e.e1+2]);

    glVertex3d(vertices[3*e.e2],
	       vertices[3*e.e2+1],
	       vertices[3*e.e2+2]);


  }
  glEnd();

  glFlush();

  SDL_GL_SwapWindow(surface);
}


std::vector<double> forces;
void simulateFrame(std::vector<double>& vertices, std::vector<double>& velocities, std::vector<Edge>& edges, double dt){

  
  //gravity
  //assume everything has mass 1

  forces.assign(vertices.size(), 0.0);

  for(auto i : range(vertices.size()/3)){
    forces[3*i +1] -= 9.81;
  }

  double kStretch = 20000;
  double kDamp = 100;
  double kBend = 20000;

  for(auto& e : edges){
    
    double direction[3];
    auto currentDistance = l2Norm(vertices.begin() + 3*e.e1, vertices.begin() + 3*e.e2);
    
    direction[0] = vertices[3*e.e1] - vertices[3*e.e2];
    direction[1] = vertices[3*e.e1 + 1] - vertices[3*e.e2 + 1];
    direction[2] = vertices[3*e.e1 + 2] - vertices[3*e.e2 + 2];
    
    direction[0] /= currentDistance;
    direction[1] /= currentDistance;
    direction[2] /= currentDistance;

    auto magnitude = kStretch*(e.edgeLength - currentDistance);
    forces[3*e.e1] += direction[0]*magnitude;
    forces[3*e.e1 + 1] += direction[1]*magnitude;
    forces[3*e.e1 + 2] += direction[2]*magnitude;

    forces[3*e.e2] -= direction[0]*magnitude;
    forces[3*e.e2 + 1] -= direction[1]*magnitude;
    forces[3*e.e2 + 2] -= direction[2]*magnitude;


    //damping forces:
    
    double deltaV[3];
    deltaV[0] = velocities[3*e.e1] - velocities[3*e.e2];
    deltaV[1] = velocities[3*e.e1 + 1] - velocities[3*e.e2 + 1];
    deltaV[2] = velocities[3*e.e1 + 2] - velocities[3*e.e2 + 2];
    
    magnitude = -kDamp*(deltaV[0]*direction[0] + deltaV[1]*direction[1] + deltaV[2]*direction[2]);
    
    //forces[3*e.e1] += deltaV[0]*directio 

    //magnitude = -kDamp*velocityMag;
    forces[3*e.e1] += direction[0]*magnitude;
    forces[3*e.e1 + 1] += direction[1]*magnitude;
    forces[3*e.e1 + 2] += direction[2]*magnitude;
    
    forces[3*e.e2] -= direction[0]*magnitude;
    forces[3*e.e2 + 1] -= direction[1]*magnitude;
    forces[3*e.e2 + 2] -= direction[2]*magnitude;
    

    //bending forces:
    if(e.t1 < (vertices.size()/3) && e.t2 < (vertices.size()/3)){
      
      currentDistance = l2Norm(vertices.begin() + 3*e.t1, vertices.begin() + 3*e.t2);
    
      direction[0] = vertices[3*e.t1] - vertices[3*e.t2];
      direction[1] = vertices[3*e.t1 + 1] - vertices[3*e.t2 + 1];
      direction[2] = vertices[3*e.t1 + 2] - vertices[3*e.t2 + 2];
    
      direction[0] /= currentDistance;
      direction[1] /= currentDistance;
      direction[2] /= currentDistance;
      
      auto magnitude = kBend*(e.dihedralLength - currentDistance);
      forces[3*e.t1] += direction[0]*magnitude;
      forces[3*e.t1 + 1] += direction[1]*magnitude;
      forces[3*e.t1 + 2] += direction[2]*magnitude;
      
      forces[3*e.t2] -= direction[0]*magnitude;
      forces[3*e.t2 + 1] -= direction[1]*magnitude;
      forces[3*e.t2 + 2] -= direction[2]*magnitude;


    }
    
    

  }

  



  //timestep velocities
  for(auto i : range(velocities.size())){
    velocities[i] += forces[i]*dt;
  }
  //timestep positions
  for(auto i : range(vertices.size())){
    vertices[i] += velocities[i]*dt;
  }
  for(auto i : range(vertices.size()/3)){
    if(vertices[3*i +1] < 0){
      vertices[3*i +1] = 0;
      velocities[3*i + 1] = 0;
    }
  }

}

void simulationLoop(std::vector<double>& vertices, std::vector<Edge>& edges, double dt){


  std::vector<double> velocities(vertices.size(), 0);


  std::chrono::time_point<std::chrono::high_resolution_clock> start, end, renderEnd;
  //size_t start, end, renderEnd;
  renderEnd = std::chrono::high_resolution_clock::now();
  while(true){
  
    SDL_Event event;
    
    while(SDL_PollEvent(&event)){

      switch(event.type){
      case SDL_KEYDOWN:
	exit(0);
      default:
	0;// do nothing
      }
    }

    start = std::chrono::high_resolution_clock::now();
    
    simulateFrame(vertices, velocities, edges, dt);
    end = std::chrono::high_resolution_clock::now();
    
    //std::chrono::microseconds elapsed_u_seconds = end- renderEnd;//start;
    
    auto microseconds_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - renderEnd);
    auto dtInMicroseconds = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::duration<double>(dt));

    //std::cout << "elapsed: " << microseconds_elapsed.count() << "  dt: "  << dtInMicroseconds.count() << std::endl;;

    /*if(elapsed_seconds.count() < dt){
      std::this_thread::sleep_for(std::chrono::duration<double>(dt - elapsed_seconds.count()
								));
      std::cout << "sim took: " << 1.0/elapsed_seconds.count() << " hz" << std::endl;
      std::cout << "slept for : " << std::chrono::duration<double>(dt - elapsed_seconds.count()).count() << std::endl;
      } */
    //while( (std::chrono::high_resolution_clock::now() - renderEnd).count() < dt){
      //spin
    //}
    //std::cout << "waited " << (std::chrono::high_resolution_clock::now() - end).count() << " seconds" << std::endl;

    /*if((end - renderEnd) < (1000*dt)){
      auto toWait = (1000*dt) - (end - renderEnd);
      std::cout << "waiting " << toWait << "ms" << std::endl;
      SDL_Delay(toWait);

      }*/

    while(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - renderEnd).count() <
	  dtInMicroseconds.count()){
      //std::this_thread::sleep_for(dtInMicroseconds - microseconds_elapsed);
    }
    //std::cout << "simulation took " << microseconds_elapsed.count() << "us" << std::endl;
    //std::cout << "waited: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - end).count() << "us" << std::endl;
    drawFrame(vertices, edges);


    renderEnd = std::chrono::high_resolution_clock::now();

    
    
  }
}
