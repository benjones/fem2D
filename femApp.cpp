//#include "readNodeAndEle.h"
//#include "femMesh.h"
#include "world.h"

#include <SDL.h>
#include <SDL_opengl.h>
#include <OpenGL/glu.h>

#include <iostream>

#include "utils.h"


SDL_Window* surface;

void simulationLoop(World& mesh);


int main(int argc, char** argv){


  //std::vector<double> vertices;
  //std::vector<unsigned> triangles;
  
  //readNodeAndEle(argv[1], vertices, triangles);

  //FemMesh mesh;
  //mesh.initialize(vertices, triangles, 1);

  auto world = World(argv[1]);

  std::cout << "initialized world" << std::endl;

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
  
  std::cout << "about to loop" << std::endl;

  simulationLoop(world);
  

  


  return 0;
}



void drawFrame(World& world){

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-2, 2, -2, 2);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();


  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);

  world.renderOpenGL();

  glFlush();
  SDL_GL_SwapWindow(surface);
}


void simulateFrame(World& world){

  world.integrate();
  
}

void simulationLoop(World& world){

  while(true){
    SDL_Event event;
    
    while(SDL_PollEvent(&event)){
      
      switch(event.type){
      case SDL_KEYDOWN:
	exit(0);
	//default:
	//0;// do nothing
      }
    }
    simulateFrame(world);
    drawFrame(world);
  }

}
