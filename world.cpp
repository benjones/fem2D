#include "world.h"

#include <fstream>
#include "json/json.h"
#include "range.h"
#include "readNodeAndEle.h"


World::World(std::string filename){

  std::ifstream ins(filename.c_str());

  Json::Value root;
  Json::Reader jReader;
  
  bool jsonSuccess = jReader.parse(ins, root);
  ins.close();
  if(!jsonSuccess){
    std::cout << "couldn't read input file: " << filename << std::endl
	      << jReader.getFormattedErrorMessages();
    exit(1);
  }
  auto femObjectsIn = root["femObjects"];
  for(auto i = 0; i < femObjectsIn.size(); ++i){
    auto femJson = femObjectsIn[i];
    femObjects.emplace_back();
    auto& f = femObjects.back();
    std::string meshName = femJson["meshName"].asString();
    std::cout << "added FEM object from file: " << femJson["meshName"].asString() << std::endl;
    
    std::vector<double> vertices;
    std::vector<unsigned> triangles;
    
    readNodeAndEle(meshName, vertices, triangles);
    
    auto materialParameters = MaterialParameters 
      { femJson.get("density", 1000).asDouble(),
	femJson.get("lambda", 1e4).asDouble(),
	femJson.get("mu", 1e4).asDouble(),
	femJson.get("dampLambda", 1).asDouble(),
	femJson.get("dampMu", 1).asDouble()};

    if( !femJson["centerTo"].isNull()){

      vec2 COM = vec2::Zero();
      for(auto i : range(vertices.size()/2)){
	COM.x() += vertices[2*i];
	COM.y() += vertices[2*i +1];
      }
      COM /= (vertices.size()/2);
      auto delta = vec2(femJson["centerTo"]["x"].asDouble(),
			femJson["centerTo"]["y"].asDouble()) -
	COM;

      for(auto i : range(vertices.size()/2)){
	vertices[2*i] += delta.x();
	vertices[2*i +1] += delta.y();
      }
    }
    
    if( !femJson["desiredScale"].isNull()){
      auto bbMin = vec2(1e10, 1e10);
      auto bbMax = vec2(-1e10, -1e10);

      for(auto i : range(vertices.size()/2)){
	bbMin(0) = std::min(bbMin(0), vertices[2*i]);
	bbMin(1) = std::min(bbMin(0), vertices[2*i + 1]);
	bbMax(0) = std::max(bbMin(0), vertices[2*i]);
	bbMax(1) = std::max(bbMin(0), vertices[2*i + 1]);
      }

      auto scale = femJson["desiredScale"].asDouble()/(std::max(bbMax(0) - bbMin(0),
							   bbMax(1) - bbMin(1)));
      for(auto& v : vertices){
	v*= scale;
      }
    }


    f.initialize(vertices, triangles, materialParameters);
    
  }

  auto planeObstaclesIn = root["planeObstacles"];
  for(auto i = 0; i < planeObstaclesIn.size(); ++i){
    auto femPlane = planeObstaclesIn[i];
    planeObstacles.emplace_back(vec2(femPlane["x"].asDouble(),
				     femPlane["y"].asDouble()),
				femPlane["offset"].asDouble());
    std::cout << "added: " << planeObstacles.back() << std::endl;
  }

  dt = root.get("dt", 0.001).asDouble();

}


void World::integrate(){


  for(auto & f : femObjects){
    f.integrate(dt);
    f.collidePlanes(planeObstacles);
  }
}


void World::renderOpenGL(){

  for(auto & f : femObjects){
    f.renderOpenGL();
  }
  for(auto & p : planeObstacles){
    p.renderOpenGL();
  }
}
