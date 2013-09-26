#include "femMesh.h"
#include "range.cpp"
#include "utils.h"
#include "planeObstacle.h"
#include <SDL_opengl.h>


void FemMesh::initialize(const std::vector<double>& vertices,
			 const std::vector<unsigned>& triangles,
			 MaterialParameters _materialParameters){
  materialParameters = _materialParameters;

  nodes.resize(vertices.size()/2);
  for(auto i : range(nodes.size())){
    nodes[i].worldPos(0) = vertices[2*i];
    nodes[i].worldPos(1) = vertices[2*i+1];

    nodes[i].referencePos = nodes[i].worldPos;
    nodes[i].velocity.setZero();
    nodes[i].mass = 0;

  }

  elements.resize(triangles.size()/3);
  for(auto i : range(elements.size())){
    elements[i].nodes[0] = triangles[3*i];
    elements[i].nodes[1] = triangles[3*i + 1];
    elements[i].nodes[2] = triangles[3*i + 2];
  }

  //compute mass and beta matrices for each element

  for(auto &e : elements){

    auto& u0 = nodes[e.nodes[0]].referencePos;
    auto& u1 = nodes[e.nodes[1]].referencePos;
    auto& u2 = nodes[e.nodes[2]].referencePos;
    
    e.beta.col(0) = u1 - u0;
    e.beta.col(1) = u2 - u0;

    //area (length) weighted face normals
    e.faceNormals[0] = (u1 - u0).norm()*vec2(u1.y() - u0.y(), u0.x() - u1.x());
    e.faceNormals[1] = (u2 - u1).norm()*vec2(u2.y() - u1.y(), u1.x() - u2.x());
    e.faceNormals[2] = (u0 - u2).norm()*vec2(u0.y() - u2.y(), u2.x() - u0.x());

    
    double area = 0.5*(e.beta.col(0)(0)*e.beta.col(1)(1) - 
		       e.beta.col(0)(1)*e.beta.col(1)(0));
    
    for(auto i : range(3)){
      nodes[e.nodes[i]].mass += area*materialParameters.density;
    }
    
    e.beta = e.beta.inverse().eval();
  }

  for(auto & n : nodes){
    n.mass /= 3;
  }

}



void FemMesh::renderOpenGL(){

  
  glColor3f(0,1,0);
  glBegin(GL_LINES);
  for(auto& e : elements){
    for(auto i : range(3)){
      
      glVertex2d(nodes[e.nodes[i]].worldPos(0),
		 nodes[e.nodes[i]].worldPos(1));
      glVertex2d(nodes[e.nodes[(i+1)%3]].worldPos(0),
		 nodes[e.nodes[(i+1)%3]].worldPos(1));
    }
  }
  glEnd();

}


void FemMesh::integrate(double dt){
  //clear forces:
  for(auto& n : nodes){
    if(true){ //gravity
      n.force = vec2(0, -9.81)* n.mass;
    } else {
      n.force.setZero();
    }
  }
  
  //compute elastic/damping forces
  for(auto & e : elements){

    auto& x0 = nodes[e.nodes[0]].worldPos;
    auto& x1 = nodes[e.nodes[1]].worldPos;
    auto& x2 = nodes[e.nodes[2]].worldPos;
    
    auto& v0 = nodes[e.nodes[0]].velocity;
    auto& v1 = nodes[e.nodes[1]].velocity;
    auto& v2 = nodes[e.nodes[2]].velocity;

    mat2 xMatrix;
    xMatrix.col(0) = x1 - x0;
    xMatrix.col(1) = x2 - x0;

    mat2 vMatrix;
    vMatrix.col(0) = v1 - v0;
    vMatrix.col(1) = v2 - v0;

    auto deformationGradient = xMatrix*e.beta;
    auto velocityGradient = vMatrix*e.beta;

    mat2 R, S;
    Utils::RSDecomp(deformationGradient, R, S);
    
    auto fTwiddle = R.transpose()*deformationGradient;
    auto vTwiddle = R.transpose()*velocityGradient;
    
    auto strain = 0.5*(fTwiddle + fTwiddle.transpose()) - mat2::Identity();
    auto strainRate = 0.5*(vTwiddle +
			   vTwiddle.transpose());

    auto stress = materialParameters.lambda*strain.trace()*mat2::Identity() +
      2*materialParameters.mu*strain;
    
    auto viscousStress = 
      materialParameters.dampLambda*strainRate.trace()*mat2::Identity() +
      2*materialParameters.dampMu*strainRate;

    auto forceProduct = R*(stress + viscousStress);
    for(auto i : range(3)){
      nodes[e.nodes[i]].force += forceProduct*e.faceNormals[(i+1)%3];
    }
  }
  
  //compute deformation gradient

  //forward symplectic euler
  //update velocities
  for(auto& n : nodes){
    n.velocity += n.force*(dt/n.mass);
    n.worldPos += n.velocity*dt;
  }
  
  //fix 2 vertices of the triangle for testing
  for(auto i : {2,3}){
    nodes[i].worldPos = nodes[i].referencePos;
    nodes[i].velocity.setZero();
  }
}

void FemMesh::collidePlanes(const std::vector<PlaneObstacle>& planes){

  for(const auto& p : planes){
    for(auto& n : nodes){
      
      if(n.worldPos.dot(p.normal) < p.offset){

	//BAD BAD BAD
	n.worldPos.y() *= -1;
	n.velocity.y() *= -1;
      }
    }
  }

}
