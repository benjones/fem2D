cmake_minimum_required(VERSION 2.6)
# Some flags for Freeglut and GLUI.
add_definitions(-D_CRT_SECURE_NO_WARNINGS )
SET(CMAKE_CXX_COMPILER clang++)
#ADD_DEFINITIONS(" -g -O0 ")
#ADD_DEFINITIONS("-std=c++11 -stdlib=libc++ -g  -O -Wall -Wuninitialized -Weverything -Wno-c++98-compat")
#ADD_DEFINITIONS("-std=c++11 -pg -stdlib=libc++ -static -O2")
ADD_DEFINITIONS("-std=c++11 -stdlib=libc++ -O2 -Wall -Wuninitialized")

# Define the framework files.
set(SRCS
  
  )


find_package(OpenGL REQUIRED)
#find_package(SDL2 REQUIRED)


include_directories (
  /usr/include
  /usr/local/include
  /usr/local/include/eigen3
  ${OPENGL_INCLUDE_DIR}
  /Library/Frameworks/SDL2.framework/Headers
  /System/Library/Frameworks/vecLib.framework/Headers
  )



#if(APPLE)
# We are not using the Apple's framework version, but X11's
#include_directories( /usr/include
#				/usr/local/include
#				/usr/X11/include )
#link_directories( /usr/X11/lib )
#set (OPENGL_LIBRARIES GL GLU GLUT X11)
#set (LINK_FLAGS "-framework OpenGL -framework GLUT")
#elseif(WIN32)
#	set (ADDITIONAL_LIBRARIES winmm)
#endif(APPLE)



add_executable(femApp
  femApp.cpp
  femMesh.cpp
  readObj.cpp
  readNodeAndEle.cpp
  range.cpp
  jsoncpp.cpp
  world.cpp
  planeObstacle.cpp
  )

#add_executable(dogRollerskatesGame
#	dogRollerskatesGame.cpp
#	${SRCS}
#)
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -framework SDL2 -framework vecLib" )

target_link_libraries ( 
  femApp
  ${ADDITIONAL_LIBRARIES}
  ${OPENGL_LIBRARIES}
  #mosek64
  c++
  stdc++
  ) 


