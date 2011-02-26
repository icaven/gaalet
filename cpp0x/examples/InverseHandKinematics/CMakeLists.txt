cmake_minimum_required (VERSION 2.6.3)
project (InverseHandKinematics)

include_directories ("../../../include/cpp0x")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

find_package(OpenSceneGraph REQUIRED osg osgViewer)
include_directories (${OPENSCENEGRAPH_INCLUDE_DIRS})
add_executable (InverseHandKinematics InverseHandKinematics.cpp)
target_link_libraries(InverseHandKinematics ${OPENSCENEGRAPH_LIBRARIES})

add_executable (InverseHandKinematicsScrew InverseHandKinematicsScrew.cpp)
target_link_libraries(InverseHandKinematicsScrew ${OPENSCENEGRAPH_LIBRARIES})

add_executable (InverseHandKinematicsLagrange InverseHandKinematicsLagrange.cpp)
target_link_libraries(InverseHandKinematicsLagrange ${OPENSCENEGRAPH_LIBRARIES})

add_executable (InverseHandKinematicsMetric InverseHandKinematicsMetric.cpp)
target_link_libraries(InverseHandKinematicsMetric ${OPENSCENEGRAPH_LIBRARIES})

add_executable (SymbolicEquationMetric SymbolicEquationMetric.cpp)
target_link_libraries(SymbolicEquationMetric ${OPENSCENEGRAPH_LIBRARIES})
