cmake_minimum_required (VERSION 2.6.3)
project (ClosedLoopKinematics)

include_directories ("../../../include/cpp0x")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

# find_package(OpenSceneGraph REQUIRED osg osgViewer osgGA)
find_package(OpenSceneGraph REQUIRED osg osgViewer)
include_directories (${OPENSCENEGRAPH_INCLUDE_DIRS})
add_executable (ClosedLoopKinematics ClosedLoopKinematics.cpp)
set_property(TARGET ClosedLoopKinematics PROPERTY CXX_STANDARD 14)
target_link_libraries(ClosedLoopKinematics ${OPENSCENEGRAPH_LIBRARIES})
