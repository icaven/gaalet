cmake_minimum_required (VERSION 2.6.3)
project (PointsAndLines)

include_directories ("../../../../include/cpp0x")
include_directories ("../")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

find_package(OpenSceneGraph REQUIRED osg osgViewer osgGA)
include_directories (${OPENSCENEGRAPH_INCLUDE_DIRS})
add_executable (PointsAndLines PointsAndLines.cpp)
set_property(TARGET PointsAndLines PROPERTY CXX_STANDARD 14)
target_link_libraries(PointsAndLines ${OPENSCENEGRAPH_LIBRARIES})
