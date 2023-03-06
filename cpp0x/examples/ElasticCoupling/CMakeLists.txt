cmake_minimum_required (VERSION 2.8.12)
project (ElasticCoupling)

include_directories ("../../../include/cpp0x")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

find_package(OpenSceneGraph REQUIRED osg osgViewer osgGA)
include_directories (${OPENSCENEGRAPH_INCLUDE_DIRS})
add_executable (ElasticCoupling ElasticCoupling.cpp)
set_property(TARGET ElasticCoupling PROPERTY CXX_STANDARD 14)
target_link_libraries(ElasticCoupling ${OPENSCENEGRAPH_LIBRARIES})
