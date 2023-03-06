cmake_minimum_required (VERSION 2.8.12)
project (MonteCarloPi)

include_directories ("../../../include/cpp0x")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

add_executable (MonteCarloPi MonteCarloPi.cpp)
set_property(TARGET MonteCarloPi PROPERTY CXX_STANDARD 14)
