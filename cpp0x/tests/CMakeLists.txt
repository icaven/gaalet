cmake_minimum_required (VERSION 2.6.3)
project (tests)

include_directories ("../../include/cpp0x")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif (NOT CMAKE_BUILD_TYPE)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

file (GLOB test_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")
foreach (test_source ${test_sources})
  string (REPLACE ".cpp" "" test_prefix ${test_source})
  add_executable (${test_prefix} ${test_source})
endforeach (test_source)
