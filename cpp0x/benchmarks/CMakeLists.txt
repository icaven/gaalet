cmake_minimum_required (VERSION 2.6.3)
project (benchmarks)

INCLUDE_DIRECTORIES ("../../include/cpp0x")
set (CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_FLAGS "-std=c++0x")

file (GLOB benchmark_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")
foreach (benchmark_source ${benchmark_sources})
  string (REPLACE ".cpp" "" benchmark_prefix ${benchmark_source})
  add_executable (${benchmark_prefix} ${benchmark_source})
endforeach (benchmark_source)
