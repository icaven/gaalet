cmake_minimum_required(VERSION 3.15)
project(gaalet)

set(CMAKE_CXX_STANDARD 14)

# Use pkgconfig to find the paths for OpenCV 4.x (OpenCV has been configured with OPENCV_GENERATE_PKGCONFIG=1)
FIND_PACKAGE(PkgConfig)

IF(PKG_CONFIG_FOUND)
    PKG_CHECK_MODULES(PC_OPENCV IMPORTED_TARGET opencv4)
ENDIF()

include_directories(
        include/cpp0x
        ${PC_OPENCV_INCLUDE_DIRS}
        )

add_subdirectory(cpp0x/tests)
add_subdirectory(cpp0x/examples)
add_subdirectory(cpp0x/benchmarks)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME}
        ${LINK_OPTIONS}
        ${PC_OPENCV_LINK_LIBRARIES}
        )
