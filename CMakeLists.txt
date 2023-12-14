cmake_minimum_required(VERSION 3.15)
project(gaalet)

set(CMAKE_CXX_STANDARD 14)

# --- Use CPM - CMake's missing package manager (https://github.com/cpm-cmake/CPM.cmake/tree/master#adding-cpm)
# download CPM.cmake
file(
        DOWNLOAD
        https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.38.3/CPM.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
        EXPECTED_HASH SHA256=cc155ce02e7945e7b8967ddfaff0b050e958a723ef7aad3766d368940cb15494
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

# -- Eigen is only used for vector and matrix classes in the test main.cpp
CPMAddPackage(
        NAME Eigen
        VERSION 3.4.0
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        # Eigen's CMakelists are not intended for library use
        DOWNLOAD_ONLY YES
)

include_directories(
        include/cpp0x
)

add_subdirectory(cpp0x/tests)
add_subdirectory(cpp0x/examples)
add_subdirectory(cpp0x/benchmarks)
add_subdirectory(python)

add_executable(${PROJECT_NAME} main.cpp)
if(Eigen_ADDED)
    target_include_directories(${PROJECT_NAME} PUBLIC ${Eigen_SOURCE_DIR})
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE
        ${LINK_OPTIONS}
        )

