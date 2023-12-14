cmake_minimum_required(VERSION 3.22)

# --- Project

project(
        clifford_gaalet_example
        VERSION 1.0.0
        LANGUAGES CXX
)

# ---- Include guards

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(
            FATAL_ERROR
            "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
    )
endif()

# --- Use CPM - CMake's missing package manager (https://github.com/cpm-cmake/CPM.cmake.git)

include(${CMAKE_BINARY_DIR}/cmake/CPM.cmake)

# --- Add dependencies

# Define the Python3_ROOT_DIR variable to python3 interpreter grandparent directory - set this as a define to cmake
if(NOT DEFINED Python3_ROOT_DIR)
    message("The Python interpreter hasn't been specified so it will be searched for...")
    set (Python3_FIND_VIRTUALENV FIRST)
else(NOT DEFINED Python3_ROOT_DIR)
    message("Using Python3_ROOT_DIR: " ${Python3_ROOT_DIR})
endif(NOT DEFINED Python3_ROOT_DIR)
set(Python3_FIND_STRATEGY LOCATION)
find_package (Python3 COMPONENTS Interpreter Development NumPy)

if(NOT Python3_FOUND)
    message("Failed to find the Python3 interpreter")
else()
    message("Python3: " ${Python3_EXECUTABLE})
endif(NOT Python3_FOUND)

# The Python3_SOABI variable doesn't have the extension on the library name, so get the extension suffix with it
execute_process(COMMAND python3-config --extension-suffix
        RESULT_VARIABLE extension_suffix_result
        OUTPUT_VARIABLE Python3_EXTENSION_SUFFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

CPMAddPackage("gh:pybind/pybind11@2.9.1")

# pybind11-stubgen is used to generate python stubs for the extension module
CPMAddPackage(
        NAME pybind11-stubgen
        GIT_REPOSITORY https://github.com/sizmailov/pybind11-stubgen.git
        GIT_TAG 2536dbc3cad276c9edea86a3af61de9a57ab644d
        DOWNLOAD_ONLY YES
)

set(stubgen ${pybind11-stubgen_SOURCE_DIR}/pybind11_stubgen/__init__.py)


# --- Create library

# If not debugging, remove debugging statements
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    add_definitions(-DNDEBUG)                   # Disable asserts
endif ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")

file(GLOB sources CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.h
        )
pybind11_add_module(${PROJECT_NAME} ${sources})

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

target_include_directories(
        ${PROJECT_NAME} PUBLIC
        ../include/cpp0x
)

target_link_directories(
        ${PROJECT_NAME} PUBLIC
)

target_link_libraries(
        ${PROJECT_NAME} PUBLIC
)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
target_compile_definitions(${PROJECT_NAME} PRIVATE
        -DMODULE_VERSION="${PROJECT_VERSION}"
        )

set_target_properties(${PROJECT_NAME} PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION ON
        CXX_VISIBILITY_PRESET hidden
        VISIBLITY_INLINES_HIDDEN ON
        )

get_target_property(TARGET_NAME ${PROJECT_NAME} NAME)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${TARGET_NAME} )

# Run a command to generate the python stubs
# The library must be installed before these stubs can be properly generated
# TODO: Add a dependency to the installed library
add_custom_target(${TARGET_NAME}_interface_file
        COMMAND "${Python3_EXECUTABLE}" ${stubgen} --no-setup-py ${TARGET_NAME} -o ${CMAKE_CURRENT_BINARY_DIR}
        BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}-stubs/__init__.pyi"
        )
add_dependencies(${TARGET_NAME}_interface_file ${TARGET_NAME})

# Create a bootstrap file for the library
file (
        WRITE ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.py
        "def __bootstrap__():
    global __bootstrap__, __loader__, __file__
    import sys, pkg_resources, imp
    __file__ = pkg_resources.resource_filename(__name__, '${TARGET_NAME}${Python3_EXTENSION_SUFFIX}')
    __loader__ = None; del __bootstrap__, __loader__
    imp.load_dynamic(__name__,__file__)
__bootstrap__()")

# Install a .pth file for the package
file (
        WRITE ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pth
        "./${TARGET_NAME}
        ")


# Install in the site-packages directory
set(CMAKE_INSTALL_PREFIX ${Python3_SITEARCH} CACHE PATH "..." FORCE)
install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${Python3_SITEARCH}/${TARGET_NAME})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.py DESTINATION ${Python3_SITEARCH}/${TARGET_NAME})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pth DESTINATION ${Python3_SITEARCH})
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}-stubs/__init__.pyi" DESTINATION ${Python3_SITEARCH}/${TARGET_NAME}/)

