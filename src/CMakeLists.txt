# Written by J. Houdayer

cmake_minimum_required(VERSION 3.11)

project(Shapes)

# Making Release build type the default to get optimization

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE
        STRING "Choose the type of build" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#let cmake know about the cmake directory

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


# Compilation settings

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_11)

set(gcc_cxx "$<CXX_COMPILER_ID:ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<CXX_COMPILER_ID:MSVC>")
target_compile_options(compiler_flags INTERFACE
    "$<${gcc_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wfatal-errors;-pedantic;-Wunused;$<$<CONFIG:Debug>:-ggdb>>>"
    "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)

link_libraries(compiler_flags)


# Threads

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads)

if (Threads_FOUND)
    option(USE_THREADS "Compile with threads" ON)
endif()

if (NOT USE_THREADS)
    add_definitions(-DNO_THREADS)
endif()


# Compute version number (from "git describe")

find_package(Git)
if (GIT_FOUND)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE version
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    set_property(GLOBAL APPEND
        PROPERTY CMAKE_CONFIGURE_DEPENDS
        "${CMAKE_SOURCE_DIR}/../.git/index"
        )
else()
    set(version "version unknown")
    message("You need git installed to produce version numbers")
endif()
        
include_directories("${PROJECT_BINARY_DIR}")
configure_file(version.hpp.in version.hpp)


# Documentation (using Doxygen)

find_package(Doxygen)
if (DOXYGEN_FOUND)
    configure_file(Doxyfile.in Doxyfile @ONLY)
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        BYPRODUCTS doc
        )
else()
    message("You need Doxygen installed to produce documentation")
endif()

# GMP

find_package(GMP)


# Add libraries

add_subdirectory(tools)
add_subdirectory(geom)
add_subdirectory(zernike)

if (GMP_FOUND)
    add_subdirectory(invariants)
endif()

# Add executables

add_executable(MakeShape MakeShape.cpp)
target_link_libraries(MakeShape geom)

add_executable(Zernike2Shape Zernike2Shape.cpp)
target_link_libraries(Zernike2Shape zernike)

add_executable(Shape2Zernike Shape2Zernike.cpp)
target_link_libraries(Shape2Zernike zernike)

add_executable(Shape2Invariant Shape2Invariant.cpp)
target_link_libraries(Shape2Invariant invariants)

# install executables

install(TARGETS MakeShape Zernike2Shape Shape2Zernike DESTINATION bin)

# testing

enable_testing()

# Testing MakeShape

add_test(NAME RunMakeShape COMMAND MakeShape)

add_test(NAME BlorkMakeShape COMMAND MakeShape -l ${CMAKE_SOURCE_DIR}/testdata/blork.off -i)
set_tests_properties(BlorkMakeShape PROPERTIES
    PASS_REGULAR_EXPRESSION "Center of mass: 0.2 0 0.225.*Radius from center of mass: 1.2412.*Area: 2.0061.*Volume: 0.106667"
)

add_test(NAME ManyMakeShape COMMAND MakeShape --sphere -s2 -r2 -t "1 2 3" -i)
set_tests_properties(ManyMakeShape PROPERTIES
    PASS_REGULAR_EXPRESSION "Center of mass: 1 2 3.*Radius from center of mass: 2.*Area: 49.31.*Volume: 32.37"
)

# Testing Shape2Zernike

add_test(NAME RunShape2Zernike COMMAND Shape2Zernike 0 ${CMAKE_SOURCE_DIR}/testdata/blork.off)

add_test(NAME BlorkShape2Zernike COMMAND Shape2Zernike 5 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
add_test(NAME BlorkApproxShape2Zernike COMMAND Shape2Zernike -a6 5 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
set_tests_properties(BlorkShape2Zernike BlorkApproxShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "# Mesh: 4 vertices, 4 facets, radius: 1.*0 0 0 0.052117.*3 3 2 0.00313953.*0.000457147"
)

add_test(NAME CubeShape2Zernike COMMAND Shape2Zernike -t0 -rd12 20 ${CMAKE_SOURCE_DIR}/testdata/cube.off)
add_test(NAME CubeApproxShape2Zernike COMMAND Shape2Zernike -t0 -rd12 -a10 20 ${CMAKE_SOURCE_DIR}/testdata/cube.off)
set_tests_properties(CubeShape2Zernike CubeApproxShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "0 0 0 0.7522527780.*10 4 0 -0.0679492421.*20 16 12 0.0046415317"
)

add_test(NAME NanShape2Zernike COMMAND Shape2Zernike 10 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
set_tests_properties(NanShape2Zernike PROPERTIES
    FAIL_REGULAR_EXPRESSION "nan;NAN;Nan"
)

# Testing Zernike2Shape

add_test(NAME RunZernike2Shape COMMAND Zernike2Shape -t0 10 30 ${CMAKE_SOURCE_DIR}/testdata/blork.zm)
set_tests_properties(RunZernike2Shape PROPERTIES
    PASS_REGULAR_EXPRESSION "OFF.*553 422 0"
)

