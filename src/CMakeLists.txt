cmake_minimum_required(VERSION 3.0)

project(Shapes)

# Making Release build type the default to get optimization

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE
        STRING "Choose the type of build" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


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

# Add libraries

add_subdirectory(tools)
add_subdirectory(geom)
add_subdirectory(zernike)

# Add executables

add_executable(MakeShape MakeShape.cpp)
target_link_libraries(MakeShape geom)

add_executable(Zernike2Shape Zernike2Shape.cpp)
target_link_libraries(Zernike2Shape zernike)

add_executable(Shape2Zernike Shape2Zernike.cpp)
target_link_libraries(Shape2Zernike zernike)

# install executables

install(TARGETS MakeShape Zernike2Shape Shape2Zernike DESTINATION bin)