cmake_minimum_required(VERSION 3.1)

project(Shapes)


# Global settings

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -Wfatal-errors -pedantic -O3)


# Threads

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads)


# Compute version number (from "git describe")

find_package(Git QUIET)
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

    configure_file(version.hpp.in version.hpp)
    include_directories("${PROJECT_BINARY_DIR}")
else()
    message("You need git installed to produce version numbers")
endif()

# Documentation (using Doxygen)

find_package(Doxygen QUIET)
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
link_libraries(tools)

add_subdirectory(geom)
link_libraries(geom)

add_subdirectory(zernike)


# Add executables

add_executable(MakeShape MakeShape.cpp)

link_libraries(zernike)

add_executable(Zernike2Shape Zernike2Shape.cpp)
add_executable(Shape2Zernike Shape2Zernike.cpp)
