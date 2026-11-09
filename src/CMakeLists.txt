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

add_test(NAME MirrorMakeShape COMMAND MakeShape --icosahedron -e "-1 1 1" -i)
set_tests_properties(MirrorMakeShape PROPERTIES
    PASS_REGULAR_EXPRESSION "Volume: 2.53615"
)

add_test(NAME NegRadiusMakeShape COMMAND MakeShape --cube -r -1)
set_tests_properties(NegRadiusMakeShape PROPERTIES
    PASS_REGULAR_EXPRESSION "-r: R must be positive"
)

add_test(NAME ZeroAxisMakeShape COMMAND MakeShape --cube -a "0 0 0 30")
set_tests_properties(ZeroAxisMakeShape PROPERTIES
    PASS_REGULAR_EXPRESSION "-a: the rotation axis must be non-zero"
)

add_test(NAME EmptyCenterMakeShape COMMAND MakeShape -c)
set_tests_properties(EmptyCenterMakeShape PROPERTIES
    PASS_REGULAR_EXPRESSION "-c: cannot center a shape of zero volume"
)

add_test(NAME BadOutputMakeShape COMMAND MakeShape --cube -o /nonexistent_dir/x.off)
set_tests_properties(BadOutputMakeShape PROPERTIES
    PASS_REGULAR_EXPRESSION "Cannot open output file: /nonexistent_dir/x.off"
)

add_test(NAME FlatCenterMakeShape COMMAND MakeShape --cube -e "1 1 1e-17" -c)
set_tests_properties(FlatCenterMakeShape PROPERTIES
    PASS_REGULAR_EXPRESSION "-c: cannot center a shape of zero volume"
)

add_test(NAME FlatInfoMakeShape COMMAND MakeShape --cube -e "1 1 0" -i)
set_tests_properties(FlatInfoMakeShape PROPERTIES
    PASS_REGULAR_EXPRESSION "Center of mass: undefined \\(zero volume\\).*Radius from origin: 0.81649"
)

# Testing Shape2Zernike

add_test(NAME RunShape2Zernike COMMAND Shape2Zernike 0 ${CMAKE_SOURCE_DIR}/testdata/blork.off)

add_test(NAME ZeroApproxShape2Zernike COMMAND Shape2Zernike -a6 0 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
set_tests_properties(RunShape2Zernike ZeroApproxShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "ORTHO 0 COMPLEX\n0 0 0 0.052117"
)

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

add_test(NAME ZmLargeNShape2Zernike COMMAND Shape2Zernike 12 ${CMAKE_SOURCE_DIR}/testdata/blork.zm)
set_tests_properties(ZmLargeNShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "N larger than maximum moment available.*ORTHO 10 COMPLEX\n0 0 0 0.052117"
)

add_test(NAME DiffShape2Zernike COMMAND Shape2Zernike -n -p --diff ${CMAKE_SOURCE_DIR}/testdata/blork.zm 3 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
set_tests_properties(DiffShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "ORTHO_N 3 COMPLEX_P\n0 0 0 [0-9.]+e-1[0-9]\n"
)

add_test(NAME BadIndexShape2Zernike COMMAND Shape2Zernike 1 ${CMAKE_SOURCE_DIR}/testdata/badindex.off)
set_tests_properties(BadIndexShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "Cannot read file .*badindex.off at line 7"
)

add_test(NAME InwardShape2Zernike COMMAND Shape2Zernike 0 ${CMAKE_SOURCE_DIR}/testdata/blork_inward.off)
set_tests_properties(InwardShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "volume is negative.*0 0 0 -0.052117"
)

add_test(NAME ApproxWarnShape2Zernike COMMAND Shape2Zernike -a14 102 ${CMAKE_SOURCE_DIR}/testdata/flat.off)
set_tests_properties(ApproxWarnShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "Allowed error by facet: 1.0"
)
add_test(NAME ApproxNoWarnShape2Zernike COMMAND Shape2Zernike -a14 2 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
set_tests_properties(ApproxNoWarnShape2Zernike PROPERTIES
    FAIL_REGULAR_EXPRESSION "Allowed error by facet"
)

add_test(NAME PolygonShape2Zernike COMMAND Shape2Zernike -d10 0 ${CMAKE_SOURCE_DIR}/testdata/quadcube.off)
set_tests_properties(PolygonShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "Polygons: 6 facets.*0 0 0 0.1055381426"
)

add_test(NAME BadParityShape2Zernike COMMAND Shape2Zernike 3 ${CMAKE_SOURCE_DIR}/testdata/badparity.zm)
set_tests_properties(BadParityShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "Cannot read file .*badparity.zm at line 4"
)

add_test(NAME InvariantsShape2Zernike COMMAND Shape2Zernike -i -d10 6 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
add_test(NAME RotatedInvariantsShape2Zernike COMMAND Shape2Zernike -i -d10 6 ${CMAKE_SOURCE_DIR}/testdata/blork_rotated.off)
set_tests_properties(InvariantsShape2Zernike RotatedInvariantsShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "ZRI\nORTHO 6\n0 0 0 0.002716244362\n1 1 1 0.001230798227\n2 0 0 -0.003889811061"
)
add_test(NAME SignaturesShape2Zernike COMMAND Shape2Zernike -s -d10 6 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
add_test(NAME RotatedSignaturesShape2Zernike COMMAND Shape2Zernike -s -d10 6 ${CMAKE_SOURCE_DIR}/testdata/blork_rotated.off)
set_tests_properties(SignaturesShape2Zernike RotatedSignaturesShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "ZSI\nORTHO 6\n0 0.002716244362\n1 0.001230798227"
)

add_test(NAME NanShape2Zernike COMMAND Shape2Zernike 10 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
set_tests_properties(NanShape2Zernike PROPERTIES
    FAIL_REGULAR_EXPRESSION "nan;NAN;Nan"
)

add_test(NAME MissingShape2Zernike COMMAND Shape2Zernike 5 ${CMAKE_SOURCE_DIR}/testdata/blork.off)
set_tests_properties(MissingShape2Zernike PROPERTIES
    PASS_REGULAR_EXPRESSION "5 3 0.*5 3 1.*5 3 2.*5 3 3"
)

# Testing Zernike2Shape

add_test(NAME RunZernike2Shape COMMAND Zernike2Shape -t0 10 30 ${CMAKE_SOURCE_DIR}/testdata/blork.zm)
set_tests_properties(RunZernike2Shape PROPERTIES
    PASS_REGULAR_EXPRESSION "OFF.*213 422 0"
)

add_test(NAME BadResZernike2Shape COMMAND Zernike2Shape 10 0 ${CMAKE_SOURCE_DIR}/testdata/blork.zm)
set_tests_properties(BadResZernike2Shape PROPERTIES
    PASS_REGULAR_EXPRESSION "RES must be positive"
)
