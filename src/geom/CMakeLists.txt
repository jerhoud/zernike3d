# Written by J. Houdayer

add_library(geom vec.cpp triangle.cpp mesh.cpp neldermead.cpp)
target_include_directories(geom INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(geom tools)
