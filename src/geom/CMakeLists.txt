add_library(geom vec.cpp triangle.cpp mesh.cpp)
target_include_directories(geom INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
