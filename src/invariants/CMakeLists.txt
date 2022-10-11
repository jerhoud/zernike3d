#Written by J. Houdayer

add_library(invariants invariants.cpp)
target_include_directories(invariants INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(invariants zernike)