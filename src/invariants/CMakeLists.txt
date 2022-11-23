#Written by J. Houdayer

add_library(invariants invariants.cpp coefs.cpp sphericalbessel.cpp s_root_data.cpp k_distance.cpp)
target_include_directories(invariants INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(invariants zernike GMP::GMPXX)