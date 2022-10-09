add_library(zernike zernike.cpp moments.cpp)
target_include_directories(zernike INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
