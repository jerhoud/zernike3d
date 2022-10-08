add_library(tools iotools.cpp vec.cpp triangle.cpp)
target_include_directories(tools INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(tools Threads::Threads)
