# Written by J. Houdayer

add_library(tools iotools.cpp)
target_include_directories(tools INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
if (USE_THREADS)
    target_link_libraries(tools INTERFACE Threads::Threads)
endif()
