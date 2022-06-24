cmake_minimum_required(VERSION 2.4)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

SET(CMAKE_CXX_FLAGS_DEBUG " -O0 -g" )
SET(CMAKE_C_FLAGS_DEBUG " -O0 -g" )

SET(CMAKE_CXX_FLAGS_RELEASE " -O3 -msse2 -mssse3 -fexpensive-optimizations -Wuninitialized" )
SET(CMAKE_C_FLAGS_RELEASE " -O3 -msse2 -mssse3 -fexpensive-optimizations -Wuninitialized" )

# SET(CMAKE_CXX_FLAGS_RELEASE " -O0 -g" )
# SET(CMAKE_C_FLAGS_RELEASE " -O0 -g" )
# set(CMAKE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-static")
# set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")

project(polysim)

add_library(r250 r250.cpp)
add_library(numerics numerics.cpp)
add_executable(polysim main.cpp)

target_link_libraries(polysim r250 numerics)
