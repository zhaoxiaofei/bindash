cmake_minimum_required (VERSION 2.6)
project (binhash)
add_executable (binhash ./src/binhash.cpp)
set (CMAKE_CXX_FLAGS "-std=c++11")
target_link_libraries(binhash -lm)
target_link_libraries(binhash -lz)

