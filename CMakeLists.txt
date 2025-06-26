cmake_minimum_required (VERSION 3.5)
project (bindash)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fopenmp")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -fopenmp") # fast floating point arithmetics
add_executable (bindash ./src/bindash.cpp)
target_link_libraries(bindash -lm)
target_link_libraries(bindash -lz)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the "not version-controlled" changes with respect to the working branch
execute_process(
  COMMAND git diff HEAD --shortstat
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_DIFF_SHORTSTAT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
add_definitions("-DGIT_DIFF_SHORTSTAT=${GIT_DIFF_SHORTSTAT}")

