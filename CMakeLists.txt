cmake_minimum_required(VERSION 3.13)
project(NLFEM)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

# set(CMAKE_BUILD_TYPE Debug)

# fetch Eigen - used for matrix/vector math and solving linear systems
message(CHECK_START "Fetching Eigen3")
list(APPEND CMAKE_MESSAGE_INDENT "  ")
include(FetchContent)
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
FetchContent_MakeAvailable(Eigen)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")

add_executable(NLFEM main.cpp element.cpp solver.cpp)
target_link_libraries(NLFEM Eigen3::Eigen)

add_executable(RadialReturn radial_return.cpp)
target_link_libraries(RadialReturn Eigen3::Eigen)