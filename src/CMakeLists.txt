cmake_minimum_required(VERSION 3.17)
project(GraphCut)

set(CMAKE_CXX_STANDARD 14)

include_directories(include/)

add_executable(GraphCut main.cpp include/GraphCutTexture.hpp include/MInCut.hpp include/MaxFlow.hpp include/Graph.hpp)

find_package(OpenCV)

target_link_libraries(GraphCut ${OpenCV_LIBS})