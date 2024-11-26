cmake_minimum_required(VERSION 3.10)

project(DG VERSION 1.0 LANGUAGES CXX)

include(FetchContent)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -funroll-loops")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(src)
