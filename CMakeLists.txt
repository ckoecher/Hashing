cmake_minimum_required(VERSION 3.3)
project(Hashing)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp ./inih/ini.c ./inih/cpp/INIReader.cpp PerfectHashFunction.cpp InputData.cpp)
add_executable(Hashing ${SOURCE_FILES})