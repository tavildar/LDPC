cmake_minimum_required(VERSION 3.5)
project(LdpcC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp LdpcCode.cpp LdpcCode.h WiFiLDPC.h Constellation.cpp Constellation.h)
add_executable(LdpcC ${SOURCE_FILES})