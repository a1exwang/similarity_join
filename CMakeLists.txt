CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
PROJECT(similarity_join)

#SET (CMAKE_CXX_FLAGS "-std=c++11 -Wall -O3 -MMD -MP")
ADD_EXECUTABLE(main main.cpp SimJoiner.cpp)
