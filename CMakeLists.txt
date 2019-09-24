cmake_minimum_required(VERSION 3.1)
project(exact-ccd)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_library(exact-ccd
expansion.cpp  interval.cpp  rootparitycollisiontest.cpp bsc.cpp
    )
target_include_directories(exact-ccd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/)