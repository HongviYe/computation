cmake_minimum_required(VERSION 3.1) # 更新CMake最低版本要求
project(Remesh LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置C++编译器优化级别为O2
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")

# 设置C编译器优化级别为O2
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast")

add_executable(project main.cpp)