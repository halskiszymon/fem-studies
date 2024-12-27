cmake_minimum_required(VERSION 3.29)
project(mes)

set(CMAKE_CXX_STANDARD 20)

add_executable(mes main.cpp
        structures.h
        metodaGaussa.h
        jakobian.h
        pomocnicze.h)
