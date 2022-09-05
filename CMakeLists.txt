cmake_minimum_required(VERSION 3.22)
project(final_project C)

set(CMAKE_C_STANDARD 90)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES})

add_executable(final_project spkmeans.c spkmeans.h kmeans.c matrixOperations.c
        matrixOperations.h vectorOperations.c vectorOperations.h spkmeansmodule.c)

