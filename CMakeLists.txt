cmake_minimum_required(VERSION 3.14)
project(StaticBTESolver VERSION 0.1 LANGUAGES CXX C Fortran)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake-modules")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)


set(PETSC_DIR "/home/yuehu/Downloads/petsc-3.14.2")
set(PETSC_ARCH "arch-linux-c-debug")


option(USE_GPU "compile as GPU version" off)
option(USE_MRM "compile as Memory_reduced_version" off)
option(USE_TIME "time each process" off)

if(DEFINED CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
else()
    SET (CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    include(CTest)
endif()

add_subdirectory(src)





