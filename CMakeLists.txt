cmake_minimum_required(VERSION 2.8)
project(simulator)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES -O2 -std=c++11")

include_directories(include)

set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FireSim.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FireSim.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/kernel_BD.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/include/kernel_BD.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/kernel_IMT.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/include/kernel_MT.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/include/kernel_common.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/include/vector.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MT.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/IMT.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/IMT.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/include/BD.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/include/BD.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MT.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MT.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/SEQ.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/SEQ.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/propagation.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/include/propagation.h
)

find_package(CUDA QUIET REQUIRED)

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};)
include_directories(${CUDA_INCLUDE_DIRS})
set(CUDALIBS ${CUDALIBS} ${CUDA_LIBRARIES})

find_package(GDAL REQUIRED)
include_directories(${GDAL_INCLUDE_DIR})
set(GDALLIBS ${GDALLIBS} ${GDAL_LIBRARY})

cuda_add_executable(simulator ${SOURCE_FILES})

target_link_libraries(simulator ${CUDALIBS} ${GDALLIBS})