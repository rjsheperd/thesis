#ifndef KERNEL_COMMON_H
#define KERNEL_COMMON_H
#include <math.h>
#include <stdio.h>
#include "Ember.h"
#include <curand.h>
#include <curand_kernel.h>

const int INF = 32767;

const int MAX = 100;

extern __device__ int spock = 0; // for testing purposes

extern __device__ bool dynamicMositure;

extern __device__ bool precipitation;

__global__ void InitialSetup(int*, float4*, float4*, float4*, float4*, float4*, 
							 float4*, float4*, float4*, float4*, float4*, float4*,
							 float2*, float3*, float3*, float2*, int*, int*, int, int*, int*);

__global__ void UpdateSpread(int*, float4*, float2*, float*, float*, float*, float*, int);

__global__ void Accelerate(float*, float*, int, float);

__global__ void TestCrownRate_NoSpot(float*, float*, float*, int, float*, float*);

__global__ void TestCrownRate_WithSpot(float*, float*, float*, int, float*,
                               float*, bool*, float*, float*, Ember*, int*, int*, int*);

__global__ void TestSpotting(bool*, float*, Ember*, float, float, int, int*, int*, float3*);

__global__ void findMax(float*, float*, float*, float*, int);

__device__ float findMaxHelper(float*, float);

__device__ float Clamp(float, float, float);

// __global__ int rand(curanStat)
// http://cs.umw.edu/~finlayson/class/fall16/cpsc425/notes/cuda-random.html
// https://stackoverflow.com/questions/18501081/generating-random-number-within-cuda-kernel-in-a-varying-range
// http://aresio.blogspot.com/2011/05/cuda-random-numbers-inside-kernels.html
// http://www.fbfrg.org/fuel-moisture/1-hour-moisture-content
#endif