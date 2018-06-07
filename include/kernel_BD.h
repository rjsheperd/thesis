#ifndef KERNEL_BD_H_
#define KERNEL_BD_H_
#include <math.h>
#include <stdio.h>
#include "kernel_common.h"
//#include "device_launch_parameters.h"
//__device__ int end;
//__device__ int syncCounter;

__global__ void BurnDist(float*, float*,float*,
                         float*,float*, int,
                         int, int, float, float);

__global__ void copyKernelBD(float* input, float* output, int size);

__global__ void updateTOA(Ember*, float*, int);

#endif // KERNEL_BD_H_
