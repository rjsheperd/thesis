
#ifndef KERNEL_H_
#define KERNEL_H_
#include <math.h>
#include <stdio.h>
#include "kernel_common.h"

extern __device__ int end;


__global__ void MinTime(int*, float*, int*, float*,
                        int, int, int);
__global__ void TimeKernelMT(int* times);

#endif // KERNEL_H_


