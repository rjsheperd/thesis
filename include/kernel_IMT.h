
#ifndef KERNEL_H_
#define KERNEL_H_
#include <math.h>
#include <stdio.h>
#include <kernel_common.h>

extern __device__ int end;
//const int INF = 32767;

__global__ void ItMinTime(int* , int* , int* ,
                          float* , int* , float* , bool* ,
                          int , int , int );

//__global__ void copyKernelIMT(int* , int* , bool* , int );


#endif // KERNEL_H_
