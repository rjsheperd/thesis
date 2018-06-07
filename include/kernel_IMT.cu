#include "kernel_IMT.h"
/////////////////////////////////////////////////////////////////////////////
//                          Iterative Minimal Time
/////////////////////////////////////////////////////////////////////////////
__global__ void ItMinTime(int* ignTimeIn, int* ignTimeOut, int* ignTimeStep,
                          float* rothData, int* times, float* L_n, bool* check,
                          int size, int rowSize, int colSize){
   /* neighbor's address*/     /* N  NE   E  SE   S  SW   W  NW  NNW NNE NEE SEE SSE SSW SWW NWW*/
   int nCol[16] =        {  0,  1,  1,  1,  0, -1, -1, -1, -1, 1, 2, 2, 1, -1, -2, -2};
   int nRow[16] =        {  1,  1,  0, -1, -1, -1,  0,  1, 2, 2, 1, -1, -2, -2, -1, 1};
   float ignCell = 0;
   float ignCellNew = 0;
   float ignTimeMin = INF;

   int cell = blockIdx.x * blockDim.x + threadIdx.x;
   int ncell, nrow, ncol, row, col;
   float ignTimeNew, ROS, ROS_Update;

   while(cell < size){
      row = cell / colSize;
      col = cell % colSize;

      // Do nothing if converged
      if(check[cell] == true){
         cell += blockDim.x * gridDim.x;
         continue;
      }

      // Check for simulation completion
      ignCell = ignTimeIn[cell];
      ignCellNew = ignTimeOut[cell];
      // Convergence Test
      if(fabs(ignCell - ignCellNew) < 2 && ignCell != INF
         && ignCellNew != INF && check[cell] != true){
         check[cell] = true;
         cell += blockDim.x * gridDim.x;
         continue;
      }
      if(ignCell > 0){
         ignTimeMin = INF;
         // Loop through neighbors
         for(int n = 0; n < 16; n++){
            // find neighbor cell index
            nrow = row + nRow[n];
            ncol = col + nCol[n];
            if ( nrow<0 || nrow>= rowSize || ncol<0 || ncol>=  colSize ){
               continue;
            }
            ncell = ncol + nrow*colSize;

            ROS = rothData[ncell * 16 + n];
            ignTimeNew = ignTimeIn[ncell] + L_n[n] / ROS;// * 100;

            ignTimeMin = ignTimeNew*(ignTimeNew < ignTimeMin) + ignTimeMin*(ignTimeNew >= ignTimeMin);
            ROS_Update = ROS*(ignTimeNew < ignTimeMin) + rothData[cell]*(ignTimeNew >= ignTimeMin);
         }
         if(ignTimeMin >0){
            ignTimeOut[cell] = (int)ignTimeMin;
            rothData[cell] = ROS_Update;
         }
      }
      cell += blockDim.x * gridDim.x;
   }
   if(blockIdx.x * blockDim.x + threadIdx.x == 0)
      end = 0;
}