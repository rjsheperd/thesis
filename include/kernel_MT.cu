#include "kernel_MT.h"

/////////////////////////////////////////////////////////////////////////////
//                              Minimal Time
/////////////////////////////////////////////////////////////////////////////
__global__ void MinTime(int* ignTime, float* rothData, int* times,
                        float* L_n, int size, int rowSize,
                        int colSize){
   /* neighbor's address*/     /* N  NE   E  SE   S  SW   W  NW  NNW NNE NEE SEE SSE SSW SWW NWW*/
   int nCol[16] =        {  0,  1,  1,  1,  0, -1, -1, -1, -1, 1, 2, 2, 1, -1, -2, -2};
   int nRow[16] =        {  1,  1,  0, -1, -1, -1,  0,  1, 2, 2, 1, -1, -2, -2, -1, 1};

   // Calculate ThreadID
   int cell = blockIdx.x * blockDim.x + threadIdx.x;
   int ncell, nrow, ncol, row, col;
   float ROS;
   int ignCell, ignCellN, timeNow, timeNext;
   int rothCell = 0;

   timeNow = times[0]; // timeNow = time_next_
   timeNext = INF;

   while(cell < size){
      row = cell / colSize;
      col = cell % colSize;
      ignCell = ignTime[cell];

      // Do atomic update of TimeNext Var (atomicMin)
      if(timeNext > ignCell && ignCell > timeNow){
         atomicMin(&times[1], ignCell);
         timeNext = ignCell;
      }
      else if(ignCell == timeNow){ // I am on fire now, and will propagate
         // Check through neighbors
         for(int n = 0; n < 16; n++){
            // // Propagate from burning cells
            rothCell = cell * 16;
            nrow = row + nRow[n];
            ncol = col + nCol[n];
            // printf("nrow: %d ncol: %d\n",nrow,ncol);
            if ( nrow<0 || nrow>= rowSize || ncol<0 || ncol>=  colSize ){
               continue;
            }
            ncell = ncol + nrow*colSize;
            ignCellN = ignTime[ncell];

            // If neighbor is unburned in this timestep
            if(ignCellN > timeNow){
               // compute ignition time
               ROS = rothData[rothCell + n];
               float ignTimeNew = (timeNow + (L_n[n] / ROS));
               if(ignTimeNew <=0)
                  break;
               // Update Output TOA Map
               int old = atomicMin(&ignTime[ncell], (int)ignTimeNew);
               if(old != ignTimeNew)
                  rothData[ncell] = ROS;
               // Local time_next_ update
               if((int)ignTimeNew < timeNext && ignTimeNew > 0){ // #thisisnotahacklol
                  timeNext = (int)ignTimeNew;
               }
            }
         }
         // Perform global time_next_ update
         atomicMin(&times[1], timeNext);
      }
      // Do striding
      cell += blockDim.x * gridDim.x;
   }

}

/////////////////////////////////////////////////////////////////////////////
//                             Time Update (MT)
/////////////////////////////////////////////////////////////////////////////
__global__ void TimeKernelMT(int* times){
   times[0] = times[1];
//   if(times[0] == 0)
   times[1] = INF;
}