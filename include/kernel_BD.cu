#include <curand_mtgp32_kernel.h>
#include "kernel_BD.h"

/////////////////////////////////////////////////////////////////////////////
//                            Burning Distances
/////////////////////////////////////////////////////////////////////////////
__global__ void BurnDist(float* ignTimeIn, float* ignTimeOut, float* burnDist,
                         float* currentspreadrate, float* L_n, int size,
                         int rowSize, int colSize, float timeStep, float t)
{

                     /* g_toa_map_in_     =  ignTimeIn 
                        g_toa_map_out_    =  ignTimeOut
                        g_loc_burndist_   =  burnDist
                        g_maxspreadrate_  =  maxspreadrate
                        g_curspreadrate_  =  currentspreadrate
                        g_l_n_            =  L_n
                        sim_size_         =  size
                        sim_rows_         =  rowSize
                        sim_cols_         =  colSize
                        timestep_         =  timeStep
                        current_time_     =  t
                     */

   /* neighbor's address*/     /* N  NE   E  SE   S  SW   W  NW  NNW NNE NEE SEE SSE SSW SWW NWW*/
   int nCol[8] =        {  0,  1,  1,  1,  0, -1, -1, -1};
   int nRow[8] =        {  1,  1,  0, -1, -1, -1,  0,  1};
  
   float ignTime, ignTimeN;
   float dist;

   int cell = blockIdx.x * blockDim.x + threadIdx.x;
   
   int ncell, nrow, ncol, row, col, distCell, rothCell;
   float ROS;

   while(cell < size){
      row = cell / colSize;
      col = cell % colSize;

      ignTime = ignTimeIn[cell];

      if(ignTime == INF){
         cell += blockDim.x * gridDim.x;
         // printf("g_fuel_t_: %i\n", g_fuel_t_[cell]);
         continue;
      }

      // check neighbors for ignition
      for(int n = 0; n < 8; n++){
         // int current_val = atomicAdd(&spock, 1);
         distCell = cell * 8;
         rothCell = cell * 16;

         nrow = row + nRow[n];
         ncol = col + nCol[n];
         if ( nrow<0 || nrow>=rowSize || ncol<0 || ncol>=colSize ) {
            continue;
         }
         ncell = ncol + nrow*colSize;

         // check for already lit or its toa is greater than the current timestep
         ignTimeN = ignTimeIn[ncell];
         // printf("ignTimeN: %f\n", ignTimeN);
         if(ignTimeN < INF /*&& ignTimeN <= timeStep*/){
            continue;
         }

         // if(isnan(currentspreadrate[rothCell + n]) || isinf(currentspreadrate[rothCell + n]) || currentspreadrate[rothCell + n] < 0)
         //    printf("currentspreadrate[rothCell + n]: %f\n", currentspreadrate[rothCell + n]);

         // if(isnan(currentspreadrate[rothCell + n]))
         //    printf("currentspreadrate[rothCell + n]: %f\n", currentspreadrate[rothCell + n]);

         // Calc roth values
         ROS = currentspreadrate[rothCell + n];

         // if(isnan(ROS) || isinf(ROS)){
         //    printf("KBD1 before atomicExch--currentspreadrate[cell]: %f\n;", ROS);
         //    asm("trap;");
         // }

         // if(isnan(currentspreadrate[rothCell + n]) || isinf(currentspreadrate[rothCell + n])){
         //    printf("KBD2 before atomicExch--currentspreadrate[rothCell + n]: %f\n;", currentspreadrate[rothCell + n]);
         //    asm("trap;");
         // }

         // Burn distance
         dist = burnDist[distCell+n];
         dist = dist - ROS*timeStep;
         burnDist[distCell+n] = dist;
         
         // Propagate fire
         if(dist <= 0){
            dist *= -1;
            float step_time = dist / ROS;
            step_time += t;
            float old = atomicExch(&ignTimeOut[ncell], step_time);
            if(old < step_time){
               atomicExch(&ignTimeOut[ncell], old);
               currentspreadrate[ncell] = ROS;
               // if(isnan(currentspreadrate[ncell]) || isinf(currentspreadrate[ncell])){
               //    printf("KBD3 After atomicExch--currentspreadrate[cell]: %f\n;", currentspreadrate[ncell]);
               //    asm("trap;");
               // }
            }
         }
      }
      cell += (blockDim.x * gridDim.x);
   }
   // printf("spock: %i\n", spock);
}

/////////////////////////////////////////////////////////////////////////////
//                 Update TOA with Spotting Ember TOA values
/////////////////////////////////////////////////////////////////////////////
__global__ void updateTOA(Ember* ember_list, float* g_toa_map_out_, int size)
{
   int dirs = 16;
   int cell = (blockIdx.x * blockDim.x + threadIdx.x) * dirs;
   int cell_im = (blockIdx.x * blockDim.x + threadIdx.x);

   while(cell < size*dirs)
   {
      for(int i = 0; i < dirs; i++)
      {
         if(ember_list[cell+i].ember_toa < g_toa_map_out_[cell_im] && g_toa_map_out_[cell_im] != 32767)
         {
            g_toa_map_out_[cell_im] = ember_list[cell+i].ember_toa;
         }
      }
      cell_im += (blockDim.x * gridDim.x);
      cell += (blockDim.x * gridDim.x)*dirs;
   }
}

/////////////////////////////////////////////////////////////////////////////
//                             Copy Kernel (BD)
/////////////////////////////////////////////////////////////////////////////
__global__ void copyKernelBD(float* input, float* output, int size)
{
   int cell_im = (blockIdx.x * blockDim.x + threadIdx.x);

   while(cell_im < size)
   {
      input[cell_im] = output[cell_im];
      cell_im += (blockDim.x * gridDim.x);
   }
}