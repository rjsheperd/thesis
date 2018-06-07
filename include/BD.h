//
// Created by jsmith on 10/2/15.
//

#ifndef SIMULATOR_BD_H
#define SIMULATOR_BD_H

#include "propagation.h"
#include "kernel_BD.h"

class BD : Propagation
{
public:
   BD(int,int, std::string, std::string);
   ~BD();
   bool Init(std::string,std::string,std::string,std::string,
             std::string,std::string,std::string,std::string,
             std::string,std::string,std::string,int,int,int);
   bool CopyToDevice(int,int,int);
   bool RunKernel(int,int,int,bool,bool,int);
   bool CopyFromDevice();
   bool UpdateCell(int,int,int);
   bool WriteToFile(std::string);
   bool WritePauseToFile(std::string, std::string* metaptr, std::string, std::string);

protected:
   // Host Variables
   float*   toa_map_;
   int*     wind_x_map_;
   int*     wind_y_map_;
   float*   loc_burndist_;
   float    current_time_;
   float    time_step_;
   float*   emberValues_;
   Ember*   embers_;
   int*     deadHr1_;
   int*     liveH_;
   float*   maxSpreadRates_;
   float*   currentSpreadRates_;


   // Device Variables
   float* g_toa_map_in_;
   float* g_toa_map_out_;
   float* g_loc_burndist_;
   Ember* g_embers_;
};
#endif //SIMULATOR_BD_H
