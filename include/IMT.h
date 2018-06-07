//
// Created by jsmith on 10/12/15.
//

#include "propagation.h"
#include "kernel_IMT.h"

#ifndef SIMULATOR_IMT_H
#define SIMULATOR_IMT_H

class IMT : Propagation
{
public:
   IMT(int,int, std::string, std::string);
   ~IMT();
   bool Init(std::string,std::string,std::string,std::string,
             std::string,std::string,std::string,std::string,
             std::string,std::string,std::string,int,int,int);
   bool CopyToDevice(int,int,int);
   bool RunKernel(int,int,int,bool,bool,int);
   bool CopyFromDevice();
   bool WriteToFile(std::string);
   bool WindXToFile(std::string, std::string* metaptr);
   bool WindYToFile(std::string, std::string* metaptr);
   bool UpdateCell(int,int,int);

protected:
   // Host Variables
   int* toa_map_;
   float* wind_x_map_;
   float* wind_y_map_;
   float* timesteppers_;
   bool* check_;
   //Device Variables
   int* g_toa_map_now_;
   int* g_toa_map_next_;
   int* g_toa_map_step_;
   float* g_wind_x_map_in_;
   float* g_wind_x_map_out_;
   float* g_wind_y_map_in_;
   float* g_wind_y_map_out_;
   int* g_timesteppers_;
   bool * g_check_;
};

#endif //SIMULATOR_IMT_H
