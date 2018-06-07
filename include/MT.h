//
// Created by jsmith on 10/10/15.
//

#ifndef SIMULATOR_MT_H
#define SIMULATOR_MT_H
#include "propagation.h"
#include "kernel_MT.h"

class MT : Propagation{
public:
   MT(int,int, std::string, std::string);
   ~MT();
   bool Init(std::string,std::string,std::string,std::string,
             std::string,std::string,std::string,std::string,
             std::string,std::string,std::string,int,int,int);
   bool CopyToDevice(int,int,int);
   bool RunKernel(int, int, int, bool, bool, int);
   bool CopyFromDevice();
   bool WriteToFile(std::string);
   bool WindXToFile(std::string, std::string* metaptr);
   bool WindYToFile(std::string, std::string* metaptr);
   bool UpdateCell(int,int,int);

protected:
   // Host Variabels
   int* toa_map_;
   float* wind_x_map_;
   float* wind_y_map_;
   int* timesteppers_;
   float* l_n_;
   // Device Variables
   int* g_toa_map_;
   float* g_wind_x_map_in_;
   float* g_wind_x_map_out_;
   float* g_wind_y_map_in_;
   float* g_wind_y_map_out_;
   int* g_timesteppers_;
};

#endif //SIMULATOR_MT_H
