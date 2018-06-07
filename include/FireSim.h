//
// Created by jsmith on 10/2/15.
//

//
// updated by andyg in SPRING 2018.
//

#ifndef SIMULATOR_FIRESIM_H
#define SIMULATOR_FIRESIM_H
#include <iostream>
#include <fstream>
#include <vector>
#include "FuelModel.h"
#include "FuelMoisture.h"
#include <string>
#include <gdal.h>
#include <gdal_priv.h>
#include <cpl_conv.h>
#include "vector.h"
#include <map>
#include "vector_types.h"
#include <stdlib.h>

// experimental #include "vector_types.h"
// https://devtalk.nvidia.com/default/topic/416771/copying-a-cuda-float4-array-back-to-host-using-memcpy/
// https://devtalk.nvidia.com/default/topic/397235/float4-linear-array-or-cudaarray-/

class FireSim {
public:
   FireSim(int _x,int _y, std::string, std::string);
   ~FireSim();
   void     Init(std::string, std::string, std::string, std::string,
                 std::string, std::string, std::string, std::string,
                 std::string, std::string, std::string, int, int, int);
   void     UpdateSpreadData();
   float    Clamp(float, float, float);
   bool     BurnDistance(float &, float, float);

   // Simulation Data
   // roth_data_ --> x - maxSpreadRate , y - spreadDirection , z - ellipseEccentricity , w - IntensityModifier

   // CUDA does not deal with vector stl
   float4*  roth_data_;
   int*     DynLiveH;
   int*     DynDead1;
   int*     xwind_;
   int*     ywind_;
   float*   I_o_;
   float*   RAC_;
   float*   canopy_height_;
   float*   csr_;
   float*   msr_;
   float    acceleration_constant_;
   float    foliar_moisture;

   // Simulation members from test Sim
   float    time_now_;
   float    time_next_;
   float*   ign_time_;
   float*   ign_time_new_;
   float**  burn_dist_;
   float*   l_n_;

   // Rothermel Data Members
   int*        fuel_t_;
   std::string root_path_;
   float4*     dead_sav_burnable_b_;
   float4*     dead_1h_b_;
   float4*     dead_10h_b_;
   float4*     dead_100h_b_;
   float4*     live_h_b_;
   float4*     live_w_b_;
   float4*     fine_dead_extinctions_density_b_;
   float4*     areas_reaction_factors_b_;
   float4*     slope_wind_factors_b_;
   float4*     residence_flux_live_sav_b_;
   float2*     fuel_sav_accel_b_;
   float3*     slope_aspect_elevation_t_;
   float3*     dead_moistures_t_;
   float2*     live_moistures_t_;
   float*      angles_;
   int         numModels;
   int         numMoistModels;
   int         sim_dim_x_;
   int         sim_dim_y_;
   float       cell_size_;
   float       time_step_;

   std::vector<sim::FuelModel> _models;
   std::vector<sim::FuelMoisture> _moistures;
};

#endif //SIMULATOR_FIRESIM_H