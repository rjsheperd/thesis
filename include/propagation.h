//
// Created by jsmith on 10/2/15.
//

#ifndef SIMULATOR_PROPAGATION_H
#define SIMULATOR_PROPAGATION_H
#include <iostream>
#include <fstream>
#include "FireSim.h"

// added for spotting
#include "Ember.h"

class Propagation{
public:
    Propagation(int,int, std::string, std::string);
    ~Propagation();
    virtual bool Init(std::string,std::string,std::string,std::string,
                      std::string,std::string,std::string,std::string,
                      std::string,std::string,std::string,int,int,int);
    virtual bool CopyToDevice(int,int,int);

protected:
    // Host Variables
    FireSim*    simulation_;
    float*      maxspreadrate_;
    float*      curspreadrate_;
    float*      intensity_modifier_;
    float*      acceleration_constant_;
    float*      l_n_;
    int         sim_size_;
    int         sim_rows_;
    int         sim_cols_;
    Ember*      ember_list_;
    float*      emberTestMap_;
    bool*       emberMap_;

    // Device Variables
    float*  g_maxspreadrate_;
    float*  g_curspreadrate_;
    float*  g_intensity_modifier_;
    float*  g_acceleration_constant_;
    float*  g_I_o_;
    float*  g_RAC_;
    float*  g_l_n_;
    float*  g_canopyHeight_;
    float*  g_emberTestMap_;
    bool*   g_emberMap_;
    Ember*  g_ember_list_;
    int*    g_fuel_t_;
    float4* g_roth_data_;
    float4* g_dead_sav_burnable_b_;
    float4* g_dead_1h_b_;
    float4* g_dead_10h_b_;
    float4* g_dead_100h_b_;
    float4* g_live_h_b_;
    float4* g_live_w_b_;
    float4* g_fine_dead_extinctions_density_b_;
    float4* g_areas_reaction_factors_b_;
    float4* g_slope_wind_factors_b_;
    float4* g_residence_flux_live_sav_b_;
    float2* g_fuel_sav_accel_b_;
    float3* g_slope_aspect_elevation_t_;
    float3* g_dead_moistures_t_;
    float2* g_live_moistures_t_;
    float*  g_angles_;
    float*  g_CellSize_;
    float*  g_maxSpread_;
    float*  g_timeStep_;
    int*    g_x_wind;
    int*    g_y_wind;
    int*    g_DynLiveH;
    int*    g_DynDead1;
};
#endif //SIMULATOR_PROPAGATION_H