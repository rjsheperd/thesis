//
// Created by jsmith on 10/2/15.
//

#include "propagation.h"
//#include <driver_types.h>
//#include <cuda_runtime_api.h>
/*
 *  Constructor
 */
Propagation::Propagation(int _x, int _y, std::string fuel_model_name, std::string fuel_moisture_name){
   // Point at appropriate rothdata files
   //printf("Propagation Constructor\n");
   simulation_ = new FireSim(_x, _y, fuel_model_name, fuel_moisture_name);
   // Initialize neighbor data
   sim_size_ = _x * _y;
   sim_rows_ = _y;
   sim_cols_ = _x;

   maxspreadrate_ = (float*) malloc(sim_size_*16*sizeof(float)); // don't need this anymore???
   curspreadrate_ = (float*) malloc(sim_size_*16*sizeof(float)); // don't need this anymore???
   intensity_modifier_ = (float*) malloc(sim_size_*sizeof(float)); // don't need this anymore???
   acceleration_constant_ = (float*) malloc(sim_size_*sizeof(float)); // don't need this anymore???
   // spotting additions
   emberMap_ = (bool*) malloc(sim_size_*16*sizeof(bool));
   emberTestMap_ = (float*) malloc(sim_size_*16*sizeof(float));
   ember_list_ = (Ember*) malloc(sim_size_*16*sizeof(Ember));
}

/*
 *  Destructor
 */
Propagation::~Propagation(){
   // Free CPU Memory
   free(maxspreadrate_);
   free(curspreadrate_);
   free(intensity_modifier_);
   free(acceleration_constant_);
   free(emberMap_);
   free(emberTestMap_);
   free(ember_list_);
   // Free GPU Memory
   cudaFree(g_maxspreadrate_);
   cudaFree(g_curspreadrate_);
   cudaFree(g_intensity_modifier_);
   cudaFree(g_acceleration_constant_);
   cudaFree(g_I_o_);
   cudaFree(g_RAC_);
   cudaFree(g_l_n_);
   cudaFree(g_emberMap_);
   cudaFree(g_emberTestMap_);
   cudaFree(g_canopyHeight_);
   cudaFree(g_ember_list_);
   cudaFree(g_x_wind);
   cudaFree(g_y_wind);
   cudaFree(g_DynLiveH);
   cudaFree(g_DynDead1);
   cudaFree(g_fuel_t_);
   cudaFree(g_roth_data_);
   cudaFree(g_dead_sav_burnable_b_);
   cudaFree(g_dead_1h_b_);
   cudaFree(g_dead_10h_b_);
   cudaFree(g_dead_100h_b_);
   cudaFree(g_live_h_b_);
   cudaFree(g_live_w_b_);
   cudaFree(g_fine_dead_extinctions_density_b_);
   cudaFree(g_areas_reaction_factors_b_);
   cudaFree(g_slope_wind_factors_b_);
   cudaFree(g_residence_flux_live_sav_b_);
   cudaFree(g_fuel_sav_accel_b_);
   cudaFree(g_slope_aspect_elevation_t_);
   cudaFree(g_dead_moistures_t_);
   cudaFree(g_live_moistures_t_);
   cudaFree(g_angles_);
   cudaFree(g_CellSize_);
   // Delete Simulation Memory
   delete simulation_;
}

/*
 *  Init()
 *  Functionality: Initializes all the values needed to be passed to the GPU.
 *                 This step is necessary for updating the simulation during the
 *                 runtime execution of the simulation.
 */
bool Propagation::Init(std::string fuel_file, std::string terrain_file,
                       std::string canopy_height_file, std::string crown_base_height_file,
                       std::string crown_bulk_density_file, std::string wind_x, std::string wind_y, std::string dynLiveH, std::string dynDead1,
                       std::string csr, std::string msr, int M_FLAG, int C_FLAG, int S_FLAG){
   /*
      cells formatting will be as follows:
      cell[x*8] is reference cell, following 8/16 vals are direction data
      N  NE   E  SE   S  SW   W  NW  NNW NNE NEE SEE SSE SSW SWW NWW
   */
   printf("Propgation Init\n");
   // Initialize Simulation
   simulation_->Init(fuel_file, terrain_file,
                     canopy_height_file, crown_base_height_file,
                     crown_bulk_density_file, wind_x, wind_y, dynLiveH, dynDead1, csr, msr,
                     M_FLAG, C_FLAG, S_FLAG);

   if(C_FLAG == 1 && S_FLAG == 1) // if crowning is off. spotting is off, if crowning is off spotting can still be off
   {                              // so check for both and if both are on then we do this
      int cell = 0;
      for(int row = 0; row < sim_rows_; row++)
      {
         for (int col = 0; col < sim_cols_; col++)
         {
            for(int dirs = 0; dirs < 16; dirs++, cell++)
            {
               // Assign default values
               emberMap_[cell] = false;
               emberTestMap_[cell] = 0.f;
               Ember launch(col, row, 0.f, 0.f, 0.0, 0.0);
               ember_list_[cell] = launch;
            }
         }
      }
   }
   return true;
}

/*
 *  CopyToDevice()
 *  Functionality: Copies memory from host to device for propagation
 */
bool Propagation::CopyToDevice(int M_FLAG, int C_FLAG, int S_FLAG) {
 
   printf("Propgation Copy To Device\n");

   // use this to catch erros and is good programming
   cudaError_t err;

   // ALLOCATING MEMORY ALLOCATING MEMORY ALLOCATING MEMORY ALLOCATING MEMORY ALLOCATING MEMORY ALLOCATING MEMORY
   if(M_FLAG == 1) // MOISTURE THINGS
   {
      err = cudaMalloc((void**) &g_DynLiveH, sim_size_ * sizeof(int));
      if (err != cudaSuccess) {
         std::cerr << "Error Allocating in propagation.cu: g_DynLiveH " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMalloc((void**) &g_DynDead1, sim_size_ * sizeof(int));
      if (err != cudaSuccess) {
         std::cerr << "Error Allocating in propagation.cu: g_DynDead1 " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
   }  
   err = cudaMalloc((void**) &g_maxspreadrate_, sim_size_ * 16 * sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_maxspreadrate_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_curspreadrate_, sim_size_ * 16 * sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_curspreadrate_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_intensity_modifier_, sim_size_ * sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_intensity_modifier_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_acceleration_constant_, sim_size_ * sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_acceleration_constant_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   if(C_FLAG == 1) // CROWNING THINGS
   {
      err = cudaMalloc((void**) &g_I_o_, sim_size_ * sizeof(float));
      if (err != cudaSuccess) {
         std::cerr << "Error Allocating in propagation.cu: g_I_o_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMalloc((void**) &g_RAC_, sim_size_ * sizeof(float));
      if (err != cudaSuccess) {
         std::cerr << "Error Allocating in propagation.cu: g_RAC_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMalloc((void**) &g_canopyHeight_, sim_size_ * sizeof(float));
      if (err != cudaSuccess) {
         std::cerr << "Error Allocating in propagation.cu: g_canopyHeight_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
   }
   err = cudaMalloc((void**) &g_l_n_, 16 * sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_l_n_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   if(C_FLAG == 1 && S_FLAG == 1) // SPOTTING THINGS
   {
      err = cudaMalloc((void**) &g_emberMap_, sim_size_ * 16 * sizeof(bool));
      if (err != cudaSuccess) {
         std::cerr << "Error Allocating in propagation.cu: g_emberMap_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMalloc((void**) &g_emberTestMap_, sim_size_ * 16 * sizeof(float));
      if (err != cudaSuccess) {
         std::cerr << "Error Allocating in propagation.cu: g_emberTestMap_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMalloc((void**) &g_ember_list_, sim_size_ * 16 * sizeof(Ember));
      if (err != cudaSuccess) {
         std::cerr << "Error Allocating in propagation.cu: g_ember_list_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
   }
   err = cudaMalloc((void**) &g_x_wind, sim_size_ * sizeof(int));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_x_wind " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_y_wind, sim_size_ * sizeof(int));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_y_wind " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   // added for moisture manipulation
   err = cudaMalloc((void**) &g_fuel_t_, sim_size_ * sizeof(int));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_fuel_t_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_roth_data_, sim_size_ * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_roth_data_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_dead_sav_burnable_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_dead_sav_burnable_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_dead_1h_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_dead_1h_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_dead_10h_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_dead_10h_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_dead_100h_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_dead_100h_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_live_h_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_live_h_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_live_w_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_live_w_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_fine_dead_extinctions_density_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_fine_dead_extinctions_density_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_areas_reaction_factors_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_areas_reaction_factors_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_slope_wind_factors_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_slope_wind_factors_b_" << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_residence_flux_live_sav_b_, simulation_->numModels * sizeof(float4));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_residence_flux_live_sav_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_fuel_sav_accel_b_, simulation_->numModels * sizeof(float2));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_fuel_sav_accel_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_slope_aspect_elevation_t_, sim_size_ * sizeof(float3));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_slope_aspect_elevation_t_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_dead_moistures_t_, simulation_->numMoistModels * sizeof(float3));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_dead_moistures_t_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_live_moistures_t_, simulation_->numMoistModels * sizeof(float2));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating in propagation.cu: g_live_moistures_t_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_angles_, sim_size_ * sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating Memory in propagation.cu: g_angles_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_CellSize_, sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating Memory in propagation.cu: g_CellSize_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_maxSpread_, sizeof(float));
   // end additions for moisture manipulation
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating Memory in propagation.cu: g_maxSpread_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMalloc((void**) &g_timeStep_, sizeof(float));
   // end additions for moisture manipulation
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating Memory in propagation.cu: g_timeStep_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }

   // TRANSFERRING MEMORY TRANSFERRING MEMORY TRANSFERRING MEMORY TRANSFERRING MEMORY TRANSFERRING MEMORY

   err = cudaMemcpy(g_maxspreadrate_, simulation_->msr_, sim_size_* 16 * sizeof(float), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_maxspreadrate_" << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_curspreadrate_, simulation_->csr_, sim_size_* 16 *sizeof(float), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_curspreadrate_" << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_intensity_modifier_, intensity_modifier_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_intensity_modifier_" << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_acceleration_constant_, acceleration_constant_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_acceleration_constant_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_l_n_, simulation_->l_n_, 16*sizeof(float), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_l_n_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   if(C_FLAG == 1)
   {
      err = cudaMemcpy(g_I_o_, simulation_->I_o_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
      if (err != cudaSuccess) {
         std::cerr << "Error Copying in propagation.cu: g_I_o_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMemcpy(g_RAC_, simulation_->RAC_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
      if (err != cudaSuccess) {
         std::cerr << "Error Copying in propagation.cu: g_RAC_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMemcpy(g_canopyHeight_, simulation_->canopy_height_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
      if (err != cudaSuccess) {
         std::cerr << "Error Copying in propagation.cu: g_canopyHeight_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
   }
   if(C_FLAG == 1 && S_FLAG == 1)
   {
      err = cudaMemcpy(g_emberMap_, emberMap_, sim_size_*16*sizeof(bool), cudaMemcpyHostToDevice);
      if (err != cudaSuccess) {
         std::cerr << "Error Copying in propagation.cu: g_emberMap_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMemcpy(g_emberTestMap_, emberTestMap_, sim_size_*16*sizeof(float), cudaMemcpyHostToDevice);
      if (err != cudaSuccess) {
         std::cerr << "Error Copying in propagation.cu: g_emberTestMap_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMemcpy(g_ember_list_, ember_list_, sim_size_*16*sizeof(Ember), cudaMemcpyHostToDevice);
      if (err != cudaSuccess) {
         std::cerr << "Error Copying in propagation.cu: g_ember_list_ " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
   }
   err = cudaMemcpy(g_x_wind, simulation_->xwind_, sim_size_*sizeof(int), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_x_wind " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_y_wind, simulation_->ywind_, sim_size_*sizeof(int), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_y_wind " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   // added for moisture manipulation
   err = cudaMemcpy(g_fuel_t_, simulation_->fuel_t_, sim_size_*sizeof(int), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_fuel_t_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_roth_data_, simulation_->roth_data_, sim_size_*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_roth_data_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_dead_sav_burnable_b_, simulation_->dead_sav_burnable_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_dead_sav_burnable_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_dead_1h_b_, simulation_->dead_1h_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_dead_1h_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_dead_10h_b_, simulation_->dead_10h_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_dead_10h_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_dead_100h_b_, simulation_->dead_100h_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_dead_100h_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_live_h_b_, simulation_->live_h_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_live_h_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_live_w_b_, simulation_->live_w_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_live_w_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_fine_dead_extinctions_density_b_, simulation_->fine_dead_extinctions_density_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_fine_dead_extinctions_density_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_areas_reaction_factors_b_, simulation_->areas_reaction_factors_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_areas_reaction_factors_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_slope_wind_factors_b_, simulation_->slope_wind_factors_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_slope_wind_factors_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_residence_flux_live_sav_b_, simulation_->residence_flux_live_sav_b_, simulation_->numModels*sizeof(float4), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_residence_flux_live_sav_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_fuel_sav_accel_b_, simulation_->fuel_sav_accel_b_, simulation_->numModels*sizeof(float2), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_fuel_sav_accel_b_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_slope_aspect_elevation_t_, simulation_->slope_aspect_elevation_t_, sim_size_*sizeof(float3), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_slope_aspect_elevation_t_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_dead_moistures_t_, simulation_->dead_moistures_t_, simulation_->numMoistModels*sizeof(float3), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_dead_moistures_t_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_live_moistures_t_, simulation_->live_moistures_t_, simulation_->numMoistModels*sizeof(float2), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_live_moistures_t_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_angles_, simulation_->angles_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   // end additions for moisture manipulation
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_angles_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   err = cudaMemcpy(g_CellSize_, &simulation_->cell_size_, 1*sizeof(float), cudaMemcpyHostToDevice);
   // end additions for moisture manipulation
   if (err != cudaSuccess) {
      std::cerr << "Error Copying in propagation.cu: g_CellSize_ " << cudaGetErrorString(err) << std::endl;
      exit(1);
   }
   if(M_FLAG == 1)
   {
      err = cudaMemcpy(g_DynLiveH, simulation_->DynLiveH, sim_size_*sizeof(int), cudaMemcpyHostToDevice);
      // end additions for moisture manipulation
      if (err != cudaSuccess) {
         std::cerr << "Error Copying in propagation.cu: g_DynLiveH " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
      err = cudaMemcpy(g_DynDead1, simulation_->DynDead1, sim_size_*sizeof(int), cudaMemcpyHostToDevice);
      // end additions for moisture manipulation
      if (err != cudaSuccess) {
         std::cerr << "Error Copying in propagation.cu: g_DynDead1 " << cudaGetErrorString(err) << std::endl;
         exit(1);
      }
   }
   return true;
}