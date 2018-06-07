//
// Created by jsmith on 10/2/15.
//
#include "BD.h"

BD::BD(int _x, int _y, std::string fuel_model_name, std::string fuel_moisture_name)
        : Propagation(_x, _y, fuel_model_name, fuel_moisture_name){
   toa_map_ = (float*) malloc(sim_size_ * sizeof(float));
   wind_x_map_ = (int*) malloc(sim_size_ * sizeof(int));
   wind_y_map_ = (int*) malloc(sim_size_ * sizeof(int));
   loc_burndist_ = (float*) malloc(8*sim_size_*sizeof(float));
   emberValues_ = (float*) malloc(sim_size_*16*sizeof(float));
   embers_ = (Ember*) malloc(sim_size_*16*sizeof(Ember));
   maxSpreadRates_ = (float*) malloc(sim_size_*16*sizeof(float));
   currentSpreadRates_ = (float*) malloc(sim_size_*16*sizeof(float));
}

BD::~BD(){
   // Free host memory
   free(toa_map_);
   free(wind_x_map_);
   free(wind_y_map_);
   free(emberValues_);
   free(embers_);
   // Free device memory
   cudaFree(g_toa_map_in_);
   cudaFree(g_toa_map_out_);
}

bool BD::Init(std::string fuel_file, std::string terrain_file,
              std::string canopy_height_file, std::string crown_base_height_file,
              std::string crown_bulk_density_file, std::string wind_x, std::string wind_y, std::string dynLiveH,
              std::string dynDead1, std::string csr, std::string msr, int M_FLAG, int C_FLAG, int S_FLAG){
   Propagation::Init(fuel_file, terrain_file,
                     canopy_height_file, crown_base_height_file,
                     crown_bulk_density_file, wind_x, wind_y, dynLiveH, dynDead1, csr, msr,
                     M_FLAG, C_FLAG, S_FLAG);

   for(unsigned int i = 0; i < sim_size_; i++){
      toa_map_[i] = simulation_->ign_time_[i];
      // std::cout<<toa_map_[i]<<std::endl;
   }
   // Populate Burn Distances
   for(unsigned int i = 0; i < sim_size_; i++){
      for(int j = 0; j < 8; j++){
         loc_burndist_[i*8+j] = simulation_->l_n_[j];
      }
   }
   current_time_ = 0.0f;
   return true;
}

bool BD::CopyToDevice(int M_FLAG, int C_FLAG, int S_FLAG){
   Propagation::CopyToDevice(M_FLAG, C_FLAG, S_FLAG);
   // Create memory on device
   cudaError_t err = cudaMalloc( (void**) &g_toa_map_in_, sim_size_*sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating Memory in BD.cu g_toa_map_in_ : " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }
   err = cudaMalloc( (void**) &g_toa_map_out_, sim_size_*sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating Memory in BD.cu g_toa_map_out_ : " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }
   err = cudaMalloc( (void**) &g_loc_burndist_, sim_size_*8*sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating Memory in BD.cu g_loc_burndist_ : " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }
   if(C_FLAG == 1 && S_FLAG == 1)
   {
      err = cudaMalloc( (void**) &g_embers_, sim_size_*16*sizeof(Ember));
      if (err != cudaSuccess) {
         std::cerr << "Error Allocating Memory in BD.cu g_embers_ : " << cudaGetErrorString(err) << std::endl;
         exit(1);
         return false;
      }
   }
   err = cudaMemcpy(g_toa_map_in_, toa_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Transferring Memory in BD.cu g_toa_map_in_: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }
   err = cudaMemcpy(g_toa_map_out_, toa_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Transferring Memory in BD.cu g_toa_map_out_: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }
   err = cudaMemcpy(g_loc_burndist_, loc_burndist_, sim_size_*8*sizeof(float), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Transferring Memory in BD.cu g_loc_burndist_: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }
   if(C_FLAG == 1 && S_FLAG == 1)
   {
      err = cudaMemcpy(g_embers_, embers_, sim_size_*16*sizeof(float), cudaMemcpyHostToDevice);
      if (err != cudaSuccess) {
         std::cerr << "Error Transferring Memory in BD.cu g_embers_" << cudaGetErrorString(err) << std::endl;
         exit(1);
         return false;
      }
   }
   return true;
}

bool BD::RunKernel(int start_tick, int B, int T, bool crowning_test, bool spotting_test, int stop_tick) {

   int counter = start_tick;

   InitialSetup<<<B,T>>>(g_fuel_t_, g_roth_data_, g_dead_sav_burnable_b_, g_dead_1h_b_, g_dead_10h_b_,
                         g_dead_100h_b_, g_live_h_b_, g_live_w_b_, g_fine_dead_extinctions_density_b_,
                         g_areas_reaction_factors_b_, g_slope_wind_factors_b_, g_residence_flux_live_sav_b_,
                         g_fuel_sav_accel_b_, g_slope_aspect_elevation_t_, g_dead_moistures_t_, g_live_moistures_t_,
                         g_x_wind, g_y_wind, sim_size_, g_DynLiveH, g_DynDead1);

   // UPDATE THE SPREAD
   UpdateSpread<<<B,T>>>(g_fuel_t_, g_roth_data_, g_fuel_sav_accel_b_, g_maxspreadrate_,
                         g_acceleration_constant_, g_intensity_modifier_, g_angles_, sim_size_);

   findMax<<<B,T,T*sizeof(int)>>>(g_maxspreadrate_, g_maxSpread_, g_CellSize_, g_timeStep_, sim_size_*16);

   cudaError_t err =  cudaMemcpy(&time_step_, g_timeStep_, sizeof(float), cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error cudaMemcpyDeviceToHost " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }
   std::cout<<"time_step_: "<<time_step_<<std::endl; // need to make this the true max right now if by chance cause
   // we need to do a blockwise reduction
   while(counter < stop_tick)
   {
      counter++;
      // if(counter >= 100)
      //    break;
      // std::cout<<"Here is counter: "<<counter<<std::endl;

      // BURN DISTANCE METHOD
      BurnDist<<<B,T>>>(g_toa_map_in_, g_toa_map_out_, g_loc_burndist_,
                        g_curspreadrate_, g_l_n_, sim_size_,
                        sim_rows_, sim_cols_, time_step_, current_time_);

      copyKernelBD<<<B,T>>>(g_toa_map_in_, g_toa_map_out_, sim_size_);

      if(crowning_test == 1 && spotting_test == 1)
      {
         TestCrownRate_WithSpot<<<B,T>>>(g_curspreadrate_, g_maxspreadrate_, g_intensity_modifier_, sim_size_,
                                         g_I_o_, g_RAC_, g_emberMap_, g_emberTestMap_, g_canopyHeight_, g_ember_list_,
                                         g_x_wind, g_y_wind, g_fuel_t_);

         TestSpotting<<<B,T>>>(g_emberMap_, g_canopyHeight_, g_ember_list_, current_time_, time_step_,
                               sim_size_, g_fuel_t_, g_DynDead1, g_dead_moistures_t_);

         updateTOA<<<B,T>>>(g_ember_list_, g_toa_map_in_, sim_size_);
      }
      else
      {
         TestCrownRate_NoSpot<<<B,T>>>(g_curspreadrate_, g_maxspreadrate_, g_intensity_modifier_,
                                       sim_size_, g_I_o_, g_RAC_);
      }

      Accelerate<<<B,T>>>(g_curspreadrate_, g_maxspreadrate_, sim_size_, time_step_);

      cudaDeviceSynchronize();

      current_time_ += time_step_;
   }
   return true;
}

bool BD::CopyFromDevice() {
   cudaError_t err =  cudaMemcpy(toa_map_, g_toa_map_in_, sim_size_*sizeof(float), cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error BD::CopyFromDevice() toa_map_, g_toa_map_in_: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   err =  cudaMemcpy(wind_x_map_, g_x_wind, sim_size_*sizeof(int), cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error BD::CopyFromDevice() wind_x_map_, g_x_wind: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   err =  cudaMemcpy(wind_y_map_, g_y_wind, sim_size_*sizeof(int), cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error BD::CopyFromDevice() wind_y_map_, g_y_wind: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   err =  cudaMemcpy(maxSpreadRates_, g_maxspreadrate_, sim_size_*16*sizeof(float), cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error BD::CopyFromDevice() maxSpreadRates_, g_maxspreadrate_: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   err =  cudaMemcpy(currentSpreadRates_, g_curspreadrate_, sim_size_*16*sizeof(float), cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error BD::CopyFromDevice() currentSpreadRates_, g_curspreadrate_: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   return true;
}

bool BD::WriteToFile(std::string filename)
{
   std::ofstream fout;
   fout.open(filename.c_str());
   int index = 0;

   for(unsigned int y = 0; y < simulation_->sim_dim_y_; y++)
   {
      for(unsigned int x = 0; x < simulation_->sim_dim_x_; x++)
      {
         if(x == simulation_->sim_dim_x_ - 1)
            fout << (int) toa_map_[index];
         else
            fout << (int) toa_map_[index] << ",";
         index++;
      }
      fout << '\n';
   }
   fout.close();
   return true;
}

bool BD::WritePauseToFile(std::string filename, std::string* metaptr, std::string msr, std::string csr)
{
   // Initialize fstream and index
   std::ofstream fout;
   int index = 0;

   // This writes the TOA to paused.csv
   fout.open(filename.c_str());

   for(int x = 0; x < 8; x++)
   {
      fout << metaptr[x];
   }
   fout << '\n';

   for(unsigned int y = 0; y < simulation_->sim_dim_y_; y++)
   {
      for(unsigned int x = 0; x < simulation_->sim_dim_x_; x++)
      {
         if((int) toa_map_[index] == 32767)
         {
            if(x == simulation_->sim_dim_x_ - 1)
               fout << 0;
            else
               fout << 0 << ",";
         }
         else
         {
            if(x == simulation_->sim_dim_x_ - 1)
               fout << 1;
            else
               fout << 1 << ",";
         }
         index++;
      }
      fout << '\n';
   }
   fout.close();

   // This writes the Max Spread Rates to msr.csv
   int dirs = 16;
   int cell = 0; // sim_x * sim_y * 16
   fout.open(msr.c_str());

   for(int row = 0; row < simulation_->sim_dim_y_; row++)
   {
      for (int col = 0; col < simulation_->sim_dim_x_-1; col++)
      {
         for (unsigned int i = 0; i < dirs; i++, cell++)
         {
            fout << maxSpreadRates_[cell] << ",";
         }
      }
      for (unsigned int i = 0; i < dirs-1; i++, cell++)
      {
         fout << maxSpreadRates_[cell] << ",";
      }
      fout << 0;
      fout << '\n';
   }
   fout.close();

   // This writes the Current Spread Rates to csr.csv

   cell = 0; // sim_x * sim_y * 16
   fout.open(csr.c_str());

   for(int row = 0; row < simulation_->sim_dim_y_; row++)
   {
      for (int col = 0; col < simulation_->sim_dim_x_-1; col++)
      {
         for (unsigned int i = 0; i < dirs; i++, cell++)
         {
            fout << currentSpreadRates_[cell] << ",";
         }
      }
      for (unsigned int i = 0; i < dirs-1; i++, cell++)
      {
         fout << currentSpreadRates_[cell] << ",";
      }
      fout << 0;
      fout << '\n';
   }
   fout.close();
   return true;
}

bool BD::UpdateCell(int _x, int _y, int val)
{
   int cell = _x * sim_cols_ + _y;
   if(cell < 0 || cell > sim_size_)
   {
      return false;
   }
   toa_map_[cell] = val;
   // std::cout<<toa_map_[cell]<<std::endl;
   return true;
}