//
// Created by jsmith on 10/12/15.
//
//#include <driver_types.h>
#include "IMT.h"

IMT::IMT(int _x, int _y, std::string fuel_model_name, std::string fuel_moisture_name)
        : Propagation(_x, _y, fuel_model_name, fuel_moisture_name){
   toa_map_ = (int*)malloc(sim_size_ * sizeof(int));
   timesteppers_ = (float*)malloc(2 * sizeof(float));
   check_ = (bool*)malloc(sim_size_ * sizeof(bool));
//   l_n_ = (float*)malloc(16 * sizeof(float));
}

IMT::~IMT(){
   // Free Host Memory
   free(toa_map_);
   free(timesteppers_);
   free(check_);
//   free(l_n_);
   // Free Device Memory
   cudaFree(g_toa_map_now_);
   cudaFree(g_toa_map_next_);
   cudaFree(g_toa_map_step_);
   cudaFree(g_timesteppers_);
   cudaFree(g_l_n_);
}

bool IMT::Init(std::string fuel_file, std::string terrain_file,
               std::string canopy_height_file, std::string crown_base_height_file,
               std::string crown_bulk_density_file, std::string wind_x, std::string wind_y,
               std::string dynLiveH, std::string dynDead1, std::string csr, std::string msr,
               int M_FLAG, int C_FLAG, int S_FLAG) {
   Propagation::Init(fuel_file, terrain_file,
                     canopy_height_file, crown_base_height_file,
                     crown_bulk_density_file, wind_x, wind_y, dynLiveH,
                     dynDead1, csr, msr, M_FLAG, C_FLAG, S_FLAG);
   // Initialize TOA Map
   for(unsigned int i = 0; i < sim_size_; i++){
      toa_map_[i] = simulation_->ign_time_[i];
      check_[i] = false;
   }
   // Initialize TimeNow and TimeNext
   timesteppers_[1] = timesteppers_[0] = 0;
   // Populate lengths
//   for(unsigned int i = 0; i < 16; i++){
//      l_n_[i] = simulation_->l_n_[i];
//   }
   l_n_ = simulation_->l_n_;
   return true;
}

bool IMT::CopyToDevice(int M_FLAG, int C_FLAG, int S_FLAG) {
   Propagation::CopyToDevice(M_FLAG, C_FLAG, S_FLAG);
   // Create Memory on Device
   cudaError_t err;
   err = cudaMalloc((void**) &g_toa_map_now_, sim_size_*sizeof(int));
   err = cudaMalloc((void**) &g_toa_map_next_, sim_size_*sizeof(int));
   err = cudaMalloc((void**) &g_toa_map_step_, sim_size_*sizeof(int));
   err = cudaMalloc( (void**) &g_wind_x_map_in_, sim_size_*sizeof(float));
   err = cudaMalloc( (void**) &g_wind_x_map_out_, sim_size_*sizeof(float));
   err = cudaMalloc( (void**) &g_wind_y_map_in_, sim_size_*sizeof(float));
   err = cudaMalloc( (void**) &g_wind_y_map_out_, sim_size_*sizeof(float));
   err = cudaMalloc((void**) &g_timesteppers_, 2*sizeof(float));
   err = cudaMalloc((void**) &g_l_n_, 16 * sizeof(float));
   err = cudaMalloc((void**) &g_check_, sim_size_ * sizeof(bool));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating Memory in IMT Class: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   // Copy data to device
   err = cudaMemcpy(g_toa_map_now_, toa_map_, sim_size_*sizeof(int), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_toa_map_next_, toa_map_, sim_size_*sizeof(int), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_toa_map_step_, toa_map_, sim_size_*sizeof(int), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_wind_x_map_in_, wind_x_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_wind_x_map_out_, wind_x_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_wind_y_map_in_, wind_y_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_wind_y_map_out_, wind_y_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_timesteppers_, timesteppers_, 2 * sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_l_n_, l_n_, 16*sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_check_, check_, sim_size_*sizeof(bool), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying Memory in IMT Class: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }
   return true;
}

bool IMT::RunKernel(int start_tick, int B, int T, bool crowning_flag, bool spotting_flag, int stop_tick) {
   int counter = start_tick;
   int terminate = -1;
//   while(terminate <= 0){
   while(counter < stop_tick){
      counter++;
      // ITERATIVE MINIMAL TIME
      // Do calculations
      ItMinTime<<<B,T>>>(g_toa_map_now_,g_toa_map_next_, g_toa_map_step_, g_curspreadrate_,
                         g_timesteppers_, g_l_n_, g_check_, sim_size_,
                         sim_rows_, sim_cols_);
      // cout << "step caclulated\n";
      // Copy from output to write
//      copyKernelIMT<<<B,T>>>(g_toa_map_now_, g_toa_map_step_,
//                             g_check_, sim_size_);
      if(crowning_flag)
         TestCrownRate_NoSpot<<<B,T>>>(g_curspreadrate_, g_maxspreadrate_, g_intensity_modifier_, sim_size_, g_I_o_, g_RAC_);
      // Accelerate Fire
      Accelerate<<<B,T>>>(g_curspreadrate_, g_maxspreadrate_, sim_size_ * 16, simulation_->time_step_);

      cudaDeviceSynchronize();

      if(terminate < sim_size_)
         terminate = -1;

      // cout << counter <<endl;
      // Swap Pointers for loop
      int *swap = g_toa_map_now_;
      g_toa_map_now_ = g_toa_map_next_;
      g_toa_map_next_ = swap;

   }
   return true;
}

bool IMT::CopyFromDevice() {
   cudaError_t err = cudaMemcpy(toa_map_, g_toa_map_now_, sim_size_ * sizeof(int), cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error copying data from GPU: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   printf("IMT Copy WINDX From Device\n");
   err =  cudaMemcpy(wind_x_map_, g_wind_x_map_in_,
                                 sim_size_*sizeof(float),
                                 cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error copying from GPU: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   printf("IMT Copy WINDY From Device\n");
   err =  cudaMemcpy(wind_y_map_, g_wind_y_map_in_,
                                 sim_size_*sizeof(float),
                                 cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error copying from GPU: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   return true;
}

bool IMT::WriteToFile(std::string filename) {
   std::ofstream fout;
//   std::string filename;
//   filename += simulation_->root_path_;
//   filename += "out/IMT_test.csv";
   fout.open(filename.c_str());
   for(unsigned int i = 0; i < sim_size_; i++){
      if(i % simulation_->sim_dim_x_ == 0 && i !=0){
         fout << '\n';
      }
      fout << (int) toa_map_[i] << ",";
   }
   fout.close();
   return true;
}

bool IMT::WindXToFile(std::string filename, std::string* metaptr) {
   std::ofstream fout;
   fout.open(filename.c_str());

   // add metadata to the first eight lines
   for(int x = 0; x < 8; x++)
   {
      fout << metaptr[x];
   }
   fout << '\n';

   for(unsigned int i = 0; i < sim_size_; i++){
      if(i % simulation_->sim_dim_x_ == 0 && i !=0){
         fout << '\n';
      }
      fout << (float) wind_x_map_[i] << ",";
   }
   fout.close();
   return true;
}

bool IMT::WindYToFile(std::string filename, std::string* metaptr) {
   std::ofstream fout;
   fout.open(filename.c_str());

   // add metadata to the first eight lines
   for(int x = 0; x < 8; x++)
   {
      fout << metaptr[x];
   }
   fout << '\n';

   for(unsigned int i = 0; i < sim_size_; i++){
      if(i % simulation_->sim_dim_x_ == 0 && i !=0){
         fout << '\n';
      }
         fout << (float) wind_y_map_[i] << ",";
   }
   fout.close();
   return true;
}

bool IMT::UpdateCell(int _x, int _y, int val){
//   if(_x < 0 || _y < 0 || _x > sim_rows_ || _y > sim_cols_)
//      return false;
   int cell = _x * sim_cols_ + _y;
   toa_map_[cell] = val;
   return true;
}
