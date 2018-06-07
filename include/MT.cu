//#include <driver_types.h>
#include "MT.h"

MT::MT(int _x, int _y, std::string fuel_model_name, std::string fuel_moisture_name)
            : Propagation(_x, _y, fuel_model_name, fuel_moisture_name) {
   toa_map_ = (int*) malloc(sim_size_ * sizeof(int));
   timesteppers_ = (int*) malloc(2* sizeof(int));
//   l_n_ = (float*) malloc(16 * sizeof(float));
}

MT::~MT(){
   // Free Host Memory
   free(toa_map_);
   free(timesteppers_);
   free(l_n_);
   // Free Device Memory
   cudaFree(g_toa_map_);
   cudaFree(g_timesteppers_);
   cudaFree(g_l_n_);
}

bool MT::Init(std::string fuel_file, std::string terrain_file,
              std::string canopy_height_file, std::string crown_base_height_file,
              std::string crown_bulk_density_file, std::string wind_x, std::string wind_y,
              std::string dynLiveH, std::string dynDead1, std::string csr, std::string msr,
              int M_FLAG, int C_FLAG, int S_FLAG) {
   // Call Parent Init Fcn
   Propagation::Init(fuel_file, terrain_file,
           canopy_height_file, crown_base_height_file,
           crown_bulk_density_file, wind_x, wind_y, dynLiveH,
           dynDead1, csr, msr, M_FLAG, C_FLAG, S_FLAG);
   // Initialize TOA Map
   for(unsigned int i = 0; i < sim_size_; i++){
      toa_map_[i] = simulation_->ign_time_[i];
   }
   // Initialize TimeNow and TimeNext
   timesteppers_[0] = timesteppers_[1] = 0;
   // Populate lengths
//   for(unsigned int i = 0; i < 16; i++){
//      l_n_[i] = simulation_->l_n_[i];
//   }
   l_n_ = simulation_->l_n_;
   return true;
}

bool MT::CopyToDevice(int M_FLAG, int C_FLAG, int S_FLAG) {
   Propagation::CopyToDevice(M_FLAG, C_FLAG, S_FLAG);
   // Create memory on device
   cudaError_t err = cudaMalloc((void**) &g_toa_map_, sim_size_*sizeof(int));
   err = cudaMalloc( (void**) &g_wind_x_map_in_, sim_size_*sizeof(float));
   err = cudaMalloc( (void**) &g_wind_x_map_out_, sim_size_*sizeof(float));
   err = cudaMalloc( (void**) &g_wind_y_map_in_, sim_size_*sizeof(float));
   err = cudaMalloc( (void**) &g_wind_y_map_out_, sim_size_*sizeof(float));
   err = cudaMalloc((void**) &g_timesteppers_, 2*sizeof(int));
   err = cudaMalloc((void**) &g_l_n_, 16*sizeof(float));
   if (err != cudaSuccess) {
      std::cerr << "Error Allocating Memory in MT Class: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   // Copy data to device
   err = cudaMemcpy(g_toa_map_, toa_map_, sim_size_*sizeof(int), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_wind_x_map_in_, wind_x_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_wind_x_map_out_, wind_x_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_wind_y_map_in_, wind_y_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_wind_y_map_out_, wind_y_map_, sim_size_*sizeof(float), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_timesteppers_, timesteppers_, 2*sizeof(int), cudaMemcpyHostToDevice);
   err = cudaMemcpy(g_l_n_, l_n_, 16* sizeof(float), cudaMemcpyHostToDevice);
   if (err != cudaSuccess) {
      std::cerr << "Error Copying Memory in MT Class: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }
   return true;
}

bool MT::RunKernel(int start_tick, int B, int T, bool crowning_flag, bool spotting_flag, int stop_tick) {
//   printf("Kicking off Kernels\n");
   int counter = start_tick;
   //int terminate = -1;
//   int B = 64;
//   int T = 64;
   //   while(terminate <= 0){
   while(counter < stop_tick){
      counter++;
      // Do calculations
      if(crowning_flag)
         TestCrownRate_NoSpot<<<B,T>>>(g_curspreadrate_, g_maxspreadrate_, g_intensity_modifier_, sim_size_, g_I_o_, g_RAC_);
      // Accelerate Fire
      Accelerate<<<B,T>>>(g_curspreadrate_,g_maxspreadrate_, sim_size_*16, simulation_->time_step_);
//      MinTime<<<B,T>>>(g_toa_map_, g_maxspreadrate_,
      MinTime<<<B,T>>>(g_toa_map_, g_curspreadrate_,
                       g_timesteppers_, g_l_n_, sim_size_,
                       sim_rows_, sim_cols_);
      // Update Time Kernel
      TimeKernelMT<<<1,1>>>(g_timesteppers_);


      cudaDeviceSynchronize();
//      cudaError_t err = cudaMemcpyFromSymbol(&terminate, end, sizeof(end), 0,
//                              cudaMemcpyDeviceToHost);
//      if (err != cudaSuccess) {
//         std::cerr << "Error copying from GPU: " << cudaGetErrorString(err) << std::endl;
//         exit(1);
//         return false;
//      }
//      printf("end: %d\n", terminate);
//      if(terminate < 4)
//      terminate = -1;
//   }
//   int finishCount = 0;
//   // Catch last corner to terminate simulation
//   while(finishCount <= 3){
//      counter++;
//      finishCount++;
//      // Do calculations
//      MinTime<<<B,T>>>(g_toa_map_, g_maxspreadrate_,
//                       g_timesteppers_, g_l_n_, sim_size_,
//                       sim_rows_, sim_cols_);
//      // Update Time Kernel
//      timeKernelMT<<<1,1>>>(g_timesteppers_);
   }
   return true;
}

bool MT::CopyFromDevice() {
   cudaError_t err = cudaMemcpy(toa_map_, g_toa_map_, sim_size_ * sizeof(int), cudaMemcpyDeviceToHost);
   if(err != cudaSuccess){
      std::cerr << "Error copying from GPU: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   printf("MT Copy WINDX From Device\n");
   err =  cudaMemcpy(wind_x_map_, g_wind_x_map_in_,
                                 sim_size_*sizeof(float),
                                 cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) {
      std::cerr << "Error copying from GPU: " << cudaGetErrorString(err) << std::endl;
      exit(1);
      return false;
   }

   printf("MT Copy WINDY From Device\n");
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

bool MT::WriteToFile(std::string filename) {
   std::ofstream fout;
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

bool MT::WindXToFile(std::string filename, std::string* metaptr) {
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

bool MT::WindYToFile(std::string filename, std::string* metaptr) {
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

bool MT::UpdateCell(int _x, int _y, int val){
   if(_x < 0 || _y < 0 || _x > sim_rows_ || _y > sim_cols_)
      return false;
   int cell = _x * sim_cols_ + _y;
   toa_map_[cell] = val;
   return true;
}
