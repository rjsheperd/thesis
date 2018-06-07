#include "FireSim.h"
#include <algorithm>
#include <cmath>
#include <typeinfo>

// these are used to figure out moisture
// ./simulator /home/andyh/Desktop/spotting/data/fixed.fuel /home/andyh/Desktop/spotting/data/fire_info.csv /home/andyh/Desktop/spotting/out/final_tests.csv
#include <fstream>

const int INF = 32767;

float* GISToFloatArray(const char*, int, int, float&);
int* GISToIntArray(const char*, int, int, float&);

/*
Constructor: builds simplest test case for testing code
*/
FireSim::FireSim(int _x, int _y, std::string fuel_model_name, std::string fuel_moistures_name)
{
  // sim dim x and y
  sim_dim_x_ = _x;
  sim_dim_y_ = _y;

  // read in the fuel model properties file, this file is table 7 from the farsite.pdf
  // it is the properties such as fuel bed depth, and how hot types of vegetation burn
  const char * fname_tmp = fuel_model_name.c_str();
  _models = sim::readFuelModels(fname_tmp);

  // this is hand made I assume because landfire.gov doesn't have any files that match
  // I believe it was made to best represent the vegetation in kyle canyon, refer to thesis
  // for a more in depth description
  fname_tmp = fuel_moistures_name.c_str();
  _moistures = sim::readFuelMoistures(fname_tmp);

  numModels = _models.size();
  numMoistModels = _moistures.size();
  foliar_moisture = 1.0f; // default, literally all literature assumes 1.0

  // roth_data is used to calculate max spread rate for every cell
  roth_data_ = new float4[sim_dim_x_ * sim_dim_y_];

  // angles are needed for all cells in 16 directions, this probably could have been a
  // size 16 array and just looped repeatedly, but this made more since to me at the time
  angles_ = new float[sim_dim_x_ * sim_dim_y_];

  // rothermel properties and values most of these are calculated using fire science equations
  // or they are the properties of the terrain and or fuel model for the corresponding cell
  slope_aspect_elevation_t_ = new float3[sim_dim_x_ * sim_dim_y_];
  dead_sav_burnable_b_ = new float4[numModels];
  dead_1h_b_ = new float4[numModels];
  dead_10h_b_ = new float4[numModels];
  dead_100h_b_ = new float4[numModels];
  live_h_b_ = new float4[numModels];
  live_w_b_ = new float4[numModels];
  fine_dead_extinctions_density_b_ = new float4[numModels];
  areas_reaction_factors_b_ = new float4[numModels];
  slope_wind_factors_b_ = new float4[numModels];
  residence_flux_live_sav_b_ = new float4[numModels];
  fuel_sav_accel_b_ = new float2[numModels];
  dead_moistures_t_ = new float3[numMoistModels];
  live_moistures_t_ = new float2[numMoistModels];
  I_o_ = new float[sim_dim_x_ * sim_dim_y_];
  RAC_ = new float[sim_dim_x_ * sim_dim_y_];
  canopy_height_ = new float[sim_dim_x_ * sim_dim_y_];
  acceleration_constant_ = 1.0; // default

  // ignite time old and ignite time new
  ign_time_ = new float[sim_dim_x_ * sim_dim_y_];
  ign_time_new_ = new float[sim_dim_x_ * sim_dim_y_];

  // current and max spread rates
  csr_ = new float[sim_dim_x_ * sim_dim_y_ * 16];
  msr_ = new float[sim_dim_x_ * sim_dim_y_ * 16];

  // burn distances between cells, used to see how much distance is present and left to burn
  // help with determining if a nieghbor cell is to be ignited yet
  burn_dist_ = new float*[sim_dim_x_ * sim_dim_y_];
  for(unsigned int i = 0; i < sim_dim_x_ * sim_dim_y_; i++){
    burn_dist_[i] = new float[8];
  }
  l_n_ = new float[16]; // this helps interpolate and populate the burn lengths/distances
}

/*
Destructor: builds simplest test case for testing code
*/
FireSim::~FireSim(){
   // Clean everything up
   delete roth_data_;
   delete slope_aspect_elevation_t_;
   delete dead_sav_burnable_b_;
   delete dead_1h_b_;
   delete dead_10h_b_;
   delete dead_100h_b_;
   delete live_h_b_;
   delete live_w_b_;
   delete fine_dead_extinctions_density_b_;
   delete areas_reaction_factors_b_;
   delete slope_wind_factors_b_;
   delete residence_flux_live_sav_b_;
   delete fuel_sav_accel_b_;
   delete dead_moistures_t_;
   delete live_moistures_t_;
   delete angles_;
   delete I_o_;
   delete RAC_;
   delete canopy_height_;
   delete csr_;
   delete msr_;
   acceleration_constant_ = 1.0;
   delete ign_time_;
   delete ign_time_new_;
   for(unsigned int i = 0; i < sim_dim_x_ * sim_dim_y_; i++){
      delete burn_dist_[i];
   }
   delete burn_dist_;
}

/*
Shader base: rothermel
Purpose: Initializes the sim.
*/
void FireSim::Init(std::string fuel_file, std::string terrain_file,
                   std::string canopy_height_file, std::string crown_base_height_file,
                   std::string crown_bulk_density_file, std::string wind_x, std::string wind_y, std::string dynLiveH, std::string dynDead1,
                   std::string csr, std::string msr, int M_FLAG, int C_FLAG, int S_FLAG){

  // These are the 16 directions
  float angles[16] = {  90,  45,  0, -45, -90, -115,  180,  115, 63.4349, 26.5651, -26.5651f,
                       -63.4349f, -116.5651f, -153.4349f, 153.4349f, 116.5651f};
                       
  int cell = 0; // used to interate
  float junk;   // used as a null argument to stop compile errors
  int* slopeTexTmp = NULL, *tempRAC = NULL, *tempIo = NULL, *tempCH = NULL; // used because GISTOFLOATARRAY isn't working
  GDALAllRegister(); // used to initialize GDAL extraction

  // this is the vegetation type map for the cells, junk is used because we don't extract cell size from this file
  const char * fname_tmp = fuel_file.c_str();
  fuel_t_ = GISToIntArray(fname_tmp, sim_dim_x_, sim_dim_y_, junk);

  // this is the terrain tif, cell_size_ is used here because this is where cell size is found
  fname_tmp = terrain_file.c_str();
  slopeTexTmp = GISToIntArray(fname_tmp, sim_dim_x_ *3, sim_dim_y_ *3, cell_size_);

  if(C_FLAG == 1) // if crowning is on we need this
  {
    // canopy height file
    fname_tmp = canopy_height_file.c_str();
    tempCH = GISToIntArray(fname_tmp, sim_dim_x_, sim_dim_y_, junk);

    for(unsigned int i = 0; i < sim_dim_x_*sim_dim_y_; i++)
    {
      // Canopy Height
      canopy_height_[i] = (float) tempCH[i] / 100; // GISTOFLOATARRAY not working but the GISTOINTARRAY was
      // That is why there are temp pointers and why we are casting ints to floats divided by 100
    }
  }

  // wind in the x-axis file
  fname_tmp = wind_x.c_str();
  xwind_ = GISToIntArray(fname_tmp, sim_dim_x_, sim_dim_y_, junk);

  // wind in the y-axis file
  fname_tmp = wind_y.c_str();
  ywind_ = GISToIntArray(fname_tmp, sim_dim_x_, sim_dim_y_, junk);

  if(M_FLAG == 1)
  {
    // dynamic moisture is on we need these files
    fname_tmp = dynLiveH.c_str();
    DynLiveH = GISToIntArray(fname_tmp, sim_dim_x_, sim_dim_y_, junk);

    fname_tmp = dynDead1.c_str();
    DynDead1 = GISToIntArray(fname_tmp, sim_dim_x_, sim_dim_y_, junk);    
  }

  std::ifstream fin(csr);

  if(!fin)
  {
   std::cerr<<"Failed to open CSR file !";
   exit(0);
  }

  int r = 0, c = 0;
  while(r < sim_dim_y_)
  {
   int value;
   char dummy;
   while(c < sim_dim_x_-1)
   {      
     fin >> value >> dummy;
     for(int i = 0; i < 16; i++)
     {
        csr_[cell+i] = value;
     }
     c++;
     cell++;
   }
   fin >> value;
   for(int i = 0; i < 16; i++)
   {
     csr_[cell+i] = value;
   }
   c = 0;
   r++;
   cell++;
  }
  fin.close();

  cell = 0;

  fin.open(msr);

  if(!fin)
  {
   std::cerr<<"Failed to open MSR file !";
   exit(0);
  }

  r = 0;
  c = 0;
  while(r < sim_dim_y_)
  {
   int value;
   char dummy;
   while(c < sim_dim_x_-1)
   {      
     fin >> value >> dummy;
     for(int i = 0; i < 16; i++)
     {
        msr_[cell+i] = value;
     }
     c++;
     cell++;
   }
   fin >> value;
   for(int i = 0; i < 16; i++)
   {
     msr_[cell+i] = value;
   }
   c = 0;
   r++;
   cell++;
  }
  fin.close();

  // reset our index
  cell = 0;

  if(C_FLAG == 1) // Don't need this if crowning is not on
  {
    // Crown base height data is only used for calculating I_o so it doesn't need to be stored separately
    // I_o is stored and transferred to the GPU
    fname_tmp = crown_base_height_file.c_str();
    tempIo = GISToIntArray(fname_tmp, sim_dim_x_, sim_dim_y_, junk);
    for(unsigned int i = 0; i < sim_dim_x_*sim_dim_y_; i++){
      I_o_[i] = (float) pow(0.01 * tempIo[i] * (460.0 + 25.9 * foliar_moisture), 1.5);
      if(I_o_[i] == 0) // sometimes there's a 0 in the file, we think this can either represent a very high value or no value
        I_o_[i] = (float) pow(0.01 * 100 * (460.0 + 25.9 * foliar_moisture), 1.5);
      // tempIo's value is changed to 100 because 100 is the max value in crown base height file, we assume high value here
      // maybe we are wrong and maybe it is a test case to tell the sim that there is no canopy here maybe check against
      // the 3 crown files to see whether the zereos are in correlation
    }

    // Crown Bulk Density data is only used for calculating RAC so it doesn't need to be stored separately
    fname_tmp = crown_bulk_density_file.c_str();
    tempRAC = GISToIntArray(fname_tmp, sim_dim_x_, sim_dim_y_, junk);
    for(unsigned int i = 0; i < sim_dim_x_*sim_dim_y_; i++){
      if(tempRAC[i] != 0)
      {
         RAC_[i] = (float) 3.0f / tempRAC[i]; // we assume the same thing here as above, high value is assumed
      }                                       // kind of have to be consistent with the assumptions
      else
      {
         RAC_[i] = 3.f/30.f; // 30 because this is the max number from crown bulk density file
      }
    }
  }

  for(unsigned int i = 0; i < sim_dim_y_; i++)
  {
    for(int j = 0; j < sim_dim_x_; j++, cell++)
    {
      roth_data_[cell].w = -1.0f; // was a value for testing purposes, will be overwritten
      roth_data_[cell].x = roth_data_[cell].y = roth_data_[cell].z = 0.f; // default value

      // Rothermel Data Members
      slope_aspect_elevation_t_[cell].x = (float) slopeTexTmp[3*cell] / 100000; // honestly no idea why I am divinding by 100,000
      slope_aspect_elevation_t_[cell].y = (float) slopeTexTmp[3*cell+1] / 100000; // the GISFLOATARRAY was giving me -1 * e^100000
      slope_aspect_elevation_t_[cell].z = (float) slopeTexTmp[3*cell+2] / 100000; // which is literally 0 so this was what I needed
      // to do to make it work. There may be more to this but for now it is uniform

      // initializing ignite times, burn distances, and at the very end we are creating the angles array
      ign_time_[cell] = INF;
      ign_time_new_[cell] = INF;
      for(int k = 0; k < 8; k++)
      {
        burn_dist_[cell][k] = l_n_[k];
      }
      angles_[cell] = angles[cell%16]; // an angle for every cell x*y
    }
  }

  // so this is storing the fuel model properties from the rothermel and anderson
  // fuel model that were created by fancy fire scientists and used for a lot of preprocessing
  unsigned int i = 0;
  for (std::vector<sim::FuelModel>::iterator it = _models.begin();
      it != _models.end(); it++)
  {
    dead_1h_b_[i].x = it->effectiveHeatingNumber[sim::Dead1h];
    dead_1h_b_[i].y = it->load[sim::Dead1h];
    dead_1h_b_[i].z = it->areaWeightingFactor[sim::Dead1h];
    dead_1h_b_[i].w = it->fuelMoisture[sim::Dead1h];

    dead_10h_b_[i].x = it->effectiveHeatingNumber[sim::Dead10h];
    dead_10h_b_[i].y = it->load[sim::Dead10h];
    dead_10h_b_[i].z = it->areaWeightingFactor[sim::Dead10h];
    dead_10h_b_[i].w = it->fuelMoisture[sim::Dead10h];

    dead_100h_b_[i].x = it->effectiveHeatingNumber[sim::Dead100h];
    dead_100h_b_[i].y = it->load[sim::Dead100h];
    dead_100h_b_[i].z = it->areaWeightingFactor[sim::Dead100h];
    dead_100h_b_[i].w = it->fuelMoisture[sim::Dead100h];

    live_h_b_[i].x = it->effectiveHeatingNumber[sim::LiveH];
    live_h_b_[i].y = it->load[sim::LiveH];
    live_h_b_[i].z = it->areaWeightingFactor[sim::LiveH];
    live_h_b_[i].w = it->fuelMoisture[sim::LiveH];

    live_w_b_[i].x = it->effectiveHeatingNumber[sim::LiveW];
    live_w_b_[i].y = it->load[sim::LiveW];
    live_w_b_[i].z = it->areaWeightingFactor[sim::LiveW];
    live_w_b_[i].w = it->fuelMoisture[sim::LiveW];

    fine_dead_extinctions_density_b_[i].x = it->fineDeadRatio;
    fine_dead_extinctions_density_b_[i].y = it->extinctionMoisture;
    fine_dead_extinctions_density_b_[i].z = it->liveExtinction;
    fine_dead_extinctions_density_b_[i].w = it->fuelDensity;

    areas_reaction_factors_b_[i].x = it->deadArea;
    areas_reaction_factors_b_[i].y = it->liveArea;
    areas_reaction_factors_b_[i].z = it->deadReactionFactor;
    areas_reaction_factors_b_[i].w = it->liveReactionFactor;

    slope_wind_factors_b_[i].x = it->slopeK;
    slope_wind_factors_b_[i].y = it->windK;
    slope_wind_factors_b_[i].z = it->windB;
    slope_wind_factors_b_[i].w = it->windE;

    residence_flux_live_sav_b_[i].x = it->residenceTime;
    residence_flux_live_sav_b_[i].y = it->propagatingFlux;
    residence_flux_live_sav_b_[i].z = it->SAV[sim::LiveH];
    residence_flux_live_sav_b_[i].w = it->SAV[sim::LiveW];

    dead_sav_burnable_b_[i].x = it->SAV[sim::Dead1h];
    dead_sav_burnable_b_[i].y = it->SAV[sim::Dead10h];
    dead_sav_burnable_b_[i].z = it->SAV[sim::Dead100h];
    dead_sav_burnable_b_[i].w = 100.0f;

    fuel_sav_accel_b_[i].x = it->fuelSAV;
    fuel_sav_accel_b_[i].y = it->accelerationConstant;
    i++;
  }

  // So, technically each fuel model has a average amount of moisture in it all year around that changes
  // but we use these averages from farsite and fire scientists to get a default state at which a place will
  // burn, then we added our dynamic moisture to enable water drops, rain, weather changes, but that needs
  // more research done cause we assume things, refer to thesis for more explanation
  i = 0;
  for (std::vector<sim::FuelMoisture>::iterator it = _moistures.begin();
      it != _moistures.end(); it++, i++)
  {
    dead_moistures_t_[i].x = it->dead1h;
    dead_moistures_t_[i].y = it->dead10h;
    dead_moistures_t_[i].z = it->dead100h;

    live_moistures_t_[i].x = it->liveH;
    live_moistures_t_[i].y = it->liveW;
  }

  // this is where we interpolate the burn distances
  float orthoSize = cell_size_;
  float diagSize = (float) (cell_size_ * sqrt(2));
  float superSize = (float) sqrt(pow(cell_size_, 2) + pow(cell_size_ * 2, 2));
  static float L_n_tmp[16] =  { orthoSize, diagSize, orthoSize, diagSize, orthoSize, diagSize,
                               orthoSize, diagSize, superSize, superSize, superSize, superSize,
                               superSize, superSize, superSize, superSize};
  for(int j = 0; j < 16; j++){
    l_n_[j] = L_n_tmp[j];
  }
}

/*
Function: UpdateSpreadData
Input: The necessary inputs are the values that are found in the textures/buffers
       in the FuelModel.h/Simulator.cpp files in Roger's code
Shader base: rothermel
Purpose: This runs rothermel's equations to initialize simulation; this is now only used for the sequential
because the preprocessing for the parallel was significant enough to move to the GPU and cut down time
*/
void FireSim::UpdateSpreadData(){
   
   std::cout << "Updating Spread Data . . ." << std::endl;
   int cell = 0;
   float dist = 10.f;

   for(unsigned int i = 0; i < sim_dim_y_; i++){
      for(int j = 0; j < sim_dim_x_; j++, cell++){

         int fuelModel = fuel_t_[cell];

         if(fuelModel == 99)
         {
            roth_data_[cell].w = 1.0f; // max spread rate
            roth_data_[cell].x = 0.0f; // ellipse eccentricity
            roth_data_[cell].y = 1.0f; // spread direction
            roth_data_[cell].z = 1.0f; // intensity modifier
            continue;
         }

         if(dead_sav_burnable_b_[fuelModel].w < 50.0){
            std::cout << "Warning: Something may have gone wrong. Check that Files were read Correctly." << std::endl;
            continue;
         }

         float maxSpreadRate = 0.f;
         float ellipseEccentricity = 0.f;
         float spreadDirection = 0.f;
         float spreadModifier = 0.f;
         float3 timeLagClass;

         if (dead_sav_burnable_b_[fuelModel].x > 192.0){
            timeLagClass.x = dead_moistures_t_[fuelModel].x;
         }
         else if (dead_sav_burnable_b_[fuelModel].x > 48.0){
            timeLagClass.x = dead_moistures_t_[fuelModel].y;
         }
         else{
            timeLagClass.x = dead_moistures_t_[fuelModel].z;
         }

         if (dead_sav_burnable_b_[fuelModel].y > 192.0){
            timeLagClass.y = dead_moistures_t_[fuelModel].x;
         }
         else if (dead_sav_burnable_b_[fuelModel].y > 48.0){
            timeLagClass.y = dead_moistures_t_[fuelModel].y;
         }
         else{
            timeLagClass.y = dead_moistures_t_[fuelModel].z;
         }

         if (dead_sav_burnable_b_[fuelModel].z > 192.0){
            timeLagClass.z = dead_moistures_t_[fuelModel].x;
         }
         else if (dead_sav_burnable_b_[fuelModel].z > 48.0){
            timeLagClass.z = dead_moistures_t_[fuelModel].y;
         }
         else{
            timeLagClass.z = dead_moistures_t_[fuelModel].z;
         }

         float weightedFuelModel =
               timeLagClass.x * dead_1h_b_[fuelModel].x * dead_1h_b_[fuelModel].y +
               timeLagClass.y * dead_10h_b_[fuelModel].x * dead_10h_b_[fuelModel].y +
               timeLagClass.z * dead_100h_b_[fuelModel].x * dead_100h_b_[fuelModel].y;

         float fuelMoistures[5];
         fuelMoistures[0] = timeLagClass.x;
         fuelMoistures[1] = timeLagClass.y;
         fuelMoistures[2] = timeLagClass.z;
         fuelMoistures[3] = live_moistures_t_[fuelModel].x;
         fuelMoistures[4] = live_moistures_t_[fuelModel].y;

         float liveExtinction = 0.0;
         if(live_h_b_[fuelModel].y > 0.0 || live_w_b_[fuelModel].y > 0.0){
            float fineDeadMoisture = 0.0;
            if (fine_dead_extinctions_density_b_[fuelModel].x > 0.0){
               fineDeadMoisture = weightedFuelModel / fine_dead_extinctions_density_b_[fuelModel].x;
            }

            liveExtinction =
                  (fine_dead_extinctions_density_b_[fuelModel].z *
                   (1.0f - fineDeadMoisture / fine_dead_extinctions_density_b_[fuelModel].y)) - 0.226f;
            liveExtinction = std::max(liveExtinction, fine_dead_extinctions_density_b_[fuelModel].y);
         }

         float heatOfIgnition =
               areas_reaction_factors_b_[fuelModel].x *
               ((250.0f + 1116.0f * fuelMoistures[0]) * dead_1h_b_[fuelModel].z * dead_1h_b_[fuelModel].x +
                (250.0f + 1116.0f * fuelMoistures[1]) * dead_10h_b_[fuelModel].z * dead_10h_b_[fuelModel].x +
                (250.0f + 1116.0f * fuelMoistures[2]) * dead_100h_b_[fuelModel].z * dead_100h_b_[fuelModel].x) +
               areas_reaction_factors_b_[fuelModel].y *
               ((250.0f + 1116.0f * fuelMoistures[3]) * live_h_b_[fuelModel].z * live_h_b_[fuelModel].x +
                (250.0f + 1116.0f * fuelMoistures[4]) * live_w_b_[fuelModel].z * live_w_b_[fuelModel].x);
         heatOfIgnition *= fine_dead_extinctions_density_b_[fuelModel].w;

         float liveMoisture = live_h_b_[fuelModel].z * fuelMoistures[3] + live_w_b_[fuelModel].z * fuelMoistures[4];
         float deadMoisture = dead_1h_b_[fuelModel].z * fuelMoistures[0] +
                              dead_10h_b_[fuelModel].z * fuelMoistures[1] +
                              dead_100h_b_[fuelModel].z * fuelMoistures[2];

         float reactionIntensity = 0.0;

         if (liveExtinction > 0.0)
         {
            float r = liveMoisture / liveExtinction;
            if (r < 1.0){
               reactionIntensity += areas_reaction_factors_b_[fuelModel].w * (1.0 -
                                                              (2.59 * r) +
                                                              (5.11 * r * r) -
                                                              (3.52 * r * r * r));
            }
         }
         if (fine_dead_extinctions_density_b_[fuelModel].y > 0.0)
         {
            float r = deadMoisture / fine_dead_extinctions_density_b_[fuelModel].y;
            if (r < 1.0){
               reactionIntensity += areas_reaction_factors_b_[fuelModel].z * (1.0 -
                                                              (2.59 * r) +
                                                              (5.11 * r * r) -
                                                              (3.52 * r * r * r));
            }
         }

         float heatPerUnitArea = reactionIntensity * residence_flux_live_sav_b_[fuelModel].x;
         float baseSpreadRate = 0.0;

         if (heatOfIgnition > 0.0){
            baseSpreadRate = reactionIntensity * residence_flux_live_sav_b_[fuelModel].y / heatOfIgnition;
         }

         float slopeFactor = slope_wind_factors_b_[fuelModel].x * slope_aspect_elevation_t_[cell].x * slope_aspect_elevation_t_[cell].x;
         float windFactor = 0.0;
         if (xwind_[cell] > 0.0){
            windFactor = slope_wind_factors_b_[fuelModel].y * pow(xwind_[cell], slope_wind_factors_b_[fuelModel].z);
         }

         spreadModifier = slopeFactor + windFactor;

         float upslope;
         if (slope_aspect_elevation_t_[cell].y >= 180.0){
            upslope = slope_aspect_elevation_t_[cell].y - 180.0f;
         }
         else{
            upslope = slope_aspect_elevation_t_[cell].y + 180.0f;
         }

         int checkEffectiveWindspeed = 0;
         int updateEffectiveWindspeed = 0;
         float effectiveWindspeed = 0.0;
         if (baseSpreadRate <= 0.0)
         {
            maxSpreadRate = 0.0;
            spreadDirection = 0.0;
         }
         else if (spreadModifier <= 0.0)
         {
            maxSpreadRate = baseSpreadRate;
            spreadDirection = 0.0;
         }
         else if (slope_aspect_elevation_t_[cell].x <= 0.0)
         {
            effectiveWindspeed = xwind_[cell];
            maxSpreadRate = baseSpreadRate * (1.0 + spreadModifier);
            spreadDirection = ywind_[cell];
            checkEffectiveWindspeed = 1;
         }
         else if (xwind_[cell] <= 0.0)
         {
            maxSpreadRate = baseSpreadRate * (1.0 + spreadModifier);
            spreadDirection = upslope;
            updateEffectiveWindspeed = 1;
            checkEffectiveWindspeed = 1;
         }
         else if (fabs(ywind_[cell] - upslope) < 0.000001)
         {
            maxSpreadRate = baseSpreadRate * (1.0 + spreadModifier);
            spreadDirection = upslope;
            updateEffectiveWindspeed = 1;
            checkEffectiveWindspeed = 1;
         }
         else
         {
            float angleDelta;
            if (upslope <= ywind_[cell]) //https://play.google.com/store/books/details?id=rlOgae6p898C&rdid=book-rlOgae6p898C&rdot=1
               angleDelta = ywind_[cell] - upslope;
            else
               angleDelta = 360.0 - upslope + ywind_[cell];
            angleDelta *= 3.14159 / 180.0;
            float slopeRate = baseSpreadRate * slopeFactor;
            float windRate = baseSpreadRate * windFactor;
            float x = slopeRate + windRate * cos(angleDelta);
            float y = windRate * sin(angleDelta);
            float addedSpeed = sqrt(x * x + y * y);
            maxSpreadRate = baseSpreadRate + addedSpeed;

            spreadModifier = maxSpreadRate / baseSpreadRate - 1.0;
            if (spreadModifier > 0.0)
               updateEffectiveWindspeed = 1;
            checkEffectiveWindspeed = 1;

            float addedAngle = 0.0;
            if (addedSpeed > 0.0)
               addedAngle = asin(Clamp(fabs(y) / addedSpeed, -1.0, 1.0));
            float angleOffset = 0.0;
            if (x >= 0.0)
            {
               if (y >= 0.0)
                  angleOffset = addedAngle;
               else
                  angleOffset = 2.0 * 3.14159 - addedAngle;
            }
            else
            {
               if (y >= 0.0)
                  angleOffset = 3.14159 + addedAngle;
               else
                  angleOffset = 3.14159 - angleOffset;
            }
            spreadDirection = upslope + angleOffset * 180.0 / 3.14159;
            if (spreadDirection > 360.0)
               spreadDirection -= 360.0;
         }

         if (updateEffectiveWindspeed == 1)
         {
            effectiveWindspeed = pow((spreadModifier * slope_wind_factors_b_[fuelModel].w), (1.0 / slope_wind_factors_b_[fuelModel].z));
         }
         if (checkEffectiveWindspeed == 1)
         {
            float maxWind = 0.9 * reactionIntensity;
            if (effectiveWindspeed > maxWind)
            {
               if (maxWind <= 0.0)
                  spreadModifier = 0.0;
               else
                  spreadModifier = slope_wind_factors_b_[fuelModel].y * pow(maxWind, slope_wind_factors_b_[fuelModel].z);
               maxSpreadRate = baseSpreadRate * (1.0 + spreadModifier);
               effectiveWindspeed = maxWind;
            }
         }
         ellipseEccentricity = 0.0;
         if (effectiveWindspeed > 0.0)
         {
            float lengthWidthRatio = 1.0 + 0.002840909 * effectiveWindspeed;
            ellipseEccentricity = sqrt(lengthWidthRatio * lengthWidthRatio - 1.0) / lengthWidthRatio;
         }

         roth_data_[cell].w =
               3.4613 * (384.0 * (reactionIntensity / 0.189275)) *
               (0.30480060960) / (60.0 * fuel_sav_accel_b_[fuelModel].x);

         roth_data_[cell].x = maxSpreadRate;
         roth_data_[cell].y = spreadDirection;
         roth_data_[cell].z = ellipseEccentricity;
      }
   }
}

/*
Function: To ensure values aren't wildly out of proportion
Input: a value
Purpose: keeps the value in a range
*/
float FireSim::Clamp(float val, float flr, float ceiling){
   if(val >= flr && val <= ceiling){
      return val;
   }
   if(val < flr){
      return flr;
   }
   return ceiling;
}

/*
Function: BurnDistance
Input: distance, rate, timestep
Purpose: reduce distance for burning over several timesteps
*/
bool FireSim::BurnDistance(float &dist, float rate, float step){
   bool torched = false;
   // lower distance based on roth rate
   // t = d / r;
   // d = d - r * time_step_
   dist = dist - rate * step;
   if( dist < 0){
      dist *= -1;
//      dist = 0;
      torched = true;
   }
   return torched;
}

// Not sure why this isn't working, I have minimal knowledge of GDAL
float* GISToFloatArray(const char* fname, int interpWidth, int interpHeight, float &cell_size)
{
   // Important note ------ Gdal considers images to be north up
   // the origin of datasets takes place in the upper-left or North-West corner.
   // Now to create a GDAL dataset
   // auto ds = ((GDALDataset*) GDALOpen(fname,GA_ReadOnly));
   GDALDataset* ds = ((GDALDataset*) GDALOpen(fname,GA_ReadOnly));
   if(ds == NULL)
   {
      return NULL;
   }

   // Creating a Raster band variable
   // A band represents one whole dataset within a dataset
   // in your case your files have one band.
   GDALRasterBand  *poBand;
   int             nBlockXSize, nBlockYSize;
   int             bGotMin, bGotMax;
   double          adfMinMax[2];

   // Assign the band
   poBand = ds->GetRasterBand( 1 );
   poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );

   // find the min and max
   adfMinMax[0] = poBand->GetMinimum( &bGotMin );
   adfMinMax[1] = poBand->GetMaximum( &bGotMax );
   if( ! (bGotMin && bGotMax) )
      GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
   int min = adfMinMax[0];
   int max = adfMinMax[1];

   // get the width and height of the band or dataset
   int width = poBand->GetXSize();
   int height = poBand->GetYSize();

   // GDAL can handle files that have multiple datasets jammed witin it
   int bands = ds->GetRasterCount();

   // the float variable to hold the DEM!
   float *pafScanline;
   // std::std::cout << "Min: " << adfMinMax[0] << " Max: " << adfMinMax[1] << std::endl;
   int dsize = 256;
   // pafScanline = (T *) CPLMalloc(sizeof(T)*width*height);
   pafScanline = (float*) CPLMalloc(sizeof(float)*interpWidth*interpHeight);

   // Lets acquire the data.  ..... this funciton will interpolate for you
   // poBand->RasterIO(GF_Read,0,0,width,height,pafScanline,width,height,GDT_Float32,0,0);
   poBand->RasterIO(GF_Read,0,0,width,height,pafScanline,interpWidth,interpHeight,GDT_Int32,0,0);
   //        chage these two to interpolate automatically ^      ^

   // The Geotransform gives information on where a dataset is located in the world
   // and the resolution.
   // for more information look at http://www.gdal.org/gdal_datamodel.html
   double geot[6];
   ds->GetGeoTransform(geot);

   // Get the x resolution per pixel(south and west) and y resolution per pixel (north and south)
   float xres = geot[1];
   float yres = geot[5];
   
   // Calc Cell Size
   cell_size = (xres * width) / interpWidth;
   return pafScanline;

}

// This is working, so I used this to read the files and then cast as floats divided by 100 if needed
int* GISToIntArray(const char* fname, int interpWidth, int interpHeight, float &cell_size)
{
   // Important note ------ Gdal considers images to be north up
   // the origin of datasets takes place in the upper-left or North-West corner.
   // Now to create a GDAL dataset
   // auto ds = ((GDALDataset*) GDALOpen(fname,GA_ReadOnly));
   GDALDataset* ds = ((GDALDataset*) GDALOpen(fname,GA_ReadOnly));
   if(ds == NULL)
   {
      return NULL;
   }

   // Creating a Raster band variable
   // A band represents one whole dataset within a dataset
   // in your case your files have one band.
   GDALRasterBand  *poBand;
   int             nBlockXSize, nBlockYSize;
   int             bGotMin, bGotMax;
   double          adfMinMax[2];

   // Assign the band
   poBand = ds->GetRasterBand( 1 );
   poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );

   // find the min and max
   adfMinMax[0] = poBand->GetMinimum( &bGotMin );
   adfMinMax[1] = poBand->GetMaximum( &bGotMax );
   if( ! (bGotMin && bGotMax) )
      GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
   int min = adfMinMax[0];
   int max = adfMinMax[1];

   // get the width and height of the band or dataset
   int width = poBand->GetXSize();
   int height = poBand->GetYSize();

   // GDAL can handle files that have multiple datasets jammed witin it
   int bands = ds->GetRasterCount();

   // the float variable to hold the DEM!
   int *pafScanline;
   // std::std::cout << "Min: " << adfMinMax[0] << " Max: " << adfMinMax[1] << std::endl;
   int dsize = 256;
   // pafScanline = (T *) CPLMalloc(sizeof(T)*width*height);
   pafScanline = (int *) CPLMalloc(sizeof(int)*interpWidth*interpHeight);

   // Lets acquire the data.  ..... this funciton will interpolate for you
   // poBand->RasterIO(GF_Read,0,0,width,height,pafScanline,width,height,GDT_Float32,0,0);
   poBand->RasterIO(GF_Read,0,0,width,height,pafScanline,interpWidth,interpHeight,GDT_Int32,0,0);
   //        chage these two to interpolate automatically ^      ^

   // The Geotransform gives information on where a dataset is located in the world
   // and the resolution.
   // for more information look at http://www.gdal.org/gdal_datamodel.html
   double geot[6];
   ds->GetGeoTransform(geot);

   // Get the x resolution per pixel(south and west) and y resolution per pixel (north and south)
   float xres = geot[1];
   float yres = geot[5];
   // Calc Cell Size
   cell_size = (xres * width) / interpWidth;

   return pafScanline;

}

/////////////////////////////////////////////////////////////////////////////////////
// DEBUGGING TOOLS DEBUGGING TOOLS DEBUGGING TOOLS DEBUGGING TOOLS DEBUGGING TOOLS //
/////////////////////////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////////////////////////////////////
// // TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
// //////////////////////////////////////////////////////////////////////////////////////
// for(int row = 0; row < sim_dim_y_; row++){
//  for(int col = 0; col < sim_dim_x_; col++, cell++){
//   if(std::isnan(ywind_[cell]) || std::isinf(ywind_[cell])){
//    std::cout<<"We messed up in FireSim.cpp ywind_: " << ywind_[cell] << std::endl;
//   }
//  }
// }
// cell = 0;
// //////////////////////////////////////////////////////////////////////////////////////
// // TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
// //////////////////////////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////////////////////////////////////
// // TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
// //////////////////////////////////////////////////////////////////////////////////////
// for(int row = 0; row < sim_dim_y_; row++){
//  for(int col = 0; col < sim_dim_x_; col++){
//   for(int dirs = 0; dirs < 16; dirs++, cell++)
//   if(std::isnan(csr_[cell]) || std::isinf(csr_[cell])){
//    std::cout<<"We messed up in FireSim.cpp csr_: " << csr_[cell] << std::endl;
//   }
//  }
// }
// cell = 0;
// //////////////////////////////////////////////////////////////////////////////////////
// // TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
// //////////////////////////////////////////////////////////////////////////////////////    
