#include <curand_mtgp32_kernel.h>
#include "kernel_common.h"

// http://learningcenter.firewise.org/Firefighter-Safety/index.php
// https://www.slideshare.net/WildlandFire/s290-ep-10-45489519?qid=2d19fd09-0dd3-4ccf-9152-747fc02febd9&v=&b=&from_search=7
// https://www.bugwood.org/pfire/weather.html

/////////////////////////////////////////////////////////////////////////////
//                             InitialSetup
/////////////////////////////////////////////////////////////////////////////
__global__ void InitialSetup(int* g_fuel_t_, float4* g_roth_data_, float4* g_dead_sav_burnable_b_, float4* g_dead_1h_b_, float4* g_dead_10h_b_,
                             float4* g_dead_100h_b_, float4* g_live_h_b_, float4* g_live_w_b_, float4* g_fine_dead_extinctions_density_b_,
                             float4* g_areas_reaction_factors_b_, float4* g_slope_wind_factors_b_, float4* g_residence_flux_live_sav_b_,
                             float2* g_fuel_sav_accel_b_, float3* g_slope_aspect_elevation_t_, float3* g_dead_moistures_t_,
                             float2* g_live_moistures_t_, int* g_x_wind, int* g_y_wind, int size, int* g_DynLiveH, int* g_DynDead1)
{
   int dirs = 16;
   int cell = (blockIdx.x * blockDim.x + threadIdx.x) * dirs;
   int cell_im = (blockIdx.x * blockDim.x + threadIdx.x);
   dynamicMositure = false;
   precipitation = false;

   while(cell < size*dirs)
   {
      int fuelModel = g_fuel_t_[cell_im];

      if(fuelModel == 99)
      {
         g_roth_data_[cell_im].w = 1.0f;
         g_roth_data_[cell_im].x = 0.0f;
         g_roth_data_[cell_im].y = 1.0f;
         g_roth_data_[cell_im].z = 1.0f;
      }
      else
      {
         if(g_dead_sav_burnable_b_[fuelModel].w < 50.0)
         {
            printf("Warning: Something may have gone wrong. Check that Files were read Correctly\n");
         }

         float maxSpreadRate = 0.f;
         float ellipseEccentricity = 0.f;
         float spreadDirection = 0.f;
         float spreadModifier = 0.f;
         float dynamicLiveHerbaceousMoisture;
         float dynamicDead1HrValue;
         float dLH, dD1, x, y;
         float3 timeLagClass;

         if(dynamicMositure)
         {
            if(precipitation)
            {
               dLH = (float) g_DynLiveH[cell_im] / 100;
               dD1 = (float) g_DynDead1[cell_im] / 100;

               if(dLH != g_live_moistures_t_[fuelModel].x)
               {
                  dynamicLiveHerbaceousMoisture = dLH;
               }
               else
               {
                  dynamicLiveHerbaceousMoisture = g_live_moistures_t_[fuelModel].x;
               }

               if(dD1 != g_dead_moistures_t_[fuelModel].x)
               {
                  dynamicDead1HrValue = dD1;
               }
               else
               {
                  dynamicDead1HrValue = g_dead_moistures_t_[fuelModel].x;
               }
            }
            else
            {
               x = dynamicDead1HrValue = g_dead_moistures_t_[fuelModel].x;
               y = dynamicLiveHerbaceousMoisture = g_live_moistures_t_[fuelModel].x;

               if(y > .30 && y < 1.20)
               {
                  dynamicLiveHerbaceousMoisture = (((-1) * x) + 1.20);
                  dynamicDead1HrValue = (-10/9) * (y-1.20);
               }
               if(y <= .30)
               {
                  dynamicLiveHerbaceousMoisture = .0f;
                  dynamicDead1HrValue = 1.0f;    
               }
               if(y >= 1.20)
               {
                  dynamicLiveHerbaceousMoisture = g_live_moistures_t_[fuelModel].x;
                  dynamicDead1HrValue = 0.0f;
               }
               g_dead_moistures_t_[fuelModel].x = dynamicDead1HrValue;
               g_live_moistures_t_[fuelModel].x = dynamicLiveHerbaceousMoisture;
            }
         }
         else
         {
            dynamicLiveHerbaceousMoisture = g_live_moistures_t_[fuelModel].x;
            dynamicDead1HrValue = g_dead_moistures_t_[fuelModel].x;
         }
         if (g_dead_sav_burnable_b_[fuelModel].x > 192.0)
         {
            timeLagClass.x = dynamicDead1HrValue;
         }
         else if (g_dead_sav_burnable_b_[fuelModel].x > 48.0)
         {
            timeLagClass.x = g_dead_moistures_t_[fuelModel].y;
         }
         else
         {
            timeLagClass.x = g_dead_moistures_t_[fuelModel].z;
         }

         if (g_dead_sav_burnable_b_[fuelModel].y > 192.0)
         {
            timeLagClass.y = dynamicDead1HrValue;
         }
         else if (g_dead_sav_burnable_b_[fuelModel].y > 48.0)
         {
            timeLagClass.y = g_dead_moistures_t_[fuelModel].y;
         }
         else
         {
            timeLagClass.y = g_dead_moistures_t_[fuelModel].z;
         }

         if (g_dead_sav_burnable_b_[fuelModel].z > 192.0)
         {
            timeLagClass.z = dynamicDead1HrValue;
         }
         else if (g_dead_sav_burnable_b_[fuelModel].z > 48.0)
         {
            timeLagClass.z = g_dead_moistures_t_[fuelModel].y;
         }
         else
         {
            timeLagClass.z = g_dead_moistures_t_[fuelModel].z;
         }

         float weightedFuelModel =
               timeLagClass.x * g_dead_1h_b_[fuelModel].x * g_dead_1h_b_[fuelModel].y +
               timeLagClass.y * g_dead_10h_b_[fuelModel].x * g_dead_10h_b_[fuelModel].y +
               timeLagClass.z * g_dead_100h_b_[fuelModel].x * g_dead_100h_b_[fuelModel].y;

         float fuelMoistures[5];
         fuelMoistures[0] = timeLagClass.x;
         fuelMoistures[1] = timeLagClass.y;
         fuelMoistures[2] = timeLagClass.z;
         fuelMoistures[3] = dynamicLiveHerbaceousMoisture;
         fuelMoistures[4] = g_live_moistures_t_[fuelModel].y;

         float liveExtinction = 0.0;
         if(g_live_h_b_[fuelModel].y > 0.0 || g_live_w_b_[fuelModel].y > 0.0)
         {
            float fineDeadMoisture = 0.0;
            if (g_fine_dead_extinctions_density_b_[fuelModel].x > 0.0)
            {
               fineDeadMoisture = weightedFuelModel / g_fine_dead_extinctions_density_b_[fuelModel].x;
            }

            liveExtinction = (g_fine_dead_extinctions_density_b_[fuelModel].z * 
                             (1.0f - fineDeadMoisture / g_fine_dead_extinctions_density_b_[fuelModel].y)) - 0.226f;
            liveExtinction = fmax(liveExtinction, g_fine_dead_extinctions_density_b_[fuelModel].y);
         }

         float heatOfIgnition =
               g_areas_reaction_factors_b_[fuelModel].x *
               ((250.0f + 1116.0f * fuelMoistures[0]) * g_dead_1h_b_[fuelModel].z * g_dead_1h_b_[fuelModel].x +
                (250.0f + 1116.0f * fuelMoistures[1]) * g_dead_10h_b_[fuelModel].z * g_dead_10h_b_[fuelModel].x +
                (250.0f + 1116.0f * fuelMoistures[2]) * g_dead_100h_b_[fuelModel].z * g_dead_100h_b_[fuelModel].x) +
               g_areas_reaction_factors_b_[fuelModel].y *
               ((250.0f + 1116.0f * fuelMoistures[3]) * g_live_h_b_[fuelModel].z * g_live_h_b_[fuelModel].x +
                (250.0f + 1116.0f * fuelMoistures[4]) * g_live_w_b_[fuelModel].z * g_live_w_b_[fuelModel].x);

         heatOfIgnition *= g_fine_dead_extinctions_density_b_[fuelModel].w;

         float liveMoisture = g_live_h_b_[fuelModel].z * fuelMoistures[3] + g_live_w_b_[fuelModel].z * fuelMoistures[4];

         float deadMoisture = g_dead_1h_b_[fuelModel].z * fuelMoistures[0] +
                              g_dead_10h_b_[fuelModel].z * fuelMoistures[1] +
                              g_dead_100h_b_[fuelModel].z * fuelMoistures[2];

         float reactionIntensity = 0.0;

         if (liveExtinction > 0.0)
         {
            float r = liveMoisture / liveExtinction;
            if (r < 1.0)
            {
               reactionIntensity += g_areas_reaction_factors_b_[fuelModel].w * (1.0 -
                                                              (2.59 * r) +
                                                              (5.11 * r * r) -
                                                              (3.52 * r * r * r));
            }
         }
         if (g_fine_dead_extinctions_density_b_[fuelModel].y > 0.0)
         {
            float r = deadMoisture / g_fine_dead_extinctions_density_b_[fuelModel].y;
            if (r < 1.0)
            {
               reactionIntensity += g_areas_reaction_factors_b_[fuelModel].z * (1.0 -
                                                              (2.59 * r) +
                                                              (5.11 * r * r) -
                                                              (3.52 * r * r * r));
            }
         }

         // float heatPerUnitArea = reactionIntensity * g_residence_flux_live_sav_b_[fuelModel].x;

         float baseSpreadRate = 0.0;

         if (heatOfIgnition > 0.0)
         {
            baseSpreadRate = reactionIntensity * g_residence_flux_live_sav_b_[fuelModel].y / heatOfIgnition;
         }
         float slopeFactor = g_slope_wind_factors_b_[fuelModel].x * g_slope_aspect_elevation_t_[cell_im].x * g_slope_aspect_elevation_t_[cell_im].x;
         float windFactor = 0.0;
         float X_wind = (float) g_x_wind[cell_im];
         float Y_wind = (float) g_y_wind[cell_im];
         if (X_wind > 0)
         {
            windFactor = g_slope_wind_factors_b_[fuelModel].y * powf(X_wind, g_slope_wind_factors_b_[fuelModel].z);
         }

         spreadModifier = slopeFactor + windFactor;

         float upslope;
         if (g_slope_aspect_elevation_t_[cell_im].y >= 180.0)
         {
            upslope = g_slope_aspect_elevation_t_[cell_im].y - 180.0f;
         }
         else
         {
            upslope = g_slope_aspect_elevation_t_[cell_im].y + 180.0f;
         }

         int checkEffectiveWindspeed = 0;
         int updateEffectiveWindspeed = 0;
         float effectiveWindspeed = 0.0;
         if (baseSpreadRate <= 0.0)
         {
            maxSpreadRate = 0.0;
            spreadDirection = 0.0;
         }
         else if (spreadModifier <= 0)
         {
            maxSpreadRate = baseSpreadRate;
            spreadDirection = 0.0;
         }
         else if (g_slope_aspect_elevation_t_[cell_im].x < 0)
         {
            effectiveWindspeed = X_wind;
            maxSpreadRate = baseSpreadRate * (1.0 + spreadModifier);
            spreadDirection = Y_wind;
            checkEffectiveWindspeed = 1;
         }
         else if (X_wind <= 0)
         {
            maxSpreadRate = baseSpreadRate * (1.0 + spreadModifier);
            spreadDirection = upslope;
            updateEffectiveWindspeed = 1;
            checkEffectiveWindspeed = 1;
         }
         else if (fabsf(Y_wind - upslope) < 0.000001)
         {
            maxSpreadRate = baseSpreadRate * (1.0 + spreadModifier);
            spreadDirection = upslope;
            updateEffectiveWindspeed = 1;
            checkEffectiveWindspeed = 1;
         }
         else
         {
            float angleDelta;
            if (upslope <= Y_wind)
            {
               angleDelta = Y_wind - upslope;
            }
            else
            {
               angleDelta = 360.0 - upslope + Y_wind;
            }
            angleDelta *= 3.14159 / 180.0;
            float slopeRate = baseSpreadRate * slopeFactor;
            float windRate = baseSpreadRate * windFactor;
            float x = slopeRate + windRate * cosf(angleDelta);
            float y = windRate * sinf(angleDelta);
            float addedSpeed = sqrtf(x * x + y * y);
            maxSpreadRate = baseSpreadRate + addedSpeed;
            spreadModifier = maxSpreadRate / baseSpreadRate - 1.0;
            if (spreadModifier > 0.0)
            {
               updateEffectiveWindspeed = 1;
            }
            checkEffectiveWindspeed = 1;

            float addedAngle = 0.0;
            if (addedSpeed > 0.0)
            {
               float val = fabsf(y) / addedSpeed;
               addedAngle = asinf(Clamp(val, -1.0f, 1.0f));
            }
            float angleOffset = 0.0;
            if (x >= 0.0)
            {
               if (y >= 0.0)
               {
                  angleOffset = addedAngle;
               }
               else
               {
                  angleOffset = 2.0 * 3.14159 - addedAngle;
               }
            }
            else
            {
               if (y >= 0.0)
               {
                  angleOffset = 3.14159 + addedAngle;
               }
               else
               {
                  angleOffset = 3.14159 - angleOffset;
               }
            }
            spreadDirection = upslope + angleOffset * 180.0 / 3.14159;
            if (spreadDirection > 360.0)
            {
               spreadDirection -= 360.0;
            }
         }
         if (updateEffectiveWindspeed == 1)
         {
            effectiveWindspeed = powf((spreadModifier * g_slope_wind_factors_b_[fuelModel].w), (1.0 / g_slope_wind_factors_b_[fuelModel].z));
         }
         if (checkEffectiveWindspeed == 1)
         {
            float maxWind = 0.9 * reactionIntensity;
            if (effectiveWindspeed > maxWind)
            {
               if (maxWind <= 0.0)
               {
                  spreadModifier = 0.0;
               }
               else
               {
                  spreadModifier = g_slope_wind_factors_b_[fuelModel].y * powf(maxWind, g_slope_wind_factors_b_[fuelModel].z);
               }
               maxSpreadRate = baseSpreadRate * (1.0 + spreadModifier);
               effectiveWindspeed = maxWind;
            }
         }
         ellipseEccentricity = 0.0;
         if (effectiveWindspeed > 0.0)
         {
            float lengthWidthRatio = 1.0 + 0.002840909 * effectiveWindspeed;
            ellipseEccentricity = sqrtf(lengthWidthRatio * lengthWidthRatio - 1.0) / lengthWidthRatio;
         }
         g_roth_data_[cell_im].w = 3.4613 * (384.0 * (reactionIntensity / 0.189275)) * 
                                   (0.30480060960) / (60.0 * g_fuel_sav_accel_b_[fuelModel].x);
         g_roth_data_[cell_im].x = maxSpreadRate;
         g_roth_data_[cell_im].y = spreadDirection;
         g_roth_data_[cell_im].z = ellipseEccentricity;
      }
      cell += (blockDim.x * gridDim.x)*dirs;
      cell_im += (blockDim.x * gridDim.x);
   }
}

/////////////////////////////////////////////////////////////////////////////
//                             UpdateSpread
/////////////////////////////////////////////////////////////////////////////
__global__ void UpdateSpread(int* g_fuel_t_, float4* g_roth_data_, float2* g_fuel_sav_accel_b_, float* g_maxspreadrate_, 
                             float* g_acceleration_constant_, float* g_intensity_modifier_, float* g_angles_, int size){
   int dirs = 16;
   int cell = (blockIdx.x * blockDim.x + threadIdx.x) * dirs;
   int cell_im = (blockIdx.x * blockDim.x + threadIdx.x);
   int fuel_model;

   while(cell < size*dirs)
   {
      for (unsigned int i = 0; i < dirs; i++)
      {
            g_maxspreadrate_[cell+i] = (float) (g_roth_data_[cell_im].x * (1.0f - g_roth_data_[cell_im].z) / 
                                            (1.0f - g_roth_data_[cell_im].z * cos(g_roth_data_[cell_im].y * 3.14159f / 180.f - g_angles_[cell_im])));
      }

      fuel_model = g_fuel_t_[cell_im];
      g_intensity_modifier_[cell_im] = g_roth_data_[cell_im].w;
      g_acceleration_constant_[cell_im] = g_fuel_sav_accel_b_[fuel_model].y;

      cell += (blockDim.x * gridDim.x)*dirs;
      cell_im += (blockDim.x * gridDim.x);
   }

}

/////////////////////////////////////////////////////////////////////////////
//                             Accelerate
/////////////////////////////////////////////////////////////////////////////
__global__ void Accelerate(float* curspreadrate, float* maxspreadrate, int size, float time_step){
   int dirs = 16;
   int cell = (blockIdx.x * blockDim.x + threadIdx.x) * dirs;
   float current, ratio;

   while(cell < size*dirs)
   {
      if(maxspreadrate[cell] <= 0.0f){ // BECAUSE IF THIS IS LESS THAN 0 THEN THE FIRE HAS NOT REACHED THESE CELLS
         cell += (blockDim.x * gridDim.x)*dirs; // SKIP THEM FOR NOW
         continue;
      }
      for (unsigned int i = 0; i < dirs; i++, cell++) // ACCELERATE EVERY CELL AND ITS NIEGHBORS
      {
         current = !(maxspreadrate[cell+i]<curspreadrate[cell+i]) ? curspreadrate[cell+i] : maxspreadrate[cell];
         ratio = current / maxspreadrate[cell];
         curspreadrate[cell] = Clamp(time_step - ratio, 0.0f,1.0f) * (maxspreadrate[cell] - current) + current;
      }
      cell += (blockDim.x * gridDim.x)*dirs;
   }
}

/////////////////////////////////////////////////////////////////////////////
//                   Test Crown Rate With Spotting Enabled
/////////////////////////////////////////////////////////////////////////////
__global__ void TestCrownRate_WithSpot(float* curspreadrate, float* maxspreadrate, float* intensity_modifier, int size,
                                       float* I_o, float* RAC, bool* emberMap, float* emberTestMap, float* canopyHeight,
                                       Ember* ember_list, int* g_x_wind, int* g_y_wind, int* g_fuel_t_){
   // Cell Id's and Locations
   int dirs = 16;
   int cell = (blockIdx.x * blockDim.x + threadIdx.x) * dirs;
   int cell_im = (blockIdx.x * blockDim.x + threadIdx.x);

   // Crowning Values
   float I_b, R_max_crown, surface_fuel_consumption, crown_coeff, CFB, crown_rate;

   // Spotting Values
   float t_t, t_0, t_1, t_2, t_3, z_F, z_o, velocity_ratio, r, maxEmberHeight;

   // Spotting Constants
   float a_x = 5.963;
   float b_x = 4.563;
   float D_p = 0.01;
   int   B = 40;
   float X_wind,Y_wind;

   while(cell < size*dirs)
   {
      if(maxspreadrate[cell] <= 0.f)
      {
         cell += (blockDim.x * gridDim.x)*dirs;
         cell_im += (blockDim.x * gridDim.x);
         continue;
      }
      for (unsigned int i = 0; i < dirs; i++)
      {
         I_b = curspreadrate[cell+i] * intensity_modifier[cell_im];
         if(I_b > I_o[cell_im])
         {
            if(!emberMap[cell+i])
            {
               z_F = (float) (0.0775 * powf(I_b,0.46) + canopyHeight[cell_im]);
               z_o = 0.4306;
               velocity_ratio = (float) (B * powf(D_p / z_F, 0.5));
               r = (float)powf((b_x+z_o/z_F)/a_x, 0.5);
               t_0 = 0.7f;
               t_1 = (float) ( 1.f - powf((z_o / z_F), 0.5) +
                     velocity_ratio*log((1.f - velocity_ratio) / (powf((z_o / z_F), 0.5) - velocity_ratio)));
               t_2 = (float) (0.2f + B*powf((D_p / z_F), 0.5) *
                                        ( 1.0f + B*powf((D_p / z_F), 0.5) *
                                        logf(1.f + 1.f/(1.f - powf((D_p / z_F), 0.5))) ));
               t_3 = (float) (a_x / (0.8*velocity_ratio) * ( logf((1.f - 0.8*velocity_ratio)/(1.f - 0.8 * r * velocity_ratio))
                                                                - 0.8*velocity_ratio*(r - 1) -0.5f * powf(.8*velocity_ratio, 2) * powf(r-1, 2) ));
               t_t = t_0 + t_1 + t_2 + t_3;
               maxEmberHeight = (b_x*z_F + z_F * a_x * powf((3 * (t_t - t_0 - 1.2) )/ a_x + 1.f, (2./3.)));
               ember_list[cell+i].x_ = cell+i;
               ember_list[cell+i].y_ = cell_im;
               ember_list[cell+i].z_ = maxEmberHeight;
               ember_list[cell+i].z_o_ = maxEmberHeight;
               X_wind = (float) g_x_wind[cell_im];
               Y_wind = (float) g_y_wind[cell_im];
               ember_list[cell+i].angle_ = atan2f( (double) Y_wind, (double) X_wind);
               ember_list[cell+i].magnitude_ = (powf((powf(X_wind,2) + powf(Y_wind,2)),0.5));
               emberMap[cell+i] = true;
               emberTestMap[cell+i] = t_t;
            }
            
            R_max_crown = 3.34f * maxspreadrate[cell+i];
            R_max_crown = Clamp(R_max_crown, 0, 10); // This is used to control a realistic crown fire rate maximum
            surface_fuel_consumption = I_o[cell_im] * curspreadrate[cell+i] / I_b;
            crown_coeff = (float) (-logf(0.1) / (0.9 * (RAC[cell_im] - surface_fuel_consumption)));
            CFB = (float) (1.0 - expf(-1*crown_coeff * (curspreadrate[cell+i] - surface_fuel_consumption)));
            CFB = Clamp(CFB, 0.0, 1.0);
            crown_rate = curspreadrate[cell] + CFB * (R_max_crown - curspreadrate[cell]);
            if(crown_rate >= RAC[cell_im]){
               maxspreadrate[cell+i] = (crown_rate > maxspreadrate[cell+i] ? crown_rate : maxspreadrate[cell+i]);
            }
         }
      }
      cell += (blockDim.x * gridDim.x)*dirs;
      cell_im += (blockDim.x * gridDim.x);
   }
}

/////////////////////////////////////////////////////////////////////////////
//                   Test Crown Rate With Spotting Disabled
/////////////////////////////////////////////////////////////////////////////
__global__ void TestCrownRate_NoSpot(float* curspreadrate, float* maxspreadrate, float* intensity_modifier,
                                     int size, float* I_o, float* RAC){
   // Cell Id's and Locations
   int dirs = 16;
   int cell = (blockIdx.x * blockDim.x + threadIdx.x) * dirs;
   int cell_im = (blockIdx.x * blockDim.x + threadIdx.x);

   // Crowning Values
   float I_b, R_max_crown, surface_fuel_consumption, crown_coeff, CFB, crown_rate;

   while(cell < size*dirs)
   {
      if(maxspreadrate[cell] <= 0.f)
      {
         cell += (blockDim.x * gridDim.x)*dirs;
         cell_im += (blockDim.x * gridDim.x);
         continue;
      }
      for (unsigned int i = 0; i < dirs; i++)
      {
         I_b = curspreadrate[cell+i] * intensity_modifier[cell_im];
         if(I_b > I_o[cell_im])
         {  
            R_max_crown = 3.34f * maxspreadrate[cell+i];
            R_max_crown = Clamp(R_max_crown, 0, 10); // This is used to control a realistic crown fire rate maximum
            surface_fuel_consumption = I_o[cell_im] * curspreadrate[cell+i] / I_b;
            crown_coeff = (float) (-logf(0.1) / (0.9 * (RAC[cell_im] - surface_fuel_consumption)));
            CFB = (float) (1.0 - expf(-1*crown_coeff * (curspreadrate[cell+i] - surface_fuel_consumption)));
            CFB = Clamp(CFB, 0.0, 1.0);
            crown_rate = curspreadrate[cell] + CFB * (R_max_crown - curspreadrate[cell]);
            if(crown_rate >= RAC[cell_im])
            {
               maxspreadrate[cell+i] = (crown_rate > maxspreadrate[cell+i] ? crown_rate : maxspreadrate[cell+i]);
            }
         }
      }
      cell += (blockDim.x * gridDim.x)*dirs;
      cell_im += (blockDim.x * gridDim.x);
   }
}

/////////////////////////////////////////////////////////////////////////////
//                             Test Spotting
/////////////////////////////////////////////////////////////////////////////
__global__ void TestSpotting(bool* emberMap, float* canopyHeight, Ember* ember_list, float currentTime, float time_step,
                             int size, int* g_fuel_t_, int* g_DynDead1, float3* g_dead_moistures_t_){
   int dirs = 16;
   int cell = (blockIdx.x * blockDim.x + threadIdx.x) * dirs;
   int cell_im = (blockIdx.x * blockDim.x + threadIdx.x);
   float p_a = 0.0012;
   float p_s = 0.3f;
   float g = 9.8f;
   float D_p = 0.01;
   float K = 0.0064;
   float C_d = 1.2f;
   float v_o = (float) powf((M_PI*g*p_s*D_p)/(2*C_d*p_a),0.5f);
   float tau = (float) ((4*C_d*v_o)/(K*M_PI*g));
   double z_o, dX;
   z_o = 0.4306;
   float dD1, chance;
   curandState_t* state;

   while(cell < size*dirs)
   {
      if(emberMap[cell] == true)
      {
         if(ember_list[cell].z_ <= canopyHeight[cell_im])
         {
            if(g_fuel_t_[cell_im] == 99)
            {
               for (unsigned int i = 0; i < dirs; i++)
               {
                  ember_list[cell+i].ember_toa = 32767;
                  emberMap[cell+i] = false;
               }
            }
            else
            {
               // http://www.fbfrg.org/ --> http://www.fbfrg.org/fuel-moisture/1-hour-moisture-content
               // Probability of Ember starting a fire
               if(precipitation)
               {
                  dD1 = g_DynDead1[cell_im];
                  if(dD1 > 17) // no chance of burning
                  {
                     for (unsigned int i = 0; i < dirs; i++)
                     {
                        ember_list[cell+i].ember_toa = 32767;
                        emberMap[cell+i] = false;
                     }
                  }
                  if(dD1 >= 15 && dD1 <= 17) // 10% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell])); // randomize between 0.0 and 1.0 to remove mod function
                     if(chance < 0.10)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 12 && dD1 <= 20) // 20% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.20)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 10 && dD1 <= 11) // 30% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.30)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 9 && dD1 <= 8) // 40% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.40)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 7 && dD1 <= 6) // 50% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.50)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 4 && dD1 <= 5) // 70% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.70)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 2 && dD1 <= 3) // 90% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.90)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
               }
               else
               {
                  dD1 = g_dead_moistures_t_[g_fuel_t_[cell_im]].x;
                  if(dD1 > 17) // no chance of burning
                  {
                     for (unsigned int i = 0; i < dirs; i++, cell++)
                     {
                        ember_list[cell+i].ember_toa = 32767;
                        emberMap[cell+i] = false;
                     }
                  }
                  if(dD1 >= 15 && dD1 <= 17) // 10% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.10)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 12 && dD1 <= 14) // 20% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.20)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 10 && dD1 <= 11) // 30% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.30)
                     {
                        for (unsigned int i = 0; i < dirs; i++, cell++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 9 && dD1 <= 8) // 40% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.40)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 7 && dD1 <= 6) // 50% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.50)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 4 && dD1 <= 5) // 70% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.70)
                     {
                        for (unsigned int i = 0; i < dirs; i++, cell++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
                  if(dD1 >= 2 && dD1 <= 3) // 90% chance of catching fire
                  {
                     chance = curand_uniform(&(state[cell]));
                     if(chance < 0.90)
                     {
                        for (unsigned int i = 0; i < dirs; i++)
                        {
                           ember_list[cell+i].ember_toa = currentTime;
                           emberMap[cell+i] = false;
                        }
                     }
                  }
               }
            }
         }

         // move ember
         dX = ((ember_list[cell].magnitude_ * log(ember_list[cell].z_ / z_o)) / log(canopyHeight[cell_im] / z_o)); // rate of change
         dX = time_step * dX; // magnitude of change -- will be different for all kernels :(
         ember_list[cell].x_ -= (int)(dX * cos(ember_list[cell].angle_)); // this will truncate the values to integers, but that's an acceptable error for now
         ember_list[cell].y_ -= (int)(dX * sin(ember_list[cell].angle_));
         ember_list[cell].z_ = (float) (ember_list[cell].z_ - v_o * ((time_step) / tau - 0.5f * powf((time_step) / tau, 2))); // update ember height
         if(ember_list[cell].z_ <= 0)
         {
            ember_list[cell].z_ = 0.0;            
         }
         if(ember_list[cell].x_ <= 0)
         {
            ember_list[cell].x_ = 0;         
         }
         if(ember_list[cell].y_ <= 0)
         {
            ember_list[cell].y_ = 0;
         }
      }  
      cell += (blockDim.x * gridDim.x)*dirs;
      cell_im += (blockDim.x * gridDim.x);
   }
}

/////////////////////////////////////////////////////////////////////////////
//                       Find Max Value in an Array printf("maxSpread[blockIdx.x]: %f\n", maxSpread[0]);
/////////////////////////////////////////////////////////////////////////////
__global__ void findMax(float* maxSpreadRates, float* maxSpread, float* cellSize, float* ts, int elements){

   extern __shared__ float sdata[];

   int tid = threadIdx.x, gid = (blockDim.x * gridDim.x) + tid;

   sdata[tid] = -1.0f;

   while(gid < elements){
      sdata[tid] = max(sdata[tid], maxSpreadRates[gid]);
      gid += gridDim.x * blockDim.x;
   }
   __syncthreads();

   gid = (blockDim.x * gridDim.x);
   // do reduction in shared mem
   for(unsigned int s = blockDim.x/2; s > 0; s >>= 1) {
      if (tid < s && gid < elements) {
         sdata[tid] = max(sdata[tid+s], sdata[tid]);
      }
      __syncthreads();
   }

   // write result for this block to global mem
   if (tid == 0){
      findMaxHelper(maxSpread, sdata[0]);
   }
   *ts = *maxSpread / *cellSize; // calculating the timestep here...
}

__device__ float findMaxHelper(float* address, float value){
 
   int* int_address = (int*) address;
   int old_value = * int_address, temp;
   while(value > __int_as_float(old_value)){
      temp = old_value;
      old_value = atomicCAS(int_address, temp, float_as_int(value));
   }
   return __int_as_float(old_value);
}

/////////////////////////////////////////////////////////////////////////////
//                               Clamp
/////////////////////////////////////////////////////////////////////////////
__device__ float Clamp(float val, float flr, float ceiling){
   if(val >= flr && val <= ceiling)
   {
      return val;
   }
   if(val < flr)
   {
      return flr;
   }
   return ceiling;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// DEBUGGING TOOLS DEBUGGING TOOLS DEBUGGING TOOLS DEBUGGING TOOLS DEBUGGING TOOLS DEBUGGING TOOLS
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

// if(isnan(curspreadrate[cell+i]) || isinf(curspreadrate[cell+i]) || curspreadrate[cell+i] > 500.0){
//    printf("A1--curspreadrate[cell+i]: %f\n;", curspreadrate[cell+i]);
//    asm("trap;");
// }

// if(isnan(maxspreadrate[cell+i]) || isinf(maxspreadrate[cell+i] )){
//    printf("A2--maxspreadrate[cell+i]: %f\n;", maxspreadrate[cell+i] || maxspreadrate[cell+i] > 500.0);
//    asm("trap;");
// }

// int current_val = atomicAdd(&spock, 1);   // put this somewhere in the code for testing purposes
// printf("spock: %i\n", spock);             // it is a very useful debugging tool
