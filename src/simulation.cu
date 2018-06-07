#include <iostream>
#include <fstream>
#include "FireSim.h"
#include <sstream>
#include <string>
#include <cstring>
#include "BD.h"
#include "MT.h"
#include "IMT.h"
#include "SEQ.h"
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#define GetDir getcwd

// ARGV #1
int _START_TICK_ = 0; // 0 is default state

// ARGV #2
int _STOP_TICK_ = 0; // 0 is a default state

// ARGV #3
int _FUEL_ = 0;      // 0 is a default state OFF

// ARGV #4
int _MOISTURE_ = 0;  // 0 is a default state OFF

// ARGV #5
int _WIND_XY_ = 0;   // 0 is a default state OFF

// ARGV #6
int _CROWNING_ = 0;  // 0 is a default state OFF

// ARGV #7
int _SPOTTING_ = 0;  // 0 is a default state OFF

// ARGV #8
int _PARALLEL_ = 0;  // 0 is sequential firesim


std::string double_to_string(double value)
{
  std::ostringstream ss;
  ss << value;
  return ss.str();
}

std::string int_to_string(int value)
{
  std::ostringstream ss;
  ss << value;
  return ss.str();
}

std::string current_data_dir_()
{
  char temp[FILENAME_MAX];
  std::string my__data_dir_ = GetDir(temp, FILENAME_MAX);
  return my__data_dir_;
}

int main(int argc, char *argv[])
{
  int step_size = 600;  // if the sim is bigger i.e. more data a higher value may be warranted
  double t_upSpread;    // used to time the simulation, for throughput and research validation
  int s = 512;          // # of blocks and threads can be changed for performance optimization

  _START_TICK_ = atoi(argv[1]); // At what time stamp/tick are we starting/resuming at
  _STOP_TICK_  = atoi(argv[2]);  // At What time stamp/tick are we stopping/pausing at
  _FUEL_       = atoi(argv[3]); // for fuel file in data
  _MOISTURE_   = atoi(argv[4]); // for moisture file in data
  _WIND_XY_    = atoi(argv[5]); // for windx and windy files, there's two of them
  _CROWNING_   = atoi(argv[6]); // to turn on or off crowning calculation
  _SPOTTING_   = atoi(argv[7]); // to turn on of off spotting, if crowning is off then this can't be on
  _PARALLEL_   = atoi(argv[8]); // to switch between parallel and sequential implementations

  if(_START_TICK_ == 0 && _STOP_TICK_ == 0)
  {
    std::cerr<<"You cannot run a simulation from 0 to 0... try it again eveyone makes mistakes";
    return 1;
  }
  if(_START_TICK_ < 0 || _STOP_TICK_ < 0)
  {
    std::cerr<<"You cannot run a simulation from a negative time tick or end on a negative... try it again dummy";
    return 1;
  }
  if(_START_TICK_ > _STOP_TICK_)
  {
    std::cerr<<"You cannot have a smaller start tick than stop tick... try it again back to the future wonder boy";
    return 1;
  }
  if(std::isnan(_START_TICK_) || std::isnan(_STOP_TICK_) || std::isinf(_START_TICK_) || std::isinf(_STOP_TICK_)) 
  {
    std::cerr<<"Uh integers please between 0 and soft cap at 600... try it again asshole";
    return 1;
  }

  // this grabs the current directory which is the build folder
  std::string _data_dir_ = current_data_dir_(); // this is for the files located in the data folder
  std::string _out_dir_; // this is for the files located in the out folder

  _out_dir_ = _data_dir_; // now data and out strings are the same follow below to see why

  int size = _data_dir_.size() - 1; // this is so we can resize out string

  _out_dir_.resize(size); // because the out directory is 1 element shorter than the data

  // pretty much replacing "/Desktop/spotting/buil" with "/Desktop/spotting/out/"
  int index = 0;
  while(index < _data_dir_.size() - 5)
  {
    index++;
  }
  _out_dir_[index] = 'o';
  _out_dir_[index+1] = 'u';
  _out_dir_[index+2] = 't';
  _out_dir_[index+3] = '/';

  // std::cout<<std::endl<<std::endl<<"checking the file directory path whatever: "<<_out_dir_<<std::endl;

  // pretty much replacing "/Desktop/spotting/build" with "/Desktop/spotting/data/"
  index = 0;
  while(index < _out_dir_.size() - 5)
  {
    index++;
  }
  _data_dir_[index+1] = 'd';
  _data_dir_[index+2] = 'a';
  _data_dir_[index+3] = 't';
  _data_dir_[index+4] = 'a';
  _data_dir_[index+5] = '/';

  // std::cout<<"checking the file directory path whatever: "<<_data_dir_<<std::endl<<std::endl<<std::endl;

  // we have 18 files total
  std::string* _files_ = new std::string[18];

  // looping through to copy the correct directory for each file
  for(int i = 0; i < 18; i++){
    if(i==14) // because the time of arrivals file is in the out folder
    {
      _files_[i] = _out_dir_;
    }
    else // and literally every other file is in the data folder
    {
      _files_[i] = _data_dir_;
    }
  }

  // here are the files, these are brief details more will be explained throughout the code and the thesis
  _files_[0]  += "default.fmd"; // this is the rothermel 13 and anderson 40 fuel model value sheet
  _files_[1]  += "kyle.fms"; // this is the fmd's in kyle canyon, this was made to best represent the vegetation
  
  if(_FUEL_ == 0)
  {
    _files_[2]  += "original.fuel"; // this is the default fmd's at a fixed resolution used to calculate the fire
  } // is Roger's and Jessie's fuel file they used and is the base case for kyle canyon
  else // if _FUEL_ == 1
  {
    _files_[2]  += "modified.fuel"; // this is the modified fmd's, i.e. bull dozer
  } // this is to show bulldozing and or firebreaks

  _files_[3]  += "fixed2.tif"; // this is the terrain data in a tif extracted by gdal
  _files_[4]  += "canopy_ht.asc"; // this is the canopy height at a fixed resolution
  _files_[5]  += "crown_base_ht.asc"; // this is the canopy base height at a fixed resolution
  _files_[6]  += "crown_bulk_density.asc"; // this is the canopy bulk density at a fixed resolution

  if(_WIND_XY_ == 0)
  {
    _files_[7]  += "originalX.fuel"; // this is original wind on the x-axis at a fixed resolution
    _files_[8]  += "originalY.fuel"; // this is original wind on the y-axis at a fixed resolution
  } // original files is all 0's, so no wind
  else
  {
    _files_[7]  += "modifiedX.fuel"; // this is modified wind on the x-axis at a fixed resolution
    _files_[8]  += "modifiedY.fuel"; // this is modified wind on the y-axis at a fixed resolution
  } // these files should have values to show wind has an affect

  _files_[9]  += "liveH.fuel";  // this is the dynamic live herbaceuous moisture at a fixed resolution
  _files_[10] += "dead1hr.fuel"; // this is the dynamic dead 1 hour moisture at a fixed resolution

  if(_START_TICK_ == 0)
  {
    _files_[11] += "csr.csv"; // this is the default current spread rates which are all set to 0
    _files_[12] += "msr.csv"; // this is the default max spread rates which are all set to 0
  }
  else
  {
    _files_[11] += "pcsr.csv"; // this is used to save and resume current spread rate values
    _files_[12] += "pmsr.csv"; // this is used to save and resume max spread rate values
  }

  _files_[13] += "fire_info.csv"; // this is the files that starts the fires initially 1 means ignite
  _files_[14] += "final_tests.csv"; // this is the time of arrivals at the end of the entire simulation
  _files_[15] += "paused.csv"; // this is used to save the time of arrivals for the paused simulation
  _files_[16] += "pmsr.csv"; // this is used to save and resume max spread rate values
  _files_[17] += "pcsr.csv"; // this is used to save and resume current spread rate values

  for(unsigned int i = 0; i < 1; i ++)
  {
    if(_PARALLEL_ == 1)
    {
      printf("---------- Running Parallel Simulation ----------\n");

      BD parallel_sim(    906,642, // sim dimension x and y
                         _files_[0],
                         _files_[1]);

      parallel_sim.Init( _files_[2],
                         _files_[3],
                         _files_[4],
                         _files_[5],
                         _files_[6],
                         _files_[7],
                         _files_[8],
                         _files_[9],
                         _files_[10],
                         _files_[11],
                         _files_[12],
                         _MOISTURE_,
                         _CROWNING_,
                         _SPOTTING_);

      // takes the fire_info.csv as an initial ignition to start the simulation 
      // or the paused.csv if we paused the simulation and are now resuming it at that point
      std::ifstream fin;
      if(_START_TICK_ == 0)
      {
        fin.open(_files_[13]); // _files_[13] += "fire_info.csv";
      }
      else
      {
        fin.open(_files_[15]); // _files_[15] += "paused.csv";
      }

      // just to be sure
      if(!fin)
      {
        std::cerr<<"Failed to open fire_info file !";
        return 1;
      }

      // this extracts the metadata needed for gdal to interpolate between various resolutions
      // please see thesis for an explanation
      char temp[255];
      double left_top_lat, left_top_long, right_bottom_lat, right_bottom_long;
      int numrows, numcols, lmaxval, notsetfire;
      std::string* metaData = new std::string[8];
      std::string spaceColon = ": ";

      fin.getline(temp, 255, ':');
      fin >> left_top_lat;
      metaData[0] = temp + spaceColon + double_to_string(left_top_lat);

      fin.getline(temp, 255, ':');
      fin >> left_top_long;
      metaData[1] = temp + spaceColon + double_to_string(left_top_long);

      fin.getline(temp, 255, ':');
      fin >> right_bottom_lat;
      metaData[2] = temp + spaceColon + double_to_string(right_bottom_lat);

      fin.getline(temp, 255, ':');
      fin >> right_bottom_long;
      metaData[3] = temp + spaceColon + double_to_string(right_bottom_long);

      fin.getline(temp, 255, ':');
      fin >> numrows;
      metaData[4] = temp + spaceColon + int_to_string(numrows);

      fin.getline(temp, 255, ':');
      fin >> numcols;
      metaData[5] = temp + spaceColon + int_to_string(numcols);

      fin.getline(temp, 255, ':');
      fin >> lmaxval;
      metaData[6] = temp + spaceColon + int_to_string(lmaxval);

      fin.getline(temp, 255, ':');
      fin >> notsetfire;
      metaData[7] = temp + spaceColon + int_to_string(notsetfire);

      // once we get all the meta data it is onto which cells are ignited and which are not
      // notice that fin is used above and below, look at the .csv files the order matters!
      int r = 0, c = 0;
      while(r < numrows)
      {
        int value;
        char dummy;
        while(c < numcols-1)
        {      
          fin >> value >> dummy;
          if(value != 0)
          { // this ignites the fire initially, this literally works no where else
            parallel_sim.UpdateCell(r,c,1);
          }
          c++;
        }
        fin >> value;
        if(value != 0)
        { // oh an another cool thing the number "1" is supposed to be a time of arrival but other values
          // don't do shit, whether the value is 1 or 1000 or anything inbetween the fire starts immediately
          parallel_sim.UpdateCell(r,c,1);
        }
        c = 0;
        r++;
      }
      fin.close();

      // time shit
      struct timeval start, finish;
      gettimeofday(&start, NULL);

      // transfer hella data to the gpu from the cpu
      parallel_sim.CopyToDevice(_MOISTURE_, _CROWNING_, _SPOTTING_);

      // pop off the kernel, the loop is for timing purposes
      for(unsigned int i = 0; i < 1; i++)
      {
        parallel_sim.RunKernel(_START_TICK_, s, s, _CROWNING_, _SPOTTING_, _STOP_TICK_);
      }

      // transer hella data from the gpu back to cpu
      parallel_sim.CopyFromDevice();

      // finish timing shit
      gettimeofday(&finish, NULL);
      t_upSpread = finish.tv_usec + finish.tv_sec * 1000000.0;
      t_upSpread -= start.tv_usec + start.tv_sec * 1000000.0;
      t_upSpread /= 1000000.0;
      std::cout << "Processing Simulation took " << t_upSpread << " seconds" << std::endl;

      // this happen whether we pause or don't because let's say we want to just run half
      // the simulation and then watch it and also if we do pause then resume it will be overwritten
      // which is good cause then no matter what we have the TOA's for every scenario
      parallel_sim.WriteToFile(_files_[14]);
      parallel_sim.WritePauseToFile(_files_[15], metaData, _files_[16], _files_[17]);

      // release memory
      delete[] metaData;
    }
    else // !_PARALLEL_ run the sequential that takes about 55 secs minimum to 90 secs
    {    // depending on what's turned on or off i.e. crowning, spotting
       printf("---------- Running Sequential Simulation ----------\n");
       SequentialSpread sequential_sim(906,642,_STOP_TICK_,_MOISTURE_,_CROWNING_,_SPOTTING_);
       sequential_sim.Init();
       sequential_sim.CalcMaxSpreadRates();
       struct timeval start, finish;
       gettimeofday(&start, NULL);
       sequential_sim.RunSimulationBD(step_size);
       gettimeofday(&finish, NULL);
       double t_upSpread = finish.tv_usec + finish.tv_sec * 1000000.0;
       t_upSpread -= start.tv_usec + start.tv_sec * 1000000.0;
       t_upSpread /= 1000000.0;
       std::cout << "Processing Simulation took " << t_upSpread << " seconds" << std::endl;
       sequential_sim.WriteToFile();
    }
  }
  std::cout << "Simulation Complete" << std::endl;
  return 0;
}