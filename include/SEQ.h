//
// Created by jsmith on 10/26/15.
//

#ifndef SIMULATOR_SEQ_H
#define SIMULATOR_SEQ_H
#include "FireSim.h"
#include <fstream>
#include "Ember.h"

class SequentialSpread{
public:
   SequentialSpread(int,int,int,int,int,int);
   ~SequentialSpread();
   bool Init();
   bool RunSimulationBD(int);
   bool RunSimulationIMT(int);
   bool RunSimulationMT(int);
   bool WriteToFile();
   bool CalcMaxSpreadRates();
   bool Accelerate();
   bool TestCrownRate(float);
   bool TestSpotting(float);
protected:
   FireSim* simulation_;
   float*   maxspreadrate_;
   float*   curspreadrate_;
   bool*    ember_map_;
   float*   test_map_;
   int sim_size_;
   int sim_rows_;
   int sim_cols_;
   std::vector<Ember> ember_list_;
};
#endif //SIMULATOR_SEQ_H
