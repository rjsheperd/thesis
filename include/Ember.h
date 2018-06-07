//
// Created by jsmith on 4/4/16.
//

#ifndef SIMULATOR_EMBER_H
#define SIMULATOR_EMBER_H


class Ember {
public:
   Ember(int _x, int _y, float _z, float _max, double _ang, double _mag)
         {x_ = _x; y_ = _y; z_ = _z; z_o_ = _max,
          angle_ = _ang; magnitude_ = _mag;};
   int x_; // col
   int y_; // row
   float z_; // current height
   float z_o_; // max height (starting point for descent)
   double angle_;
   double magnitude_;
   float ember_toa = 32767;
};


#endif //SIMULATOR_EMBER_H
