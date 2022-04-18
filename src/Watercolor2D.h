#ifndef WATERCOLOR2D_H_
#define WATERCOLOR2D_H_

#include <Eigen/Dense>
#include "Grid/StaggeredGrid.h"

class Watercolor2D
{
public:
  Watercolor2D(const int x_res, const int y_res);
  ~Watercolor2D() {};

  void step();
protected:
  void moveWater();
    void updateVelocities();
      void enforceBoundaryConditions();
    void relaxDivergence();
    void flowOutward();
  void movePigment();
  void transferPigment();
  void simulateCapillaryFlow();

  int _x_res, _y_res;
  float _dx;
  float _dt;

  // shallow-water layer
  Eigen::ArrayXXf _M;         // wet-area mask M
  StaggeredGrid _u;           // velocity in the x direction
  StaggeredGrid _v;           // velocity in the y direction
  Eigen::ArrayXXf _pressure;  // pressure
  // TODO: turn _g and _d into separate data structures for multiple pigments
  Eigen::ArrayXXf _g;         // pigment concentration (just 1)
  Eigen::ArrayXXf _h;         // paper height
  Eigen::ArrayXXf _delta_h;   // gradient of paper height
  float _viscosity;
  float _viscous_drag;

  // pigment-deposition layer
  Eigen::ArrayXXf _d;         // deposited pigment concentration (again, just 1)
  float _p_density;           // pigment density
  float _p_staining_power;    // pigment staining power
  float _p_granularity;       // pigment granularity

  // capillary layer
  Eigen::ArrayXXf _s;         // water saturation s of the paper
  Eigen::ArrayXXf _c;         // fluid-holding capacity c of the paper
};

#endif
