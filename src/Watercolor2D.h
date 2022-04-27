#ifndef WATERCOLOR2D_H_
#define WATERCOLOR2D_H_

#include <Eigen/Dense>
#include "Grid/StaggeredGrid.h"

class Watercolor2D
{
  // TODO: create a pigment class
public:
  Watercolor2D(const int x_res, const int y_res);
  ~Watercolor2D() {};

  void step();

  int& x_res()            { return _x_res; };
  const int x_res() const { return _x_res; };
  int& y_res()            { return _y_res; };
  const int y_res() const { return _y_res; };
  Eigen::ArrayXXf& g()            { return _g; };
  const Eigen::ArrayXXf g() const { return _g; };
  Eigen::ArrayXXf& d()            { return _d; };
  const Eigen::ArrayXXf d() const { return _d; };
  Eigen::ArrayXXf& M()            { return _M; };
  const Eigen::ArrayXXf M() const { return _M; };
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
  Eigen::ArrayXXf _delta_h_x; // gradient of paper height
  Eigen::ArrayXXf _delta_h_y; // gradient of paper height
  float _viscosity;           // viscosity mu
  float _viscous_drag;        // viscous drag kappa

  // pigment-deposition layer
  Eigen::ArrayXXf _d;         // deposited pigment concentration (again, just 1)
  float _p_density;           // pigment density rho
  float _p_staining_power;    // pigment staining power omega
  float _p_granularity;       // pigment granularity gamma

  // capillary layer
  Eigen::ArrayXXf _s;         // water saturation s of the paper
  Eigen::ArrayXXf _c;         // fluid-holding capacity c of the paper

  void *dhx_;
  void *dhy_;
};

#endif
