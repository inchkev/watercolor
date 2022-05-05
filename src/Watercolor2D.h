#ifndef WATERCOLOR2D_H_
#define WATERCOLOR2D_H_

#include <Eigen/Dense>
#include <vector>
#include "Grid/StaggeredGrid.h"

struct Pigment
{
  Eigen::ArrayXXf g;    // pigment concentration
  Eigen::ArrayXXf d;    // deposited pigment concentration
  Eigen::Vector3f K;    // K coefficient in KM theory
  Eigen::Vector3f S;    // S coefficient in KM theory
  float density;        // pigment density rho
  float staining_power; // pigment staining power omega
  float granularity;    // pigment granularity gamma
};

class Watercolor2D
{
  // TODO: create a pigment class
public:
  Watercolor2D(const int x_res, const int y_res);
  ~Watercolor2D()
  {
    for (Pigment* p: _pigments)
      delete p;
    _pigments.clear();
    for (int i = 0; i < _x_res; i++)
    {
      delete[] _buffer[i];
      delete[] _paper[i];
    }
    delete[] _buffer;
    delete[] _paper;
  };

  Eigen::Vector3f**& frameBuffer() { return _buffer; };

  void setPaper(float*& paper);
  void step();

  int& x_res()            { return _x_res; };
  const int x_res() const { return _x_res; };
  int& y_res()            { return _y_res; };
  const int y_res() const { return _y_res; };

  std::vector<Pigment*>& pigments()             { return _pigments; };
  const std::vector<Pigment*> pigments() const  { return _pigments; };
  Eigen::ArrayXXf& M()            { return _M; };
  const Eigen::ArrayXXf M() const { return _M; };
protected:
  void updateVelocities();
    void enforceBoundaryConditions();
  void relaxDivergence();
  void flowOutward();
  void movePigment();
  void transferPigment();
  void simulateCapillaryFlow();
  void render();

  int _x_res, _y_res;
  float _dt;

  // shallow-water layer
  Eigen::ArrayXXf _M;         // wet-area mask M
  StaggeredGrid _u;           // velocity in the x direction
  StaggeredGrid _v;           // velocity in the y direction
  Eigen::ArrayXXf _pressure;  // pressure

  Eigen::ArrayXXf _h;         // paper height
  Eigen::ArrayXXf _dhx;       // gradient of paper height
  Eigen::ArrayXXf _dhy;       // gradient of paper height

  float _viscosity;           // viscosity mu
  float _viscous_drag;        // viscous drag kappa

  // pigment-deposition layer
  std::vector<Pigment*> _pigments;

  // capillary layer
  Eigen::ArrayXXf _s;         // water saturation s of the paper
  Eigen::ArrayXXf _c;         // fluid-holding capacity c of the paper

  Eigen::Vector3f** _buffer;
  Eigen::Vector3f** _paper;
};

#endif
