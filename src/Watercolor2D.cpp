#include "Watercolor2D.h"

Watercolor2D::Watercolor2D(const int x_res, const int y_res) :
  _x_res(x_res),
  _y_res(y_res),
  _M(x_res, y_res),
  _u(x_res, y_res, 0),
  _v(x_res, y_res, 1),
  _pressure(x_res, y_res),
  _g(x_res, y_res),
  _h(x_res, y_res),
  _delta_h(x_res, y_res),
  _d(x_res, y_res),
  _s(x_res, y_res),
  _c(x_res, y_res)
{
  _dt = 0.01;
  _dx = 0.5;
}

void Watercolor2D::step()
{
  moveWater();
  movePigment();
  transferPigment();
  simulateCapillaryFlow();
}

void Watercolor2D::moveWater()
{
  updateVelocities();
  relaxDivergence();
  flowOutward();
}

void Watercolor2D::updateVelocities()
{
}

void Watercolor2D::enforceBoundaryConditions()
{
}

void Watercolor2D::relaxDivergence()
{
}

void Watercolor2D::flowOutward()
{
}

void Watercolor2D::movePigment()
{
}

void Watercolor2D::transferPigment()
{
}

void Watercolor2D::simulateCapillaryFlow()
{
}
