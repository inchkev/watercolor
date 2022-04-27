#include "Watercolor2D.h"
#include "util/GaussianBlur.h"
#include <iostream>

// http://demofox.org/gauss.html
// generated with sigma = 1.0, support = 0.995
const float GAUSSIAN_KERNEL[9][9] =
{
  {0.0000f, 0.0000f, 0.0000f, 0.0001f, 0.0001f, 0.0001f, 0.0000f, 0.0000f, 0.0000f},
  {0.0000f, 0.0000f, 0.0004f, 0.0014f, 0.0023f, 0.0014f, 0.0004f, 0.0000f, 0.0000f},
  {0.0000f, 0.0004f, 0.0037f, 0.0146f, 0.0232f, 0.0146f, 0.0037f, 0.0004f, 0.0000f},
  {0.0001f, 0.0014f, 0.0146f, 0.0584f, 0.0926f, 0.0584f, 0.0146f, 0.0014f, 0.0001f},
  {0.0001f, 0.0023f, 0.0232f, 0.0926f, 0.1466f, 0.0926f, 0.0232f, 0.0023f, 0.0001f},
  {0.0001f, 0.0014f, 0.0146f, 0.0584f, 0.0926f, 0.0584f, 0.0146f, 0.0014f, 0.0001f},
  {0.0000f, 0.0004f, 0.0037f, 0.0146f, 0.0232f, 0.0146f, 0.0037f, 0.0004f, 0.0000f},
  {0.0000f, 0.0000f, 0.0004f, 0.0014f, 0.0023f, 0.0014f, 0.0004f, 0.0000f, 0.0000f},
  {0.0000f, 0.0000f, 0.0000f, 0.0001f, 0.0001f, 0.0001f, 0.0000f, 0.0000f, 0.0000f},
};

Watercolor2D::Watercolor2D(const int x_res, const int y_res) :
  _x_res(x_res),
  _y_res(y_res),
  _u(x_res, y_res, 0),
  _v(x_res, y_res, 1)
{
  _dt = 0.001;
  _dx = 0.5;

  _M = Eigen::ArrayXXf::Zero(x_res, y_res);
  _pressure = Eigen::ArrayXXf::Zero(x_res, y_res);
  _g = Eigen::ArrayXXf::Zero(x_res, y_res);
  _h = Eigen::ArrayXXf::Zero(x_res, y_res);

  _viscosity = 0.1;
  _viscous_drag = 0.01;

  // one color model
  // TODO: struct for pigments
  _d = Eigen::ArrayXXf::Zero(x_res, y_res);
  _p_density = 0.02;
  _p_staining_power = 1.0;
  _p_granularity = 0.3;

  _s = Eigen::ArrayXXf::Zero(x_res, y_res);
  _c = Eigen::ArrayXXf::Zero(x_res, y_res);
}

void Watercolor2D::step()
{
  moveWater();
  movePigment();
  transferPigment();
  /* simulateCapillaryFlow(); */
}

void Watercolor2D::moveWater()
{
  updateVelocities();
  relaxDivergence();
  flowOutward();
}

/**
 * Update water velocities.
 */
void Watercolor2D::updateVelocities()
{
  float a, b;
  // _u -= _delta_h_x;
  // _v -= _delta_h_y;
  // for (int j = 0; j < _y_res; j++)
  //   for (int i = 0; i < _x_res; i++)
  //   {
  //     const float valx = ((float**)dhx_)[i][j] * 0.5f;
  //     const float valy = ((float**)dhy_)[i][j] * 0.5f;
  //     _u(i-0.5f,j) -= valx;
  //     _u(i+0.5f,j) -= valx;
  //     _v(i,j-0.5f) -= valy;
  //     _v(i,j+0.5f) -= valy;
  //   }

  _dt = 1.0f / std::max(std::abs(_u.max()), std::abs(_v.max()));
  StaggeredGrid u_new(_x_res, _y_res, 0);
  StaggeredGrid v_new(_x_res, _y_res, 1);
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
    {
      a = _u(i,j)*_u(i,j) - _u(i+1,j)*_u(i+1,j) +
        _u(i+0.5f,j-0.5f)*_v(i+0.5f,j-0.5f) -
        _u(i+0.5f,j+0.5f)*_v(i+0.5f,j+0.5f);
      b = _u(i+1.5f,j) + _u(i-0.5f,j) + _u(i+0.5f,j) + _u(i+0.5f,j-1) - 4*_u(i+0.5f,j);
      u_new(i+0.5f,j) = _u(i+0.5f,j) + _dt * (
          a - _viscosity*b + _pressure(i,j) - _pressure(i+1,j) - _viscous_drag*_u(i+0.5f,j));
      a = _v(i,j)*_v(i,j) - _u(i+1,j)*_u(i+1,j) +
        _u(i-0.5f,j+0.5f)*_v(i-0.5f,j+0.5f) -
        _u(i+0.5f,j+0.5f)*_v(i+0.5f,j+0.5f);
      b = _v(i+1,j+0.5f) + _v(i-1,j+0.5f) + _v(i,j+1.5f) + _v(i,j-0.5f) - 4*_v(i,j+0.5f);
      v_new(i,j+0.5f) = _u(i,j+0.5f) + _dt * (
          a - _viscosity*b + _pressure(i,j) - _pressure(i,j+1) - _viscous_drag*_u(i,j+0.5f));
    }
  _u = u_new;
  _v = v_new;
  enforceBoundaryConditions();
}

/**
 * Set velocities at the boundary of any pixel
 * not in the wet-area mask to zero.
 */
void Watercolor2D::enforceBoundaryConditions()
{
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
    {
      // if not wet, set velocities of cell to zero
      if (_M(i, j) == 1.0f)
        continue;

      _u(i+0.5f,j) = 0.0f;
      _u(i-0.5f,j) = 0.0f;
      _u(i,j+0.5f) = 0.0f;
      _u(i,j-0.5f) = 0.0f;
      _v(i+0.5f,j) = 0.0f;
      _v(i-0.5f,j) = 0.0f;
      _v(i,j+0.5f) = 0.0f;
      _v(i,j-0.5f) = 0.0f;
    }
}

/**
 * Relax the divergence of the velocity field after
 * each time step until it less than some tolerance
 * tau by redistributing the fluid into neighboring
 * grid cells.
 */
void Watercolor2D::relaxDivergence()
{
  const int N = 50;
  const float tau = 0.01f;
  const float xi = 0.1f;

  int t = 0;
  float delta_max = 0.0f;

  while (delta_max > tau || t >= N)
  {
    StaggeredGrid u_new(_u);
    StaggeredGrid v_new(_v);
    delta_max = 0.0f;
    for (int j = 0; j < _y_res; j++)
      for (int i = 0; i < _x_res; i++)
      {
        float delta = xi * (_u(i+0.5f,j) - _u(i-0.5f,j) + _v(i,j+0.5f) - _v(i,j-0.5f));
        _pressure(i,j) += delta;
        u_new(i+0.5f,j) += delta;
        u_new(i-0.5f,j) -= delta;
        v_new(i,j+0.5f) += delta;
        v_new(i,j-0.5f) -= delta;
        delta_max = std::max(std::abs(delta), delta_max);
      }
    _u = u_new;
    _v = v_new;
    t++;
  }
}

/**
 * Removes an amount of water from each cell according
 * to the cell's distance from the boundary of the
 * wet-area mask, with more water removed from cells closer
 * to the boundary.
 *
 * p <- p - η (1 - M') M
 * "In our examples, K = 10 and 0.01 <= η <= 0.05."
 */
void Watercolor2D::flowOutward()
{
  const float eta = 0.01f;
  const int gaussian_radius = 5; // K = 9 ish

  Eigen::ArrayXXf M_copy = _M;
  Eigen::ArrayXXf M_blur(_x_res, _y_res);
  approximateGaussianBlur(M_copy, M_blur, _x_res, _y_res, gaussian_radius);
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
      _pressure -= eta * (1.0f - M_blur(i, j)) * _M(i,j);
}

/**
 * Distributes pigment from each cell to its neighbors
 * according to the rate of fluid movement out of
 * the cell.
 */
void Watercolor2D::movePigment()
{
  _dt = 1.0f / std::max(std::abs(_u.max()), std::abs(_v.max()));
  // when more than one pigment, loop through each pigment
  Eigen::ArrayXXf g_new = _g;
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
    {
      const float gij = _g(i,j);
      g_new(i+1,j) += std::max(0.0f, _u(i+0.5f,j) * gij);
      g_new(i-1,j) += std::max(0.0f, -_u(i-0.5f,j) * gij);
      g_new(i,j+1) += std::max(0.0f, _v(i,j+0.5f) * gij);
      g_new(i,j-1) += std::max(0.0f, -_v(i,j-0.5f) * gij);
      g_new(i,j) += -std::max(0.0f, _u(i+0.5f,j) * gij) +
        std::max(0.0f, -_u(i-0.5f,j) * gij) +
        std::max(0.0f, _v(i,j+0.5f) * gij) +
        std::max(0.0f, -_v(i,j-0.5f) * gij);
    }
  _g = g_new;
}

/**
 * Controls pigment adsorption and desorption, and granulation.
 */
void Watercolor2D::transferPigment()
{
  float delta_up, delta_down;
  // when more than one pigment, loop through each pigment
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
    {
      if (_M(i,j) != 1.0f)
        continue;

      delta_down = _g(i,j) * (1.0f - _h(i,j) * _p_granularity) * _p_density;
      delta_up = _d(i,j) * (1.0f + (_h(i,j) - 1.0f) * _p_granularity) * _p_density / _p_staining_power;
      if ((_d(i,j) + delta_down) > 1.0f)
        delta_down = std::max(0.0f, 1.0f - _d(i,j));
      if ((_g(i,j) + delta_up) > 1.0f)
        delta_up = std::max(0.0f, 1.0f - _g(i,j));
      /* std::cout << delta_down << "," << delta_up << std::endl; */
      _d(i,j) += delta_down - delta_up;
      _g(i,j) += delta_up - delta_down;
    }
}

/**
 * Idk
 * Variation in cell capacity c from pixel to pixel results
 * in an irregular branching pattern
 */
void Watercolor2D::simulateCapillaryFlow()
{
  // what the fuck are all these thresholds and constants

  // alpha is absorption rate - I have no idea
  //
  const float alpha = 0.1;
  // epsilon, minimum saturation a pixel must have before it can diffuse
  // to its neighbors
  const float epsilon = 0.1;
  // delta, saturation value below which a pixel will not receive diffusion
  const float delta = 0.1;
  // sigma, saturation threshold which determines if wet-area mask _M is
  // expanded or not
  const float sigma = 0.1;

  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
      if (_M(i,j) > 0.0f)
        _s(i,j) += std::max(0.0f, std::min(alpha, _c(i,j) - _s(i,j)));

  Eigen::ArrayXXf s_new = _s;
  for (int j = 1; j < _y_res - 1; j++)
    for (int i = 1; i < _x_res - 1; i++)
      // for each cell (k, l) that are neijghbors of (i, j)
      for (int l = j - 1; l < j + 2; l += 2)
        for (int k = i - 1; k < i + 2; k += 2)
          if (_s(i,j) > epsilon && _s(i,j) > _s(k,l) && _s(k,l) > delta)
          {
            const float delta_s = std::max(0.0f,
              std::min(_s(i,j) - _s(k,l), _c(k,l) - _s(k,l)) / 4.0f);
            s_new(i,j) -= delta_s;
            s_new(k,l) += delta_s;
          }
  _s = s_new;

  // expand wet-area mask if cell's saturation exceeds a threshold sigma
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
      if (_s(i,j) > sigma)
        _M(i,j) = 1.0f;
}
