#include "Watercolor2D.h"
#include "util/GaussianBlur.h"
#include <iostream>

Watercolor2D::Watercolor2D(const int x_res, const int y_res) :
  _x_res(x_res),
  _y_res(y_res),
  _u(x_res, y_res, 0),
  _v(x_res, y_res, 1)
{
  _dt = 0.01;
  _dx = 0.5;

  _buffer = new float[3 * _x_res * _y_res];
  _paper = new float[3 * _x_res * _y_res];

  _M = Eigen::ArrayXXf::Zero(x_res, y_res);
  _pressure = Eigen::ArrayXXf::Zero(x_res, y_res);
  _h = Eigen::ArrayXXf::Zero(x_res, y_res);
  _dhx = Eigen::ArrayXXf::Zero(x_res, y_res);
  _dhy = Eigen::ArrayXXf::Zero(x_res, y_res);

  _viscosity = 0.1;
  _viscous_drag = 0.01;

  // TODO: struct for pigments
  Pigment* red = new Pigment();
  red->g = Eigen::ArrayXXf::Zero(x_res, y_res);
  red->d = Eigen::ArrayXXf::Zero(x_res, y_res);
  red->K = Eigen::Vector3f(0.22f, 1.47f, 0.57f);
  red->S = Eigen::Vector3f(0.05f, 0.003f, 0.03f);
  red->density = 0.02;
  red->staining_power = 5.5;
  red->granularity = 0.81;
  _pigments.push_back(red);

  _s = Eigen::ArrayXXf::Zero(x_res, y_res);
  _c = Eigen::ArrayXXf::Zero(x_res, y_res);
}

void Watercolor2D::setPaper(float*& paper)
{
  for (int x = 0; x < _x_res; x++)
    for (int y = 0; y < _y_res; y++)
    {
      float sum = 0.0f;
      const int index = 3 * (x * _y_res + y);
      for (int i = 0; i < 3; i++)
      {
        _paper[index + i] = paper[index + i];
        _buffer[index + i] = paper[index + i];
        sum += paper[index + i];
      }
      _h(x,y) = sum / 3;
    }

  // calculate gradients
  _dhx = _h.rightCols(_h.cols() - 1) - _h.leftCols(_h.cols() - 1);
  _dhy = _h.bottomRows(_h.rows() - 1) - _h.topRows(_h.rows() - 1);

  // ugly hack to make both 200x200
  // TODO: cleaner boundary condition handling
  _dhx.resize(_x_res, _y_res);
  _dhy.resize(_x_res, _y_res);
}

void Watercolor2D::step()
{
  std::cout << "-----" << std::endl;
  std::cout << "   " << _u.max() << " " << _v.max() << std::endl;
  std::cout << "-----updateVelocities" << std::endl;
  updateVelocities();
  std::cout << "   " << _u.max() << " " << _v.max() << std::endl;
  std::cout << "-----relaxDivergence" << std::endl;
  relaxDivergence();
  std::cout << "   " << _u.max() << " " << _v.max() << std::endl;
  std::cout << "-----movePigment" << std::endl;

  movePigment();
  std::cout << "   " << _u.max() << " " << _v.max() << std::endl;
  std::cout << "-----transferPigment" << std::endl;
  transferPigment();
  std::cout << "   " << _u.max() << " " << _v.max() << std::endl;
  flowOutward();
  simulateCapillaryFlow();
  render();
}

/**
 * Update water velocities.
 */
void Watercolor2D::updateVelocities()
{
  float a, b;
  _u -= _dhx;
  _v -= _dhy;

  /* _dt = std::max(1.0f / std::max(_u.absmax(), _v.absmax()), 0.01f); */
  _dt = 0.1f / std::max(_u.absmax(), _v.absmax());
  /* std::cout << "  _dt " << _dt << std::endl; */
  for (float t = 0.0f; t < 0.1f; t += _dt)
  {
    StaggeredGrid u_new(_x_res, _y_res, 0);
    StaggeredGrid v_new(_x_res, _y_res, 1);
    for (int j = 1; j < _y_res-1; j++)
      for (int i = 1; i < _x_res-1; i++)
      {
        /* std::cout << "b "; */
        a = _u.get(i,j)*_u.get(i,j) - _u.get(i+1,j)*_u.get(i+1,j) +
          _u.get(i+0.5f,j-0.5f)*_v.get(i+0.5f,j-0.5f) -
          _u.get(i+0.5f,j+0.5f)*_v.get(i+0.5f,j+0.5f);
        b = _u.get(i+1.5f,j) + _u.get(i-0.5f,j) + _u.get(i+0.5f,j+1) + _u.get(i+0.5f,j-1) - 4.0f*_u.get(i+0.5f,j);
        /* std::cout << b << std::endl; */
        u_new(i+0.5f,j) = _u.get(i+0.5f,j) + _dt * (
        /* u_new(i+0.5f,j) += _dt * ( */
            a - _viscosity*b + _pressure(i,j) - _pressure(i+1,j) - _viscous_drag*_u.get(i+0.5f,j));
        a = _v.get(i,j)*_v.get(i,j) - _v.get(i,j+1)*_v.get(i,j+1) +
          _u.get(i-0.5f,j+0.5f)*_v.get(i-0.5f,j+0.5f) -
          _u.get(i+0.5f,j+0.5f)*_v.get(i+0.5f,j+0.5f);
        b = _v.get(i+1,j+0.5f) + _v.get(i-1,j+0.5f) + _v.get(i,j+1.5f) + _v.get(i,j-0.5f) - 4.0f*_v.get(i,j+0.5f);
        /* std::cout << b << std::endl; */
        v_new(i,j+0.5f) = _v.get(i,j+0.5f) + _dt * (
        /* v_new(i,j+0.5f) += + _dt * ( */
            a - _viscosity*b + _pressure(i,j) - _pressure(i,j+1) - _viscous_drag*_v.get(i,j+0.5f));
      }
    /* std::cout << "updatevel " << u_new.absmax() << " " << v_new.absmax() << std::endl; */
    /* std::cout << _pressure.maxCoeff() << std::endl; */
    _u = u_new;
    _v = v_new;
    enforceBoundaryConditions();
  }
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
  float sum = 0.0f;
  std::cout << "pressure " << _pressure.maxCoeff() << std::endl;

  StaggeredGrid u_new(_u);
  StaggeredGrid v_new(_v);
  while (1)
  {
    delta_max = 0.0f;
    u_new = _u;
    v_new = _v;
    for (int j = 0; j < _y_res; j++)
      for (int i = 0; i < _x_res; i++)
      {
        float delta = -xi * (_u.get(i+0.5f,j) - _u.get(i-0.5f,j) +
                             _v.get(i,j+0.5f) - _v.get(i,j-0.5f));
        _pressure(i,j) += delta;
        u_new(i+0.5f,j) += delta;
        u_new(i-0.5f,j) -= delta;
        v_new(i,j+0.5f) += delta;
        v_new(i,j-0.5f) -= delta;
        sum += delta;
        delta_max = std::max(std::abs(delta), delta_max);
      }
    t++;
    _u = u_new;
    _v = v_new;

    if (delta_max <= tau || t >= N)
      break;
  }
  std::cout << "sum  " << sum << std::endl;
  std::cout << "dmax " << delta_max << std::endl;
  std::cout << "pressure " << _pressure.maxCoeff() << std::endl;
}

/**
 * Removes an amount of water from each cell according
 * to the cell's distance from the boundary of the
 * wet-area mask, with more water removed from cells closer
 * to the boundary. Darkens edges.
 *
 * For backrun effects only.
 *
 * p <- p - η (1 - M') M
 * "In our examples, K = 10 and 0.01 <= η <= 0.05."
 */
void Watercolor2D::flowOutward()
{
  const float eta = 0.01f;
  const int gaussian_radius = 4; // K = 9 ish

  Eigen::ArrayXXf M_copy = _M;
  Eigen::ArrayXXf M_blur(_x_res, _y_res);
  approximateGaussianBlur(M_copy, M_blur, _x_res, _y_res, gaussian_radius);
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
    {
      const float blur = eta * (1.0f - M_blur(i, j)) * _M(i,j);
      _pressure(i,j) -= std::max(std::min(blur, _pressure(i,j)), 0.0f);
    }
}

/**
 * Distributes pigment from each cell to its neighbors
 * according to the rate of fluid movement out of
 * the cell.
 */
void Watercolor2D::movePigment()
{
  /* _dt = 1.0f / std::max(_u.absmax(), _v.absmax()); */
  _dt = 0.1f / std::max(_u.absmax(), _v.absmax());
  std::cout << "movePigment _dt = " << _dt << std::endl;
  // when more than one pigment, loop through each pigment
  for (Pigment* pig: _pigments)
  {
    for (float t = 0.0f; t < 0.1f; t += _dt)
    {
      Eigen::ArrayXXf g_new = pig->g;
      for (int j = 1; j < _y_res-1; j++)
        for (int i = 1; i < _x_res-1; i++)
        {
          const float gij = pig->g(i,j);
          g_new(i+1,j) += std::max(0.0f, _u.get(i+0.5f,j) * gij);
          g_new(i-1,j) += std::max(0.0f, -_u.get(i-0.5f,j) * gij);
          g_new(i,j+1) += std::max(0.0f, _v.get(i,j+0.5f) * gij);
          g_new(i,j-1) += std::max(0.0f, -_v.get(i,j-0.5f) * gij);
          g_new(i,j) += -std::max(0.0f, _u.get(i+0.5f,j) * gij) +
            std::max(0.0f, -_u.get(i-0.5f,j) * gij) +
            std::max(0.0f, _v.get(i,j+0.5f) * gij) +
            std::max(0.0f, -_v.get(i,j-0.5f) * gij);
        }
      pig->g = g_new;
    }
  }
}

/**
 * Controls pigment adsorption and desorption, and granulation.
 */
void Watercolor2D::transferPigment()
{
  float delta_up, delta_down;
  for (Pigment* pig: _pigments)
  {
    for (int j = 1; j < _y_res-1; j++)
      for (int i = 1; i < _x_res-1; i++)
      {
        if (_M(i,j) != 1.0f)
          continue;

        delta_down = pig->g(i,j) * (1.0f - _h(i,j) * pig->granularity) * pig->density;
        delta_up = pig->d(i,j) * (1.0f + (_h(i,j) - 1.0f) * pig->granularity) * pig->density / pig->staining_power;
        if ((pig->d(i,j) + delta_down) > 1.0f)
          delta_down = std::max(0.0f, 1.0f - pig->d(i,j));
        if ((pig->g(i,j) + delta_up) > 1.0f)
          delta_up = std::max(0.0f, 1.0f - pig->g(i,j));
        pig->d(i,j) += delta_down - delta_up;
        pig->g(i,j) += delta_up - delta_down;
      }
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
  const float sigma = 0.01;

  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
      if (_M(i,j) > 0.0f)
        _s(i,j) += std::max(0.0f, std::min(alpha, _c(i,j) - _s(i,j)));

  Eigen::ArrayXXf s_new = _s;
  for (int j = 1; j < _y_res - 1; j++)
    for (int i = 1; i < _x_res - 1; i++)
      // for each cell (k, l) that are neighbors of (i, j)
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

void Watercolor2D::render()
{
  
}
