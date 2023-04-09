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

  _buffer = new Eigen::Vector3f*[_x_res];
  for (int i = 0; i < _x_res; ++i)
    _buffer[i] = new Eigen::Vector3f[_y_res];

  _paper = new Eigen::Vector3f*[_x_res];
  for (int i = 0; i < _x_res; ++i)
    _paper[i] = new Eigen::Vector3f[_y_res];

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
  red->color = Eigen::Vector3f(0.83f, 0.2f, 0.21f);
  red->density = 0.02f;
  red->staining_power = 5.5f;
  red->granularity = 0.81f;
  _pigments.push_back(red);

  Pigment* blue = new Pigment();
  blue->g = Eigen::ArrayXXf::Zero(x_res, y_res);
  blue->d = Eigen::ArrayXXf::Zero(x_res, y_res);
  blue->K = Eigen::Vector3f(0.86f, 0.86f, 0.06f);
  blue->S = Eigen::Vector3f(0.005f, 0.005f, 0.09f);
  blue->color = Eigen::Vector3f(0.2f, 0.3f, 0.80f);
  blue->density = 0.01f;
  blue->staining_power = 3.1f;
  blue->granularity = 0.91f;
  _pigments.push_back(blue);

  for (Pigment* pig: _pigments) {
    for (int i = 0; i < 3; i++)
    {
      pig->a[i] = 1.0f + pig->K[i]/pig->S[i];
      pig->b[i] = sqrtf(pig->a[i]*pig->a[i] - 1.0f);
    }
  }

  _s = Eigen::ArrayXXf::Zero(x_res, y_res);
  _c = Eigen::ArrayXXf::Zero(x_res, y_res);
  /* _c_max = 0.8f; */
  /* _c_min = 0.05f; */
  _c_max = 0.8f;
  _c_min = 0.05f;
}

void Watercolor2D::setPaper(float*& paper, const int paper_x, const int paper_y)
{
  for (int x = 0; x < _x_res; x++)
    for (int y = 0; y < _y_res; y++)
    {
      float sum = 0.0f;
      const int index = 3 * (x * paper_y + y);
      for (int i = 0; i < 3; i++)
      {
        _paper[x][y][i] = paper[index + i];
        _buffer[x][y][i] = paper[index + i];
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

  // update fluid capacity c
  _c = _h * (_c_max - _c_min) + _c_min;
}

void Watercolor2D::step()
{
  // updateVelocities();
  // relaxDivergence();
  // flowOutward();
  // movePigment();
  // transferPigment();
  // simulateCapillaryFlow();
  // render();
  std::cout << "=====STEP START=====" << std::endl;
  std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl;
  std::cout << "d:" << _pigments[0]->d.maxCoeff() << " ";
  std::cout << "dsum:" << _pigments[0]->d.sum() << " | ";
  std::cout << "g:" << _pigments[0]->g.maxCoeff() << " ";
  std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl;
  std::cout << "-----updateVelocities" << std::endl;
  updateVelocities();
  std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl;
  std::cout << "d:" << _pigments[0]->d.maxCoeff() << " ";
  std::cout << "dsum:" << _pigments[0]->d.sum() << " | ";
  std::cout << "g:" << _pigments[0]->g.maxCoeff() << " ";
  std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl;
  std::cout << "-----relaxDivergence" << std::endl;
  relaxDivergence();
  std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl;
  std::cout << "d:" << _pigments[0]->d.maxCoeff() << " ";
  std::cout << "dsum:" << _pigments[0]->d.sum() << " | ";
  std::cout << "g:" << _pigments[0]->g.maxCoeff() << " ";
  std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl;
  std::cout << "-----flowOutward" << std::endl;
  flowOutward();
  std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl;
  std::cout << "d:" << _pigments[0]->d.maxCoeff() << " ";
  std::cout << "dsum:" << _pigments[0]->d.sum() << " | ";
  std::cout << "g:" << _pigments[0]->g.maxCoeff() << " ";
  std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl;
  std::cout << "-----movePigment" << std::endl;
  movePigment();
  std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl;
  std::cout << "d:" << _pigments[0]->d.maxCoeff() << " ";
  std::cout << "dsum:" << _pigments[0]->d.sum() << " | ";
  std::cout << "g:" << _pigments[0]->g.maxCoeff() << " ";
  std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl;
  std::cout << "-----transferPigment" << std::endl;
  transferPigment();
  std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl;
  std::cout << "d:" << _pigments[0]->d.maxCoeff() << " ";
  std::cout << "dsum:" << _pigments[0]->d.sum() << " | ";
  std::cout << "g:" << _pigments[0]->g.maxCoeff() << " ";
  std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl;
  std::cout << "-----simulateCapillaryFlow" << std::endl;
  simulateCapillaryFlow();
  std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl;
  std::cout << "d:" << _pigments[0]->d.maxCoeff() << " ";
  std::cout << "dsum:" << _pigments[0]->d.sum() << " | ";
  std::cout << "g:" << _pigments[0]->g.maxCoeff() << " ";
  std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl;
  std::cout << "=====" << std::endl;
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

  _dt = 0.4f / ceil(std::max(_u.absmax(), _v.absmax()));

  for (float t = 0.0f; t < 0.4f; t += _dt)
  {
    StaggeredGrid u_new(_u);
    StaggeredGrid v_new(_v);

    for (int i = 1; i < _x_res - 1; ++i)
      for (int j = 1; j < _y_res - 1; ++j)
      {
        const float pressureij = _pressure(i, j);
        const float uij = _u.get(i, j);
        const float uip1j = _u.get(i+1, j);
        const float uip05j = _u.get(i+0.5f, j);
        const float uip05jp05 = _u.get(i+0.5f, j+0.5f);
        const float vij = _v.get(i, j);
        const float vijp1 = _v.get(i, j+1);
        const float vijp05 = _v.get(i, j+0.5f);
        const float vip05jp05 = _v.get(i+0.5f, j+0.5f);

        a = uij*uij - uip1j*uip1j + _u.get(i+0.5f,j-0.5f)*_v.get(i+0.5f,j-0.5f) - uip05jp05*vip05jp05;
        b = _u.get(i+1.5f,j) + _u.get(i-0.5f,j) + _u.get(i+0.5f,j+1) + _u.get(i+0.5f,j-1) - 4.0f*uip05j;

        u_new(i+0.5f,j) = uip05j + _dt * (
            a - _viscosity*b + pressureij - _pressure(i+1,j) - _viscous_drag*uip05j);
        u_new(i+0.5f,j) = std::max(std::min(u_new(i+0.5f,j), 50.0f), -50.0f);

        a = vij*vij - vijp1*vijp1 + _u.get(i-0.5f,j+0.5f)*_v.get(i-0.5f,j+0.5f) - uip05jp05*vip05jp05;
        b = _v.get(i+1,j+0.5f) + _v.get(i-1,j+0.5f) + _v.get(i,j+1.5f) + _v.get(i,j-0.5f) - 4.0f*vijp05;

        v_new(i,j+0.5f) = vijp05 + _dt * (
            a - _viscosity*b + pressureij - _pressure(i,j+1) - _viscous_drag*vijp05);
        v_new(i,j+0.5f) = std::max(std::min(v_new(i,j+0.5f), 50.0f), -50.0f);
      }

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
  for (int i = 0; i < _x_res; ++i)
    for (int j = 0; j < _y_res; ++j)
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
  const float xi = -0.1f;

  int t = 0;
  float delta_max = 0.0f;

  StaggeredGrid u_new(_u);
  StaggeredGrid v_new(_v);
  do {
    u_new = _u;
    v_new = _v;
    delta_max = 0.0f;
    for (int i = 0; i < _x_res; ++i)
      for (int j = 0; j < _y_res; ++j)
      {
        const float delta = xi * (_u.get(i+0.5f,j) - _u.get(i-0.5f,j) +
                                  _v.get(i,j+0.5f) - _v.get(i,j-0.5f));
        _pressure(i,j) += delta;
        u_new(i+0.5f,j) += delta;
        u_new(i-0.5f,j) -= delta;
        v_new(i,j+0.5f) += delta;
        v_new(i,j-0.5f) -= delta;
        delta_max = std::max(std::abs(delta), delta_max);
      }
    _u = u_new;
    _v = v_new;
    ++t;
  } while (delta_max > tau && t < N);
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

  for (int i = 1; i < _x_res - 1; ++i)
    for (int j = 1; j < _y_res - 1; ++j) {
      const float blur = eta * (1.0f - M_blur(i,j)) * _M(i,j);
      _pressure(i,j) -= blur;
    }
}

/**
 * Distributes pigment from each cell to its neighbors
 * according to the rate of fluid movement out of
 * the cell.
 */
void Watercolor2D::movePigment()
{
  _dt = 0.4f / ceil(std::max(_u.absmax(), _v.absmax()));

  // when more than one pigment, loop through each pigment
  for (Pigment* pig: _pigments)
  {
    for (float t = 0.0f; t < 0.4f; t += _dt)
    {
      Eigen::ArrayXXf g_new = pig->g;

      for (int i = 1; i < _x_res - 1; ++i)
        for (int j = 1; j < _y_res - 1; ++j)
        {
          const float gij = pig->g(i,j);
          const float uip05j = _u.get(i+0.5f,j);
          const float uim05j = -_u.get(i-0.5f,j);
          const float vijp05 = _v.get(i,j+0.5f);
          const float vijm05 = -_v.get(i,j-0.5f);

          g_new(i+1,j) += _dt*std::max(0.0f, uip05j * gij);
          g_new(i-1,j) += _dt*std::max(0.0f, uim05j * gij);
          g_new(i,j+1) += _dt*std::max(0.0f, vijp05 * gij);
          g_new(i,j-1) += _dt*std::max(0.0f, vijm05 * gij);
          g_new(i,j) -= _dt*(
              std::max(0.0f, uip05j * gij) +
              std::max(0.0f, uim05j * gij) +
              std::max(0.0f, vijp05 * gij) +
              std::max(0.0f, vijm05 * gij));
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
    for (int i = 1; i < _x_res - 1; ++i)
      for (int j = 1; j < _y_res - 1; ++j)
      {
        if (_M(i,j) != 1.0f)
          continue;

        delta_down = pig->g(i,j) * (1.0f - _h(i,j) * pig->granularity) * pig->density;
        delta_up = pig->d(i,j) * (1.0f + (_h(i,j) - 1.0f) * pig->granularity) * pig->density / pig->staining_power;

        if ((pig->d(i,j) + delta_down) > 1.0f)
          delta_down = std::max(0.0f, 1.0f - pig->d(i,j));
        if ((pig->g(i,j) + delta_up) > 1.0f)
          delta_up = std::max(0.0f, 1.0f - pig->g(i,j));

        const float delta = delta_down - delta_up;
        pig->d(i,j) += delta;
        pig->g(i,j) -= delta;
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
  // these constants are not specified in the paper :(

  // alpha is absorption rate - not sure what is appropriate
  const float alpha = 0.10;

  // epsilon, minimum saturation a pixel must have before it can diffuse
  // to its neighbors
  const float epsilon = 0.72;

  // delta, saturation value below which a pixel will not receive diffusion
  const float delta = -0.01;

  // sigma, saturation threshold which determines if wet-area mask _M is
  // expanded or not
  const float sigma = 0.1;

  for (int i = 0; i < _x_res; ++i)
    for (int j = 0; j < _y_res; ++j)
      if (_M(i,j) > 0.0f)
        _s(i,j) += std::max(0.0f, std::min(alpha, _c(i,j) - _s(i,j)));

  Eigen::ArrayXXf s_new = _s;
  for (int j = 1; j < _y_res - 1; ++j)
    for (int i = 1; i < _x_res - 1; ++i)
    {
      // for each cell (k, l) that are neighbors of (i, j)
      int k = i-1;
      int l = j;
      if (_s(i,j) > epsilon && _s(i,j) > _s(k,l) && _s(k,l) > delta)
      {
        const float delta_s = std::max(0.0f,
          std::min(_s(i,j) - _s(k,l), _c(k,l) - _s(k,l)) / 4.0f);
        s_new(i,j) -= delta_s;
        s_new(k,l) += delta_s;
      }
      k = i+1;
      l = j;
      if (_s(i,j) > epsilon && _s(i,j) > _s(k,l) && _s(k,l) > delta)
      {
        const float delta_s = std::max(0.0f,
          std::min(_s(i,j) - _s(k,l), _c(k,l) - _s(k,l)) / 4.0f);
        s_new(i,j) -= delta_s;
        s_new(k,l) += delta_s;
      }
      k = i;
      l = j-1;
      if (_s(i,j) > epsilon && _s(i,j) > _s(k,l) && _s(k,l) > delta)
      {
        const float delta_s = std::max(0.0f,
          std::min(_s(i,j) - _s(k,l), _c(k,l) - _s(k,l)) / 4.0f);
        s_new(i,j) -= delta_s;
        s_new(k,l) += delta_s;
      }
      k = i;
      l = j+1;
      if (_s(i,j) > epsilon && _s(i,j) > _s(k,l) && _s(k,l) > delta)
      {
        const float delta_s = std::max(0.0f,
          std::min(_s(i,j) - _s(k,l), _c(k,l) - _s(k,l)) / 4.0f);
        s_new(i,j) -= delta_s;
        s_new(k,l) += delta_s;
      }
    }
  _s = s_new;
  std::cout << _s.abs().maxCoeff() << std::endl;

  // expand wet-area mask if cell's saturation exceeds a threshold sigma
  for (int j = 0; j < _y_res; ++j)
    for (int i = 0; i < _x_res; ++i)
      if (_M(i,j) < 1.0f && _s(i,j) > sigma) {
        _M(i,j) = 1.0f;
      }
}

void Watercolor2D::render()
{
  for (int i = 0; i < _x_res; ++i)
    for (int j = 0; j < _y_res; ++j)
    {
      /* Eigen::Vector3f paper = _paper[j][i]; // UNCOMMENT ME */
      Eigen::Vector3f paper(1.0f, 1.0f, 1.0f); // white background
      float total_thickness = 0.0f;
      for (Pigment* pig: _pigments)
      {
        float thickness =
          pig->d(i,j) + pig->g(i,j);
        thickness = std::min(std::max(0.0f, thickness), 1.0f);
        total_thickness += thickness;

        const float opacity = thickness;
        Eigen::Vector3f newcol =
          (1.0-opacity)*paper + opacity*pig->color;
        _buffer[j][i] = newcol;
        paper = newcol;
      }
      continue;

      // if no ink, continue
      if (total_thickness == 0.0f)
        continue;

      for (Pigment* pig: _pigments)
      {
        float thickness =
          pig->d(i,j) + pig->g(i,j);
        thickness = std::min(std::max(0.0f, thickness), 1.0f);
        const float x = thickness / total_thickness;
        for (int i = 0; i < 3; i++)
        {
          const float c = pig->a[i] * sinh(pig->b[i] * pig->S[i] * x) +
            pig->b[i] * cosh(pig->b[i] * pig->S[i] * x);
          pig->R[i] = sinh(pig->b[i] * pig->S[i] * x / c);
          pig->T[i] = pig->b[i] / c;
        }
      }

      // for now, can only properly blend 2 pigments
      Eigen::Vector3f R;
      for (int i = 0; i < 3; ++i)
      {
        R[i] = _pigments[0]->R[i] + (
            (_pigments[0]->T[i] * _pigments[0]->T[i] * _pigments[1]->R[i]) /
            (1.0f - _pigments[0]->R[i] * _pigments[1]->R[i])
        );
      }
      /* R = 1 - R; */

      const float opacity = std::min(total_thickness, 1.0f);
      Eigen::Vector3f newcol = opacity * R + (1.0f - opacity) * paper;
      _buffer[j][i] = newcol;
      /* _buffer[j][i] -= R; */
      /* paper = newcol; */
    }
}
