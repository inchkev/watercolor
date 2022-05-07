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
  for (int i = 0; i < _x_res; i++)
    _buffer[i] = new Eigen::Vector3f[_y_res];

  _paper = new Eigen::Vector3f*[_x_res];
  for (int i = 0; i < _x_res; i++)
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
  red->color = Eigen::Vector3f(0.81f, 0.3f, 0.44f);
  red->density = 0.02;
  red->staining_power = 5.5;
  red->granularity = 0.61;
  _pigments.push_back(red);

  Pigment* blue = new Pigment();
  blue->g = Eigen::ArrayXXf::Zero(x_res, y_res);
  blue->d = Eigen::ArrayXXf::Zero(x_res, y_res);
  blue->K = Eigen::Vector3f(0.22f, 1.47f, 0.57f);
  blue->S = Eigen::Vector3f(0.05f, 0.003f, 0.03f);
  blue->color = Eigen::Vector3f(0.2f, 0.3f, 0.67f);
  blue->density = 0.02;
  blue->staining_power = 3.1;
  blue->granularity = 0.91;
  _pigments.push_back(blue);

  _s = Eigen::ArrayXXf::Zero(x_res, y_res);
  _c = Eigen::ArrayXXf::Zero(x_res, y_res);
  _c_max = 0.5f;
  _c_min = 0.1f;
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
  /* std::cout << "=====STEP START=====" << std::endl; */
  /* std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl; */
  /* std::cout << "d:" << _pigments[0]->d.maxCoeff() << " "; */
  /* std::cout << "dsum:" << _pigments[0]->d.sum() << " | "; */
  /* std::cout << "g:" << _pigments[0]->g.maxCoeff() << " "; */
  /* std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl; */
  /* std::cout << "-----updateVelocities" << std::endl; */
  updateVelocities();
  /* std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl; */
  /* std::cout << "d:" << _pigments[0]->d.maxCoeff() << " "; */
  /* std::cout << "dsum:" << _pigments[0]->d.sum() << " | "; */
  /* std::cout << "g:" << _pigments[0]->g.maxCoeff() << " "; */
  /* std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl; */
  /* std::cout << "-----relaxDivergence" << std::endl; */
  relaxDivergence();
  /* std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl; */
  /* std::cout << "d:" << _pigments[0]->d.maxCoeff() << " "; */
  /* std::cout << "dsum:" << _pigments[0]->d.sum() << " | "; */
  /* std::cout << "g:" << _pigments[0]->g.maxCoeff() << " "; */
  /* std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl; */
  /* std::cout << "-----flowOutward" << std::endl; */
  flowOutward();
  /* std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl; */
  /* std::cout << "d:" << _pigments[0]->d.maxCoeff() << " "; */
  /* std::cout << "dsum:" << _pigments[0]->d.sum() << " | "; */
  /* std::cout << "g:" << _pigments[0]->g.maxCoeff() << " "; */
  /* std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl; */
  /* std::cout << "-----movePigment" << std::endl; */
  movePigment();
  /* std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl; */
  /* std::cout << "d:" << _pigments[0]->d.maxCoeff() << " "; */
  /* std::cout << "dsum:" << _pigments[0]->d.sum() << " | "; */
  /* std::cout << "g:" << _pigments[0]->g.maxCoeff() << " "; */
  /* std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl; */
  /* std::cout << "-----transferPigment" << std::endl; */
  transferPigment();
  /* std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl; */
  /* std::cout << "d:" << _pigments[0]->d.maxCoeff() << " "; */
  /* std::cout << "dsum:" << _pigments[0]->d.sum() << " | "; */
  /* std::cout << "g:" << _pigments[0]->g.maxCoeff() << " "; */
  /* std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl; */
  /* std::cout << "-----simulateCapillaryFlow" << std::endl; */
  /* simulateCapillaryFlow(); */
  /* std::cout << "p:" << _pressure.maxCoeff() << " u:" << _u.absmax() << " v:" << _v.absmax() << std::endl; */
  /* std::cout << "d:" << _pigments[0]->d.maxCoeff() << " "; */
  /* std::cout << "dsum:" << _pigments[0]->d.sum() << " | "; */
  /* std::cout << "g:" << _pigments[0]->g.maxCoeff() << " "; */
  /* std::cout << "gsum: " << _pigments[0]->g.sum() << std::endl; */
  /* std::cout << "=====" << std::endl; */
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
  /* _dt = 0.10f / std::max(_u.absmax(), _v.absmax()); */
  _dt = std::min(0.01f, 0.1f / ceil(std::max(_u.absmax(), _v.absmax())));
  /* std::cout << "  _dt " << _dt << std::endl; */
  for (float t = 0.0f; t < 0.1f; t += _dt)
  {
    /* StaggeredGrid u_new(_x_res, _y_res, 0); */
    /* StaggeredGrid v_new(_x_res, _y_res, 1); */
    StaggeredGrid u_new(_u);
    StaggeredGrid v_new(_v);
    for (int i = 1; i < _x_res-1; i++)
      for (int j = 1; j < _y_res-1; j++)
      {
        /* std::cout << "b "; */
        a = _u.get(i,j)*_u.get(i,j) - _u.get(i+1,j)*_u.get(i+1,j) +
          _u.get(i+0.5f,j-0.5f)*_v.get(i+0.5f,j-0.5f) -
          _u.get(i+0.5f,j+0.5f)*_v.get(i+0.5f,j+0.5f);
        b = _u.get(i+1.5f,j) + _u.get(i-0.5f,j) + _u.get(i+0.5f,j+1) + _u.get(i+0.5f,j-1) - 4.0f*_u.get(i+0.5f,j);
        /* std::cout << b << std::endl; */
        /* u_new(i+0.5f,j) = _u.get(i+0.5f,j) + _dt * ( */
        u_new(i+0.5f,j) += _dt * (
            a - _viscosity*b + _pressure(i,j) - _pressure(i+1,j) - _viscous_drag*_u.get(i+0.5f,j));
        a = _v.get(i,j)*_v.get(i,j) - _v.get(i,j+1)*_v.get(i,j+1) +
          _u.get(i-0.5f,j+0.5f)*_v.get(i-0.5f,j+0.5f) -
          _u.get(i+0.5f,j+0.5f)*_v.get(i+0.5f,j+0.5f);
        b = _v.get(i+1,j+0.5f) + _v.get(i-1,j+0.5f) + _v.get(i,j+1.5f) + _v.get(i,j-0.5f) - 4.0f*_v.get(i,j+0.5f);
        /* std::cout << b << std::endl; */
        /* v_new(i,j+0.5f) = _v.get(i,j+0.5f) + _dt * ( */
        v_new(i,j+0.5f) += _dt * (
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
  for (int i = 0; i < _x_res; i++)
    for (int j = 0; j < _y_res; j++)
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

  StaggeredGrid u_new(_u);
  StaggeredGrid v_new(_v);
  while (1)
  {
    delta_max = 0.0f;
    u_new = _u;
    v_new = _v;
    for (int i = 1; i < _x_res-1; i++)
      for (int j = 1; j < _y_res-1; j++)
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
  const float eta = 0.05f;
  const int gaussian_radius = 4; // K = 9 ish

  Eigen::ArrayXXf M_copy = _M;
  Eigen::ArrayXXf M_blur(_x_res, _y_res);
  approximateGaussianBlur(M_copy, M_blur, _x_res, _y_res, gaussian_radius);
  for (int i = 1; i < _x_res-1; i++)
    for (int j = 1; j < _y_res-1; j++)
    {
      const float blur = eta * (1.0f - M_blur(i, j)) * _M(i,j);
      /* _pressure(i,j) -= std::max(std::min(blur, _pressure(i,j)), 0.0f); */
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
  /* _dt = 1.0f / std::max(_u.absmax(), _v.absmax()); */
  _dt = std::min(0.01f, 0.1f / ceil(std::max(_u.absmax(), _v.absmax())));
  /* std::cout << "movePigment _dt = " << _dt << std::endl; */
  // when more than one pigment, loop through each pigment
  for (Pigment* pig: _pigments)
  {
    for (float t = 0.0f; t < 0.1f; t += _dt)
    {
      Eigen::ArrayXXf g_new = pig->g;
      for (int i = 1; i < _x_res-1; i++)
        for (int j = 1; j < _y_res-1; j++)
        {
          const float gij = pig->g(i,j);
          g_new(i+1,j) += _dt*std::max(0.0f, _u.get(i+0.5f,j) * gij);
          g_new(i-1,j) += _dt*std::max(0.0f, -_u.get(i-0.5f,j) * gij);
          g_new(i,j+1) += _dt*std::max(0.0f, _v.get(i,j+0.5f) * gij);
          g_new(i,j-1) += _dt*std::max(0.0f, -_v.get(i,j-0.5f) * gij);
          g_new(i,j) -= _dt*(
              std::max(0.0f, _u.get(i+0.5f,j) * gij) +
              std::max(0.0f, -_u.get(i-0.5f,j) * gij) +
              std::max(0.0f, _v.get(i,j+0.5f) * gij) +
              std::max(0.0f, -_v.get(i,j-0.5f) * gij));
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
    for (int i = 1; i < _x_res-1; i++)
      for (int j = 1; j < _y_res-1; j++)
      {
        if (_M(i,j) != 1.0f)
          continue;

        /* delta_down = pig->g(i,j) * (1.0f - _h(i,j) * pig->granularity) * pig->density; */
        /* delta_up = pig->d(i,j) * (1.0f + (_h(i,j) - 1.0f) * pig->granularity) * pig->density / pig->staining_power; */
        delta_down = pig->d(i,j) * (1.0f + (_h(i,j) - 1.0f) * 0.2f) * 0.1f;
        delta_up = pig->g(i,j) * (1.0f - _h(i,j) * 0.2f) * 0.05f;

        /* delta_down = abs(delta_down); */
        /* delta_up = abs(delta_up); */

        if ((pig->d(i,j) + delta_down) > 1.0f)
          delta_down = std::max(0.0f, 1.0f - pig->g(i,j));
        if ((pig->g(i,j) + delta_up) > 1.0f)
          delta_up = std::max(0.0f, 1.0f - pig->d(i,j));
        /* std::cout <<"("<<i<<","<<j<<"): "<< delta_down << "," << delta_up << std::endl; */
        const float delta = delta_down - delta_up;
        /* pig->d(i,j) += delta_down - delta_up; */
        /* pig->g(i,j) += delta_up - delta_down; */
        pig->d(i,j) -= delta;
        pig->g(i,j) += delta;
        /* if (pig->g.sum() > 20) */
        /*   exit(0); */
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
  const float sigma = 0.1;

  for (int i = 0; i < _x_res; i++)
    for (int j = 0; j < _y_res; j++)
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

  std::cout << "s max " << _s.maxCoeff() << std::endl;

  // expand wet-area mask if cell's saturation exceeds a threshold sigma
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
      if (_s(i,j) > sigma)
        _M(i,j) = 1.0f;
}

void Watercolor2D::render()
{
  for (int i = 0; i < _x_res; i++)
    for (int j = 0; j < _y_res; j++)
    {
      Eigen::Vector3f paper = _paper[i][j];
      for (Pigment* pig: _pigments)
      {
        float opacity =
          pig->d(i,j) + pig->g(i,j);
        if (opacity < 0.0)
          opacity = 0.0;
        if (opacity > 1.0)
          opacity = 1.0;

        Eigen::Vector3f newcol =
          (1.0-opacity)*paper + opacity*pig->color;
        _buffer[j][i] = newcol;
        paper = newcol;

        // KS to reflectance
        /* Vector3f reflectance; */
        /* for (int x = 0; x < 3; x++) */
        /* { */
        /*   if (pig->K[i] <= 0.0f) */
        /*     reflectance[i] = 1.0f; */
        /*   else if (pig->S[i] <= 0) */
        /*     reflectance[i] = 0.0f; */
        /*   else */
        /*   { */
        /*     const float val = pig->S[i] / pig->K[i]; */
        /*     reflectance[i] = ( */
        /*         1.0 + val - sqrtf(2.0*val+1.0))/val; */
        /*   } */
        /* } */

        /* // reflectance to rgb */
        /* Vector3f rgbvec; */
        /* float sum = 0; */
        /* for (int i = 0; i < 3; i++) */
        /*   sum += reflectance[i]; */
        /* if (sum <= 0) { */
        /*     rgbvec[0] = rgbvec[1] = rgbvec[2] = 0; */
        /* } */
        /* else if (sum >= 3) { */
        /*     rgbvec[0] = rgbvec[1] = rgbvec[2] = 1; */
        /* } */
        /* else */
        /* { */
        /*   applyMatrix(3, m_wl, m_T, m_refvec, rgbvec); */
        /*   for (int i = 0; i < 3; i++) { */
        /*       if (rgbvec[i] < 0) rgbvec[i] = 0; */
        /*       if (rgbvec[i] > 1) rgbvec[i] = 1; */
        /*   } */
        /* } */


      }
    }
}
