#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
#include "StaggeredGrid.h"

StaggeredGrid::StaggeredGrid(const int x_res, const int y_res, bool axis) :
  _x_res(x_res),
  _y_res(y_res),
  _axis(axis)
{
  _data = Eigen::ArrayXXf::Zero(x_res + 1, y_res * 2 + 1);
}

StaggeredGrid::StaggeredGrid(const StaggeredGrid& g) :
  _x_res(g._x_res),
  _y_res(g._y_res),
  _axis(g._axis)
{
  _data = g._data;
}

StaggeredGrid& StaggeredGrid::operator=(const StaggeredGrid& g)
{
  assert(g.x_res() == _x_res);
  assert(g.y_res() == _y_res);
  assert(g.axis() == _axis);
  _data = g._data;
  return *this;
}

StaggeredGrid& StaggeredGrid::operator+=(const Eigen::ArrayXXf a)
{
  assert(a.rows() == _x_res);
  assert(a.cols() == _y_res);
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
    {
      const float val = a(i,j) * 0.5f;
      if (_axis)
      {
        (*this)(i-0.5f,j) -= val;
        (*this)(i+0.5f,j) -= val;
      }
      else
      {
        (*this)(i,j-0.5f) -= val;
        (*this)(i,j+0.5f) -= val;
      }
    }
  return *this;
}

StaggeredGrid& StaggeredGrid::operator-=(const Eigen::ArrayXXf a)
{
  assert(a.rows() == _x_res);
  assert(a.cols() == _y_res);
  for (int j = 0; j < _y_res; j++)
    for (int i = 0; i < _x_res; i++)
    {
      const float val = a(i,j) * 0.5f;
      if (_axis)
      {
        (*this)(i-0.5f,j) -= val;
        (*this)(i+0.5f,j) -= val;
      }
      else
      {
        (*this)(i,j-0.5f) -= val;
        (*this)(i,j+0.5f) -= val;
      }
    }
  return *this;
}

float& StaggeredGrid::operator()(float x, float y)
{
  assert (x >= -0.5f && x <= (float)_x_res + 0.5f);
  assert (y >= -0.5f && y <= (float)_y_res + 0.5f);

  int xi = (int) (x + 0.5f);
  int yi = (int) (y * 2.0f) + 1;
  bool x_border = (x - std::trunc(x) == .5 || x - std::trunc(x) == -.5);
  bool y_border = (y - std::trunc(y) == .5 || y - std::trunc(y) == -.5);

  // TODO: can i make this cleaner?
  if (x_border)
  {
    if (y_border)
    {
      if (_axis)
      {
        if (yi == 0)
          return _data(xi, 1);
        if (yi == _y_res * 2)
          return _data(xi, _y_res*2-1);
        throw std::invalid_argument("invalid argument");
      }
      if (xi == 0)
        return _data(0, yi);
      if (xi == _x_res)
        return _data(_x_res-1, yi);
      throw std::invalid_argument("invalid argument");
    }
    return _data(xi, yi);
  }
  if (!y_border)
    throw std::invalid_argument("invalid argument");
  return _data(xi, yi);
}

const float StaggeredGrid::operator()(float x, float y) const
{
  assert (x >= -0.5f && x <= (float)_x_res + 0.5f);
  assert (y >= -0.5f && y <= (float)_y_res + 0.5f);

  int xi = (int) (x + 0.5f);
  int yi = (int) (y * 2.0f) + 1;
  bool x_border = (x - std::trunc(x) == .5 || x - std::trunc(x) == -.5);
  bool y_border = (y - std::trunc(y) == .5 || y - std::trunc(y) == -.5);

  if (x_border)
  {
    if (y_border)
    {
      if (_axis)
      {
        if (yi == 0)
          return _data(xi, 1);
        else if (yi == _y_res * 2)
          return _data(xi, _y_res * 2 - 1);
        else
          return (_data(xi, yi-1) + _data(xi, yi+1)) * 0.5f;
      }
      else
      {
        if (xi == 0)
          return _data(0, yi);
        else if (xi == _x_res)
          return _data(_x_res-1, yi);
        else
          return (_data(xi-1, yi) + _data(xi, yi)) * 0.5f;
      }
    }
    else
      return _data(xi, yi);
  }
  else
  {
    if (y_border)
      return _data(xi, yi);
    else
    {
      if (_axis)
        return (_data(xi, yi-1) + _data(xi, yi+1)) * 0.5f;
      else
        return (_data(xi, yi) + _data(xi+1, yi)) * 0.5f;
    }
  }
}

const float StaggeredGrid::get(float x, float y)
{
  assert (x >= -0.5f && x <= (float)_x_res + 0.5f);
  assert (y >= -0.5f && y <= (float)_y_res + 0.5f);
  if (x < -0.5f || x > (float)_x_res + 0.5f || y < -0.5f || y > (float)_y_res + 0.5f)
  {
    std::cout << "INVALID " << x << ", " << y << std::endl;
    throw;
  }

  int xi = (int) (x + 0.5f);
  int yi = (int) (y * 2.0f) + 1;
  bool x_border = (x - std::trunc(x) == .5 || x - std::trunc(x) == -.5);
  bool y_border = (y - std::trunc(y) == .5 || y - std::trunc(y) == -.5);

  if (x_border)
  {
    if (y_border)
    {
      if (_axis)
      {
        if (yi == 0)
          return _data(xi, 1);
        if (yi == _y_res * 2)
          return _data(xi, _y_res * 2 - 1);
        return (_data(xi, yi-1) + _data(xi, yi+1)) * 0.5f;
      }
      if (xi == 0)
        return _data(0, yi);
      if (xi == _x_res)
        return _data(_x_res-1, yi);
      return (_data(xi-1, yi) + _data(xi, yi)) * 0.5f;
    }
    return _data(xi, yi);
  }
  if (y_border)
    return _data(xi, yi);
  if (_axis)
    return (_data(xi, yi-1) + _data(xi, yi+1)) * 0.5f;
  return (_data(xi, yi) + _data(xi+1, yi)) * 0.5f;
}

float StaggeredGrid::absmax() const
{
  float max_value = 0.0f;
  for (int j = 0; j < _y_res * 2 + 1; j++)
    for (int i = 0; i < _x_res + (j % 2); i++)
      max_value = std::max(max_value, std::abs(_data(i, j)));
  return max_value;
}

float StaggeredGrid::max() const
{
  float max_value = std::numeric_limits<float>::min();
  for (int j = 0; j < _y_res * 2 + 1; j++)
    for (int i = 0; i < _x_res + (j % 2); i++)
      max_value = std::max(max_value, _data(i, j));
  return max_value;
}

float StaggeredGrid::min() const
{
  float min_value = std::numeric_limits<float>::max();
  for (int j = 0; j < _y_res * 2 + 1; j++)
    for (int i = 0; i < _x_res + (j % 2); i++)
      min_value = std::min(min_value, _data(i, j));
  return min_value;
}
