#include <cmath>
#include "StaggeredGrid.h"

StaggeredGrid::StaggeredGrid(const int x_res, const int y_res, bool axis) :
  _x_res(x_res),
  _y_res(y_res),
  _axis(axis),
  _data(x_res + 1, y_res * 2 + 1)
{
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
          return _data(xi, _y_res * 2 - 1);
        assert(false);
      }
      if (xi == 0)
        return _data(0, yi);
      if (xi == _x_res)
        return _data(_x_res-1, yi);
      assert(false);
    }
    return _data(xi, yi);
  }
  if (!y_border)
    assert(false);
  return _data(xi, yi);
}

const float StaggeredGrid::operator()(float x, float y) const
{
  assert (x >= -0.5f && x <= cols + 0.5f);
  assert (y >= -0.5f && y <= rows + 0.5f);

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
