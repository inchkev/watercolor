#ifndef STAGGERED_GRID_H_
#define STAGGERED_GRID_H_

#include <Eigen/Dense>

class StaggeredGrid
{
public:
  StaggeredGrid(const int& rows, const int& cols);
  ~StaggeredGrid();

  inline float& operator()(float x, float y);
  const float operator()(float x, float y) const;
  inline float& operator()(int x, int y)          { return (*this)((float)x, (float)y); };
  const float operator()(int x, int y) const      { return (*this)((float)x, (float)y); };
  inline float& operator()(int x, float y)        { return (*this)((float)x, y); };
  const float operator()(int x, float y) const    { return (*this)((float)x, y); };
  inline float& operator()(float x, int y)        { return (*this)(x, (float)y); };
  const float operator()(float x, int y) const    { return (*this)(x, (float)y); };

  Eigen::ArrayXXf& data() { return _data; };
  int x_res() const { return _x_res; };
  int y_res() const { return _y_res; };
  int axis() const { return _axis; };
private:
  int _x_res;
  int _y_res;
  bool _axis; // false = x, true = y

  Eigen::ArrayXXf _data;
};

#endif
