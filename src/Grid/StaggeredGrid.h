#ifndef STAGGERED_GRID_H_
#define STAGGERED_GRID_H_

#include <Eigen/Dense>

class StaggeredGrid
{
public:
  StaggeredGrid(const int x_res, const int y_res, bool axis);
  StaggeredGrid(const StaggeredGrid& g);
  StaggeredGrid& operator=(const StaggeredGrid& g);
  StaggeredGrid& operator+=(const Eigen::ArrayXXf a);
  StaggeredGrid& operator-=(const Eigen::ArrayXXf a);
  ~StaggeredGrid() {};

  float& operator()(float x, float y);
  const float operator()(float x, float y) const;
  float& operator()(int x, int y)                 { return (*this)((float)x, (float)y); };
  const float operator()(int x, int y) const      { return (*this)((float)x, (float)y); };
  float& operator()(int x, float y)               { return (*this)((float)x, y); };
  const float operator()(int x, float y) const    { return (*this)((float)x, y); };
  float& operator()(float x, int y)               { return (*this)(x, (float)y); };
  const float operator()(float x, int y) const    { return (*this)(x, (float)y); };

  // read only access, interpolates at non-boundary locations
  const float get(float x, float y);
  const float get(int x, float y)    { return get((float)x, y); };
  const float get(float x, int y)    { return get(x, (float)y); };
  const float get(int x, int y)      { return get((float)x, (float)y); };

  Eigen::ArrayXXf& data() { return _data; };
  int x_res() const { return _x_res; };
  int y_res() const { return _y_res; };
  int is_y_axis() const { return _is_y_axis; };

  float absmax() const;
  float max() const;
  float min() const;
private:
  int _x_res;
  int _y_res;
  bool _is_y_axis; // false = x, true = y

  Eigen::ArrayXXf _data;
};

#endif
