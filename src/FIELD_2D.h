#ifndef FIELD_2D_H
#define FIELD_2D_H

#include <cstdlib>
#include <iostream>
#include <cassert>

using namespace std;

class FIELD_2D {

public:
  // construction and destruction
  FIELD_2D();
  FIELD_2D(int xRes, int yRes); 
  virtual ~FIELD_2D();

  // (x,y) accessor
  inline float& operator()(int x, int y) {
    if (x < 0)
      std::cout << x << "," << y << std::endl;
    assert(x >= 0);
    if (x >= _xRes)
      std::cout << x << "," << y << std::endl;
    assert(x < _xRes);
    if (y < 0)
      std::cout << x << "," << y << std::endl;
    assert(y >= 0);
    if (y >= _yRes)
      std::cout << x << "," << y << std::endl;
    assert(y < _yRes);
    return _data[x + y * _xRes];
  };

  // set one field to another
  FIELD_2D& operator=(const FIELD_2D& field);

  // scale the field by a constant
  FIELD_2D& operator*=(const float& scalar);

  // add another field to this field
  FIELD_2D& operator+=(const FIELD_2D& field);

  // subtract another field from this field
  FIELD_2D& operator-=(const FIELD_2D& field);

  // accessors
  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  int totalCells() const { return _xRes * _yRes; };
  float*& data() { return _data; };

  // set the size of the array
  void resize(int xRes, int yRes);

  // set to zero
  void clear();

  // squared sum of entries
  float squaredSum() const;
  
  // dot product with another field
  float dotProduct(FIELD_2D& field) const;

  // axpy operation with another field
  // this = this + scalar * field
  void axpy(const float& scalar, const FIELD_2D& field);

private:
  int _xRes;
  int _yRes;

  float* _data;
};

#endif
