#include "GaussianBlur.h"
#include <cmath>
#include <vector>

using std::vector;

vector<int> &boxesForGauss(int sigma, int n)
{
  const float w_ideal = std::sqrt(12.0*sigma*sigma/(float)n + 1.0);
  int wl = std::floor(w_ideal);
  if (wl % 2 == 0)
    wl--;
  const int wu = wl + 2;

  const int m = std::round((float)(12*sigma*sigma - n*wl*wl - 4*n*wl - 3*n) / (-4*wl - 4));

  vector<int> sizes;
  for (int i = 0; i < n; i++)
    sizes.push_back(i < m ? wl : wu);
  return sizes;
}

void boxBlurH_4(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r)
{
  tcl = scl;
  boxBlurH_4(tcl, scl, w, h, r);
  boxBlurT_4(scl, tcl, w, h, r);
}

void boxBlurT_4(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r)
{
  const float iarr = 1.0 / (r+r+1);
}

void boxBlur_4(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r)
{
}

// approximate gaussian blur by applying box blur 3 times
void gaussBlur_4(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r)
{
  vector<int>& boxes = boxesForGauss(r, 3);
  boxBlur_4(scl, tcl, w, h, (boxes[0] - 1) / 2);
  boxBlur_4(tcl, scl, w, h, (boxes[1] - 1) / 2);
  boxBlur_4(scl, tcl, w, h, (boxes[2] - 1) / 2);
}

void approximateGaussianBlur(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r)
{
  gaussBlur_4(scl, tcl, w, h, r);
}
