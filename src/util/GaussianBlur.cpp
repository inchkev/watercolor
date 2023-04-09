#include "GaussianBlur.h"
#include <cmath>
#include <vector>

using std::vector;

// lookup table of ideal averaging filter widths for n=3. replaces
// w_ideal = std::sqrt(12.0*sigma*sigma/(float)n + 1.0);
static const float W_IDEAL_N_EQUALS_3[10] = {
  1.0f,
  2.23606797749979f,
  4.123105625617661f,
  6.082762530298219f,
  8.06225774829855f,
  10.04987562112089f,
  12.041594578792296f,
  14.035668847618199f,
  16.0312195418814f,
  18.027756377319946f,
};

vector<int> boxesForGauss(int sigma, int n)
{
  // ideal averaging filter width
  const float w_ideal = W_IDEAL_N_EQUALS_3[sigma]; // std::sqrt(12.0*sigma*sigma/(float)n + 1.0);
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
  const float iarr = 1.0f / (r+r+1);
  for (int j = 0; j < h; j++)
  {
    const float fv = scl(0, j);
    const float lv = scl(w-1, j);
    int ti = 0; // pointer to target index
    int li = 0; // pointer to left index
    int ri = r; // pointer to right index
    float val = (r + 1) * fv;
    int i;
    for (i = 0; i < r; i++)
      val += scl(i, j);
    for (i = 0; i <= r; i++)
    {
      val += scl(ri++, j) - fv;           tcl(ti++, j) = val * iarr;
    }
    for (i = r + 1; i < w - r; i++)
    {
      val += scl(ri++, j) - scl(li++, j); tcl(ti++, j) = val * iarr;
    }
    for (i = w - r; i < w; i++)
    {
      val += lv           - scl(li++, j); tcl(ti++, j) = val * iarr;
    }
  }
}

void boxBlurT_4(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r)
{
  const float iarr = 1.0f / (r+r+1);
  for (int i = 0; i < w; i++)
  {
    const float fv = scl(i, 0);
    const float lv = scl(i, h-1);
    int ti = 0; // pointer to target index
    int li = 0; // pointer to left index
    int ri = r; // pointer to right index
    float val = (r + 1) * fv;
    int j;
    for (j = 0; j < r; j++)
      val += scl(i, j);
    for (j = 0; j <= r; j++)
    {
      val += scl(i, ri++) - fv;           tcl(i, ti++) = val * iarr;
    }
    for (j = r + 1; j < h - r; j++)
    {
      val += scl(i, ri++) - scl(i, li++); tcl(i, ti++) = val * iarr;
    }
    for (j = h - r; j < h; j++)
    {
      val += lv           - scl(i, li++); tcl(i, ti++) = val * iarr;
    }
  }
}

void boxBlur_4(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r)
{
  tcl = scl;
  boxBlurH_4(tcl, scl, w, h, r);
  boxBlurT_4(scl, tcl, w, h, r);
}

// approximate gaussian blur by applying box blur 3 times
void gaussBlur_4(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r)
{
  vector<int> boxes = boxesForGauss(r, 3);
  boxBlur_4(scl, tcl, w, h, (boxes[0] - 1) / 2);
  boxBlur_4(tcl, scl, w, h, (boxes[1] - 1) / 2);
  boxBlur_4(scl, tcl, w, h, (boxes[2] - 1) / 2);
}

void approximateGaussianBlur(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r)
{
  gaussBlur_4(scl, tcl, w, h, r);
}
