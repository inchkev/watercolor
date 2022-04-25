#ifndef GAUSSIAN_BLUR_H_
#define GAUSSIAN_BLUR_H_

#include <iostream>
#include <Eigen/Dense>

using Eigen::ArrayXXf;

/**
 * Approximate gaussian blur, from
 * http://blog.ivank.net/fastest-gaussian-blur.html.
 */
void approximateGaussianBlur(ArrayXXf& scl, ArrayXXf& tcl, int w, int h, int r);

#endif
