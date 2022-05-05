#ifndef READ_PPM_H
#define READ_PPM_H

#include <string>

using std::string;

void readPPM(const string& filename, int& xRes, int& yRes, float*& values);

void writePPM(const string& filename, int& xRes, int& yRes, const float*& values);

#endif
