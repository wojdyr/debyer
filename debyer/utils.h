/*  debyer -- program for calculation of diffration patterns
 *  Copyright (C) 2006-2007 Marcin Wojdyr
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  $Id$
 *
 *  Misc utilities (in C++). ATM all functions/methods are inlined.
 */
#ifndef DEBYER_UTILS_H_
#define DEBYER_UTILS_H_

#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <cstdlib>


// popular single-pass algorithm for calculation of variance and mean
class StdDev
{
public:
    StdDev() : n_(0), mean_(0.), S_(0.) {}

    void add_x(double x)
    {
        ++n_;
        double delta = x - mean_;
        mean_ += delta / n_;
        S_ += delta * (x - mean_);
    }

    std::string str() const { return static_cast<std::ostringstream&>(
              std::ostringstream() << mean() << " +- " << stddev()).str(); }

    double variance() const { return S_ / (n_ - 1); }
    double stddev() const { return sqrt(variance()); }
    double mean() const { return mean_; }
    int n() const { return n_; }

private:
    int n_;
    double mean_, S_;
};

inline int iround(double d) { return static_cast<int>(floor(d+0.5)); }

/// S() converts to string
template <typename T>
inline std::string S(T k) {
    return static_cast<std::ostringstream&>(std::ostringstream() << k).str();
}

inline void rodrigues(double theta, double axis[3], double mat[3][3])
{
    double alen = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    for (int i = 0; i < 3; ++i)
        axis[i] /= alen;

    double c = cos(theta);
    double s = sin(theta);

    mat[0][0] = (1-c)*axis[0]*axis[0] + c;
    mat[0][1] = (1-c)*axis[0]*axis[1] - s*axis[2];
    mat[0][2] = (1-c)*axis[0]*axis[2] + s*axis[1];

    mat[1][0] = (1-c)*axis[0]*axis[1] + s*axis[2];
    mat[1][1] = (1-c)*axis[1]*axis[1] + c;
    mat[1][2] = (1-c)*axis[1]*axis[2] - s*axis[0];

    mat[2][0] = (1-c)*axis[2]*axis[0] - s*axis[1];
    mat[2][1] = (1-c)*axis[2]*axis[1] + s*axis[0];
    mat[2][2] = (1-c)*axis[2]*axis[2] + c;
}

inline
void inverse_3x3_matrix(const double a[][3], double b[][3])
{
    int i, j;
    b[0][0] =  (a[1][1] * a[2][2] - a[1][2] * a[2][1]);
    b[0][1] = -(a[0][1] * a[2][2] - a[0][2] * a[2][1]);
    b[0][2] =  (a[0][1] * a[1][2] - a[0][2] * a[1][1]);
    b[1][0] = -(a[1][0] * a[2][2] - a[1][2] * a[2][0]);
    b[1][1] =  (a[0][0] * a[2][2] - a[0][2] * a[2][0]);
    b[1][2] = -(a[0][0] * a[1][2] - a[0][2] * a[1][0]);
    b[2][0] =  (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    b[2][1] = -(a[0][0] * a[2][1] - a[0][1] * a[2][0]);
    b[2][2] =  (a[0][0] * a[1][1] - a[0][1] * a[1][0]);

    double s = 1. / (a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0]);

    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            b[i][j] *= s;
}


#endif // DEBYER_UTILS_H_
