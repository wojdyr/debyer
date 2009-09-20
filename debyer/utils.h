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
 *  $Id: utils.h 107 2009-04-12 22:39:39Z wojdyr $
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
    double n() const { return n_; }

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

inline void rodrigues(dbr_real theta, dbr_real axis[3], dbr_real mat[3][3])
{
    dbr_real alen = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    for (int i = 0; i < 3; ++i)
        axis[i] /= alen;

    dbr_real c = cos(theta);
    dbr_real s = sin(theta);

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

// r - result
inline void matrix_dot_vec3(const dbr_real mat[3][3], const dbr_real v[3],
                            dbr_real r[3])
{
    for (int i = 0; i < 3; ++i)
        r[i] = mat[i][0] * v[0] + mat[i][1] * v[1] + mat[i][2] * v[2];
}


#endif // DEBYER_UTILS_H_
