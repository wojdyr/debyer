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
    StdDev() : n(0), mean(0.), S(0.) {}

    void add_x(double x)
    {
        ++n;
        double delta = x - mean;
        mean += delta / n;
        S += delta * (x - mean);
    }

    std::string str() const { return static_cast<std::ostringstream&>(
              std::ostringstream() << mean << " +- " << get_stddev()).str(); }

    double get_variance() const { return S / (n - 1); }
    double get_stddev() const { return sqrt(get_variance()); }
    double get_mean() const { return mean; }
    double get_n() const { return n; }

private:
    int n;
    double mean, S;
};

inline int iround(double d) { return static_cast<int>(floor(d+0.5)); }

#endif // DEBYER_UTILS_H_
