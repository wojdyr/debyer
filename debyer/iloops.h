/*  debyer -- program for calculation of diffration patterns
 *  Copyright 2006 Marcin Wojdyr
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
 */

/* This file is to be included only by debyer.c
 * and contains inner loops of ID calculation, where most of
 * computation time is spend. It is put into separate file to make
 * optimization easier.
 */

#ifndef DEBYER_ILOOPS_H_
#define DEBYER_ILOOPS_H_

#include <math.h>

#ifdef __STRICT_ANSI__
# define SQRT sqrt
#else
# define SQRT sqrtf
#endif

dbr_real get_sq_dist(const dbr_real *xyz1, const dbr_real *xyz2)
{
    dbr_real dx = xyz1[0] - xyz2[0];
    dbr_real dy = xyz1[1] - xyz2[1];
    dbr_real dz = xyz1[2] - xyz2[2];
    return dx*dx + dy*dy + dz*dz;
}

void calculate_irdf_innerloop(int n, int nbins, dbr_real rquanta,
                              const dbr_real *xyz1, const dbr_cell* c2, int *t)
{
    int j, bin;
    for (j = 0; j != n; ++j) {
        bin = (int) (SQRT(get_sq_dist(xyz1, c2->atoms[j])) / rquanta);
        if (bin < nbins)
            ++t[bin];
        /*
        else
            fprintf(stderr, "Out: %i>=%i\n", bin, nbins);
        */
    }
}

void calculate_irdf_innerloop2(int n, int nbins, dbr_real rquanta,
                               const dbr_real *xyz1, const dbr_cell* c2, int *t,
                               dbr_real rcut2)
{
    int j, bin;
    for (j = 0; j != n; ++j) {
        dbr_real d2 = get_sq_dist(xyz1, c2->atoms[j]);
        if (d2 > rcut2)
            continue;
        bin = (int) (SQRT(d2) / rquanta);
        if (bin < nbins) /* TODO this check should not be neccessary */
            ++t[bin];
    }
}


#endif // DEBYER_ILOOPS_H_
