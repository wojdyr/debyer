// cell method used for neighbour searching in dbr_extend
// Works only with orthorhombic PBC.

// Copyright (C) 2009 Marcin Wojdyr
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// $Id: $


#include "cells.h"
#include "iloops.h"
#include <algorithm>
#include <math.h>

using namespace std;


void CellMethod::init(double r)
{
    delete [] data_;

    for (int i = 0; i < 3; ++i) {
        double p = get_real_pbc(aconf_, i);
        n_[i] = (int) (p / r);
        // n_[i] must be > 0, 
        // for n_[i]==2 searching algorithm as implemented here doesn't work
        if (n_[i] < 3)
            n_[i] = 1;
        pbc_[i] = get_pbc(aconf_, i);
    }
    int number_of_cells = n_[0] * n_[1] * n_[2];

    cell_capacity_ = 1 + (int) ceil(4. * aconf_.n / number_of_cells);
    data_size_ = (number_of_cells + 1) * cell_capacity_;
    data_ = new int[data_size_];
    for (int i = 0; i != data_size_; ++i)
        data_[i] = -1;
    //if (verbosity > 1)
        printf("%d x %d x %d = %d cells created, cell capacity: %d\n",
                n_[0], n_[1], n_[2], number_of_cells, cell_capacity_);

    fill();
}

void CellMethod::fill()
{
    bool normalized = true;
    for (int i = 0; i < aconf_.n; ++i) {
        dbr_real* x = aconf_.atoms[i].xyz;
        if (x[1] < 0 || x[1] > pbc_[1]
                || x[1] < 0 || x[1] > pbc_[1]
                || x[2] < 0 || x[2] > pbc_[2])
            normalized = false;
        //assert (x[0] >= 0 && x[1] >= 0 && x[2] >= 0);
        int *p = find_cell_for_atom(x);
        while (*p != -1)
            ++p;
        *p = i;
    }
    if (!normalized)
        printf("Warning: atom coordinates are not wrapped to PBC box,\n"
               "         calculation of distances may not work.\n");
}

double CellMethod::get_dist(dbr_real const* a, dbr_real const* b) const
{
    double sum = 0.;
    for (int i = 0; i < 3; ++i) {
        double r = fabs(a[i] - b[i]);
        if (r > pbc_[i]/2.)
            r = pbc_[i] - r;
        sum += r * r;
    }
    return sqrt(sum);
}

// get vector v=b-a (in PBC) such that |v| is minimal
void CellMethod::get_shortest_vec(dbr_real const* a, dbr_real const* b,
                                  dbr_real *r) const
{
    for (int i = 0; i < 3; ++i) {
        r[i] = b[i] - a[i];
        if (r[i] >= pbc_[i] / 2.)
            r[i] -= pbc_[i];
        else if (r[i] < -pbc_[i] / 2.)
            r[i] += pbc_[i];
    }
}

dbr_real CellMethod::get_angle(const dbr_real *xyz1,
                               const dbr_real *xyz2,
                               const dbr_real *xyz3)
{
    dbr_real a[3], b[3];
    get_shortest_vec(xyz2, xyz1, a);
    get_shortest_vec(xyz2, xyz3, b);
    double t = dbr_dot3(a,b) / sqrt(dbr_dot3(a,a) * dbr_dot3(b,b));
    // it happens (very rarely) that due to rounding errors |t| > 1
    return fabs(t) < 1 ? acos(t) : 0.;
}

double CellMethod::min_cell_size() const
{
    return min(min(cell_size(0), cell_size(1)), cell_size(2));
}

// gets cell and atom coordinates; 
// if the cell indices are correct returns pointer to data of this cell,
// otherwise wraps the cell indices and shifts the atom.
int const* CellMethod::adjust_atom(int idx0, int idx1, int idx2,
                                   dbr_real* xyz) const
{
    int a[3] = { idx0, idx1, idx2 };
    for (int i = 0; i < 3; ++i) {
        int t = mod(a[i], n_[i]);
        if (t != a[i]) {
            xyz[i] += pbc_[i] * ((t - a[i]) / n_[i]);
            a[i] = t;
        }
    }
    return cell_start(a);
}

// search for atom in distance < epsilon from point xyz,
// returns pointer to the first atom found or NULL if not found.
dbr_atom const* CellMethod::find_one_in_cell(int const* idx,
                                    dbr_real const* xyz, double epsilon) const
{
    dbr_real x[3] = { xyz[0], xyz[1], xyz[2] };
    int const* p = adjust_atom(idx[0], idx[1], idx[2], x);
    while (*p != -1) {
        if (get_sq_dist(x, aconf_.atoms[*p].xyz) < epsilon * epsilon)
            return &aconf_.atoms[*p];
        ++p;
    }
    return NULL;
}

int CellMethod::find_nearest(dbr_real const* xyz) const
{
    int closest_atom = -1;
    double min_sq_dist = 1e9;
    int idx[3];
    find_cell_indices(xyz, idx);
    for (int i = idx[0]-1; i <= idx[0]+1; ++i)
        for (int j = idx[1]-1; j <= idx[1]+1; ++j)
            for (int k = idx[2]-1; k <= idx[2]+1; ++k) {
                dbr_real x[3] = { xyz[0], xyz[1], xyz[2] };
                int const* p = adjust_atom(i, j, k, x);
                while (*p != -1) {
                    double sq_dist = get_sq_dist(x, aconf_.atoms[*p].xyz);
                    if (sq_dist < min_sq_dist) {
                        closest_atom = *p;
                        min_sq_dist = sq_dist;
                    }
                    ++p;
                }
            }
    return closest_atom;
}

int CellMethod::find_nearest_to_atom(int atom_idx, double min_sq_dist) const
{
    // very similar to find_nearest()
    dbr_real const* xyz = aconf_.atoms[atom_idx].xyz; // this line was added
    int closest_atom = -1;
    int idx[3];
    find_cell_indices(xyz, idx);
    for (int i = idx[0]-1; i <= idx[0]+1; ++i)
        for (int j = idx[1]-1; j <= idx[1]+1; ++j)
            for (int k = idx[2]-1; k <= idx[2]+1; ++k) {
                dbr_real x[3] = { xyz[0], xyz[1], xyz[2] };
                int const* p = adjust_atom(i, j, k, x);
                while (*p != -1) {
                    // this condition was added
                    if (*p != atom_idx && !is_null(aconf_.atoms[*p])) {
                        double sq_dist = get_sq_dist(x, aconf_.atoms[*p].xyz);
                        if (sq_dist < min_sq_dist) {
                            closest_atom = *p;
                            min_sq_dist = sq_dist;
                        }
                    }
                    ++p;
                }
            }
    return closest_atom;
}

void CellMethod::find_all_bonds(double min_dist, vector<Bond> &bonds) const
{
    bonds.clear();
    for (int a = 0; a != aconf_.n; ++a) {
        const dbr_real* xyz = aconf_.atoms[a].xyz; // this line was added
        int idx[3];
        find_cell_indices(xyz, idx);
        for (int i = idx[0]-1; i <= idx[0]+1; ++i)
            for (int j = idx[1]-1; j <= idx[1]+1; ++j)
                for (int k = idx[2]-1; k <= idx[2]+1; ++k) {
                    dbr_real x[3] = { xyz[0], xyz[1], xyz[2] };
                    const int* p = adjust_atom(i, j, k, x);
                    while (*p != -1) {
                        if (*p > a) {
                            double d2 = get_sq_dist(x, aconf_.atoms[*p].xyz);
                            if (d2 < min_dist * min_dist) {
                                Bond bond = { a, *p, sqrt(d2) };
                                bonds.push_back(bond);
                            }
                        }
                        ++p;
                    }
                }
    }
}


double CellMethod::find_min_distance() const
{
    double min_dist = 1e9;
    for (int i = 0; i < aconf_.n; ++i) {
        int k = find_nearest_to_atom(i, min_dist*min_dist);
        if (k != -1) {
            min_dist = get_dist(aconf_.atoms[i].xyz, aconf_.atoms[k].xyz);
        }
    }
    return min_dist;
}

dbr_atom const*
CellMethod::get_atom_at(dbr_real const* xyz, double epsilon) const
{
    // pre: epsilon < r/2, where r=pbc_[i]/n_[i]
    int idx[3];
    find_cell_indices(xyz, idx);

    //int *p = cell_start(idx);
    //while (*p != -1) {
    //    if (get_sq_dist(xyz, aconf_.atoms[*p].xyz) < epsilon * epsilon)
    //        return &aconf_.atoms[*p];
    //    ++p;
    //}
    dbr_atom const* ret = find_one_in_cell(idx, xyz, epsilon);
    if (ret != NULL)
        return ret;

    int nb[3] = { 0, 0, 0 };
    for (int i = 0; i < 3; ++i) {
        double x0 = cx0(idx[i], i);
        if (xyz[i] - epsilon < x0)
            nb[i] = -1;
        else if (xyz[i] + epsilon > x0 + cell_size(i))
            nb[i] = 1;
        else
            continue;
        // side neighbours
        idx[i] += nb[i];
        dbr_atom const* ret = find_one_in_cell(idx, xyz, epsilon);
        if (ret != NULL)
            return ret;
        idx[i] -= nb[i];
    }
    if (abs(nb[0]) + abs(nb[1]) + abs(nb[2]) < 2)
        return NULL;

    for (int i = 0; i < 3; ++i)
        idx[i] += nb[i];
    // edge neighbours
    if (nb[0] != 0 && nb[1] != 0) {
        idx[2] -= nb[2];
        dbr_atom const* ret = find_one_in_cell(idx, xyz, epsilon);
        if (ret != NULL)
            return ret;
        idx[2] += nb[2];
    }
    if (nb[0] != 0 && nb[2] != 0) {
        idx[1] -= nb[1];
        dbr_atom const* ret = find_one_in_cell(idx, xyz, epsilon);
        if (ret != NULL)
            return ret;
        idx[1] += nb[1];
    }
    if (nb[1] != 0 && nb[2] != 0) {
        idx[0] -= nb[0];
        dbr_atom const* ret = find_one_in_cell(idx, xyz, epsilon);
        if (ret != NULL)
            return ret;
        idx[0] += nb[0];
    }

    // corner neighbours
    if (nb[0] != 0 && nb[1] != 0 && nb[2] != 0)
        return find_one_in_cell(idx, xyz, epsilon);

    return NULL;
}



dbr_atom const*
CellMethod::get_image(dbr_real const* r, int atom_nr, double epsilon,
                      dbr_real* delta) const
{
    dbr_atom const& atom = aconf_.atoms[atom_nr];
    dbr_real x[3] = { atom.xyz[0]+r[0], atom.xyz[1]+r[1], atom.xyz[2]+r[2] };
    //wrap_one(aconf, x);
    dbr_atom const* img = get_atom_at(x, epsilon);
    if (img == NULL || strcmp(img->name, atom.name) != 0)
        return NULL;
    else {
        if (delta != NULL)
            get_shortest_vec(x, img->xyz, delta);
        return img;
    }
}

