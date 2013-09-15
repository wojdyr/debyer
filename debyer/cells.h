// cell method used for neighbour searching in dbr_extend
// Works only with orthorhombic PBC.

// Copyright 2009 Marcin Wojdyr
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

#ifndef DEBYER_CELLS_H_
#define DEBYER_CELLS_H_

#include <vector>
#include "fileio.h"

inline int mod(int a, int n)
{
    // pre: n > 0
    int r = a % n;
    return r >= 0 ? r : r + n;
}

struct Bond
{
    int a, b; // indices of bonded atoms
    float length;
};


// supports only orthorhombic PBC
class CellMethod
{
public:
    CellMethod(dbr_aconf const& aconf, double r)
        : aconf_(aconf), data_(NULL) { init(r); }
    ~CellMethod() { delete [] data_; }
    void init(double r);

    // Gets one atom (in unspecified order) in distance r < epsilon from given
    // point xyz. Returns NULL if there are no such atoms.
    dbr_atom const* get_atom_at(dbr_real const* xyz, double epsilon) const;

    int find_nearest(dbr_real const* xyz) const;
    int find_nearest_to_atom(int atom_idx, double min_sq_dist) const;
    void find_all_bonds(double min_dist, std::vector<Bond> &bonds) const;

    double get_dist(dbr_real const* a, dbr_real const* b) const;
    void get_shortest_vec(dbr_real const* a, dbr_real const* b,
                          dbr_real *r) const;
    dbr_real get_angle(const dbr_real *xyz1, const dbr_real *xyz2,
                                                       const dbr_real *xyz3);

    // returns min(cx, cy, cz), where cx is the size of one cell in x, etc.
    double min_cell_size() const;

    double const* pbc() const { return pbc_; }
    dbr_aconf const& aconf() const { return aconf_; }

    double find_min_distance() const;

    // returns atom of the same kind as atom atom_nr located at x+r,
    // where x is the position of atom atom_nr
    dbr_atom const* get_image(dbr_real const* r, int atom_nr,
                              double epsilon, dbr_real* delta=NULL) const;

private:
    void fill();

    int* cell_start(int a0, int a1, int a2) const
        { return data_ +  cell_capacity_ * (a0*n_[1]*n_[2] + a1*n_[2] + a2); }

    int* cell_start(int *a) const { return cell_start(a[0], a[1], a[2]); }

    int* cell_start_safe(int a0, int a1, int a2) const
        { return cell_start(mod(a0, n_[0]), mod(a1, n_[1]), mod(a2, n_[2])); }

    int cpos(dbr_real const* xyz, int d) const
        { return xyz[d] / pbc_[d] * n_[d]; }

    double cell_size(int dim) const { return pbc_[dim] / n_[dim]; }

    double real_cell_size(int dim) const
        { return get_real_pbc(aconf_, dim) / n_[dim]; }

    double cx0(int cell_idx, int d) const { return cell_idx * cell_size(d); }

    int* find_cell_for_atom(const dbr_real *xyz) const
        { return cell_start_safe(cpos(xyz, 0), cpos(xyz, 1), cpos(xyz, 2)); }

    void find_cell_indices(const dbr_real *xyz, int* cell_index) const {
        for (int i = 0; i < 3; ++i)
            cell_index[i] = (int) (xyz[i] / cell_size(i));
    }

    int const* adjust_atom(int idx0, int idx1, int idx2, dbr_real* xyz) const;

    dbr_atom const* find_one_in_cell(int const* idx,
                                    dbr_real const* xyz, double epsilon) const;

    dbr_aconf const& aconf_;
    double pbc_[3];
    int cell_capacity_;
    int data_size_;
    int n_[3];
    int *data_;
};

#endif // DEBYER_CELLS_H_
