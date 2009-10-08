// dbr_extend -- tool to extend a system (configuration for atomistic
// simulation) in x, y or z direction by increasing PBC and filling the
// new space with periodic structure, a copy of the structure around it.
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
// $Id$

// TODO: 
// * it doesn't work when there are two cells in given direction,
//   ./dbr_extend -r  -e0.02 -c2.99 -M0.03  -v -v  -z0.5 mydisl.cfg.gz
// * option -cN with N>1 should change PBC only once

#include <algorithm>
#include <cstring>

#include "fileio.h"
#include "iloops.h"
#include "utils.h"
#include "extend_cmd.h"

using namespace std;

static int verbosity = 0;

const char* usage_examples[] = {
"Usage examples:",
"  dbr_extend -z41.5 -e0.2 -vvv file.cfg",
"    Try to find periodicity of the structure in z direction starting",
"    from z0=41.5. It tries to find z1 such that each atom",
"    with z0 < z < z1 has a periodic image with",
"    x'=x, y'=y, z'=z+delta, delta=z1-z0.",
"    Comparisons of coordinates are done with epsilon=0.2.",
"  dbr_extend -c10 -z41.5 -e0.2 -i file.cfg",
"    Extend PBC box by in z direction by 10 x the value reported",
"    in the previous try run, copy atoms to the created space,",
"    write the configuration back to file.cfg."
};


struct Slab
{
    int dim;
    double x0;
    double delta;
};

int mod(int a, int n)
{
    // pre: n > 0
    int r = a % n;
    return r >= 0 ? r : r + n;
}

dbr_real dist_forward(dbr_real x2, dbr_real x1, dbr_real pbc)
{
    return x2 >= x1 ? x2 - x1 : x2 - x1 + pbc;
}

double* get_pbc_ptr(dbr_aconf& aconf, int dim)
{
    switch (dim) {
        case 0: return &aconf.pbc.v00;
        case 1: return &aconf.pbc.v11;
        case 2: return &aconf.pbc.v22;
        default: return NULL;
    }
}

double get_real_pbc(dbr_aconf const& aconf, int dim)
{
    switch (dim) {
        case 0: return aconf.pbc.v00;
        case 1: return aconf.pbc.v11;
        case 2: return aconf.pbc.v22;
        default: return 0;
    }
}

void wrap_to_pbc(dbr_aconf& aconf)
{
    if (aconf.reduced_coordinates)
        for (int i = 0; i != aconf.n; ++i)
            for (int j = 0; j < 3; ++j)
                aconf.atoms[i].xyz[j] -= floor(aconf.atoms[i].xyz[j]);
    else {
        // works only for orthorhombic PBC
        double pbc[3] = { aconf.pbc.v00, aconf.pbc.v11, aconf.pbc.v22 };
        for (int i = 0; i != aconf.n; ++i)
            for (int j = 0; j < 3; ++j) {
                double n = floor(aconf.atoms[i].xyz[j] / pbc[j]);
                aconf.atoms[i].xyz[j] -= n * pbc[j];
            }
    }
}

void wrap_one(dbr_aconf const& aconf, dbr_real *xyz)
{
    if (aconf.reduced_coordinates)
        for (int j = 0; j < 3; ++j)
            xyz[j] -= floor(xyz[j]);
    else {
        // works only for orthorhombic PBC
        double pbc[3] = { aconf.pbc.v00, aconf.pbc.v11, aconf.pbc.v22 };
        for (int j = 0; j < 3; ++j)
            xyz[j] -= floor(xyz[j] / pbc[j]) * pbc[j];
    }
}

// expand PBC, but don't move atoms, i.e. change reduced coordinates
// in case of reduced coords, delta is also reduced
void expand_pbc(dbr_aconf& aconf, int dim, double delta)
{
    double* pp = get_pbc_ptr(aconf, dim);
    double old_p = *pp;
    double new_p = aconf.reduced_coordinates ? old_p * (1. + delta)
                                             : old_p + delta;

    if (aconf.reduced_coordinates)
        for (int i = 0; i != aconf.n; ++i) {
            aconf.atoms[i].xyz[dim] /= (1. + delta);
            if (aconf.atoms[i].xyz[dim] > 1)
                aconf.atoms[i].xyz[dim] -= floor(aconf.atoms[i].xyz[dim]);
        }

    *pp = new_p;

    wrap_to_pbc(aconf);

    if (verbosity > 0)
        printf("new PBC[%d]=%g (was: %g)\n", dim, new_p, old_p);
}

double get_pbc(dbr_aconf const& aconf, int dim)
{
    double p = get_real_pbc(aconf, dim);
    return aconf.reduced_coordinates && p != 0. ? 1. : p;
}

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

    double get_dist(dbr_real const* a, dbr_real const* b) const;
    void get_shortest_vec(dbr_real const* a, dbr_real const* b,
                          dbr_real *r) const;

    double half_cell() const
        { return min(min(cell_size(0), cell_size(1)), cell_size(2)); }

    double const* pbc() const { return pbc_; }
    dbr_aconf const& aconf() const { return aconf_; }

    double find_min_distance() const;

    // returns atom of the same kind as atom atom_nr located at x+r,
    // where x is the position of atom atom_nr
    dbr_atom const* get_image(dbr_real const* r, int atom_nr,
                              double epsilon, dbr_real* delta=NULL) const;

    vector<double> get_image_distances(int atom_nr,
                                       double lo_lim, double hi_lim,
                                       double epsilon) const;

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

void CellMethod::init(double r)
{
    delete [] data_;

    for (int i = 0; i < 3; ++i) {
        double p = get_real_pbc(aconf_, i);
        n_[i] = p > r ? (int) (p / r) : 1;
        pbc_[i] = get_pbc(aconf_, i);
    }
    int number_of_cells = n_[0] * n_[1] * n_[2];

    cell_capacity_ = 1 + (int) ceil(4. * aconf_.n / number_of_cells);
    data_size_ = (number_of_cells + 1) * cell_capacity_;
    data_ = new int[data_size_];
    for (int i = 0; i != data_size_; ++i)
        data_[i] = -1;
    if (verbosity > 1)
        printf("%d x %d x %d = %d cells created, cell capacity: %d\n",
                n_[0], n_[1], n_[2], number_of_cells, cell_capacity_);

    fill();
}

void CellMethod::fill()
{
    for (int i = 0; i < aconf_.n; ++i) {
        dbr_real* x = aconf_.atoms[i].xyz;
        //assert (x[0] >= 0 && x[1] >= 0 && x[2] >= 0);
        int *p = find_cell_for_atom(x);
        while (*p != -1)
            ++p;
        *p = i;
    }
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

vector<double> CellMethod::get_image_distances(int atom_nr,
                                               double lo_lim, double hi_lim,
                                               double epsilon) const
{
    vector<double> dd(aconf_.n, 0.);
    dbr_atom const& atom0 = aconf_.atoms[atom_nr];
    dbr_real const* xyz0 = atom0.xyz;
    for (int i = 0; i != aconf_.n; ++i) {
        if (i == atom_nr || strcmp(atom0.name, aconf_.atoms[i].name) != 0)
            continue;
        double dist = get_dist(xyz0, aconf_.atoms[i].xyz);
        if (dist < lo_lim || dist > hi_lim)
            continue;
        dbr_real r[3];
        dbr_diff3(aconf_.atoms[i].xyz, xyz0, r);
        bool ok = true;
        for (int j = 0; j != aconf_.n; ++j) {
            if (get_image(r, j, epsilon) == NULL) {
                ok = false;
                break;
            }
        }
        if (ok) {
            dd[i] = dist;
            //printf("debug: %d r=(%f %f %f)\n", i, r[0], r[1], r[2]);
        }
    }
    return dd;
}

void merge_atoms(dbr_aconf& aconf, gengetopt_args_info const& args)
{
    CellMethod cm(aconf, args.min_cell_arg);
    int new_n = aconf.n;
    for (int i = 0; i != aconf.n; ++i) {
        if (is_null(aconf.atoms[i]))
            continue;
        double min_sq_dist = args.epsilon_arg * args.epsilon_arg;
        int k = cm.find_nearest_to_atom(i, min_sq_dist);
        if (k != -1) {
            if (strcmp(aconf.atoms[i].name, aconf.atoms[k].name) != 0)
                printf("WARNING: atoms with different symbols were merged\n");
            // mark atoms that are to be deleted
            int t = min(i, k);
            nullify(aconf.atoms[t]);
            --new_n;
        }
        // TODO? : average coordinates of merged atoms
    }
    double ratio = double(aconf.n) / new_n;
    printf("merging, i.e. removing duplicates: %g%% (1/%g) of atoms left\n",
            100./ratio, ratio);

    int n = 0;
    for (int i = 0; i != aconf.n; ++i) {
        if (!is_null(aconf.atoms[i])) {
            dbr_copy_atom(aconf, i, n);
            if (!aconf.auxiliary.empty() && i != n)
                aconf.auxiliary[n] = aconf.auxiliary[i];
            ++n;
        }
    }
    assert (n == new_n);
    if (!aconf.auxiliary.empty())
        aconf.auxiliary.resize(new_n);

    string comment = "merging atoms, " + S(aconf.n) + " -> " + S(new_n);
    aconf.comments.insert(aconf.comments.begin(), comment);

    aconf.n = new_n;
}

void rotate_atoms(dbr_aconf& aconf, double angle, double axis[3])
{
    double rot[3][3];
    rodrigues(angle, axis, rot);
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* x = aconf.atoms[i].xyz;
        dbr_real r[3];
        matrix_dot_vec3(rot, x, r);
        for (int i = 0; i < 3; ++i)
            x[i] = r[i];
    }
}

void transform1(dbr_aconf& aconf, gengetopt_args_info const& args)
{
    wrap_to_pbc(aconf);
    aconf.pbc.v00 /= sqrt(2);
    aconf.pbc.v11 /= sqrt(2);
    double axis[3] = { 0, 0, 1 };
    rotate_atoms(aconf, M_PI / 4, axis);
    wrap_to_pbc(aconf);
    merge_atoms(aconf, args);
}

// select the last atom that in direction `dim' has position lesser than `x0'.
int select_atom(dbr_aconf const& aconf, int dim, double x0)
{
    int a0 = -1;
    double min_diff = 1e9;
    double pbcd = get_pbc(aconf, dim);
    for (int i = 0; i != aconf.n; ++i) {
        double x = aconf.atoms[i].xyz[dim];
        double diff = dist_forward(x0, x, pbcd);
        if (diff < min_diff) {
            a0 = i;
            min_diff = diff;
        }
    }
    return a0;
}

const char* atom_to_str(dbr_aconf const& aconf, int n)
{
    static char str[50];
    dbr_atom const& a = aconf.atoms[n];
    const char *fmt = aconf.reduced_coordinates ? "#%d %s (%.3f %.3f %.3f)"
                                                : "#%d %s (%.1f %.1f %.1f)";
    snprintf(str, 50, fmt, n, a.name, a.xyz[0], a.xyz[1], a.xyz[2]);
    return str;
}

// find all atoms that are on a line parallel to axis `dim' (x axis if dim=0, 
// y if dim=1, z if dim=2) as atom a0 (+/- epsilon). Returns sorted list of
// distances from a0 to these atoms, up to half of the PBC in `dim' direction.
vector<double> nabe_in_direction_distances(dbr_aconf const& aconf,
                                           int a0, int dim,
                                           double lo_lim, double hi_lim,
                                           double epsilon)
{
    vector<double> img_dist;
    double pbcd = get_pbc(aconf, dim);
    dbr_real* xyz0 = aconf.atoms[a0].xyz;
    int dim1 = (dim + 1) % 3;
    int dim2 = (dim + 2) % 3;
    if (verbosity > 2) {
        printf("searching for atoms with %c=%g and %c=%g +/- %g,",
                'x'+dim1, xyz0[dim1], 'x'+dim2, xyz0[dim2], epsilon);
        printf(" %g < |%c-%g| < %g\n", lo_lim, 'x'+dim, xyz0[dim], hi_lim);
    }
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
        if (fabs(xyz[dim1] - xyz0[dim1]) > epsilon ||
                fabs(xyz[dim2] - xyz0[dim2]) > epsilon ||
                i == a0)
            continue;
        double diff = dist_forward(xyz[dim], xyz0[dim], pbcd);
        bool ok = (diff > lo_lim && diff < hi_lim);
        if (verbosity > 2)
            printf("distance %g to %s  %s\n", diff, atom_to_str(aconf, i),
                                              ok ? "ok" : "rejected");
        if (ok)
            img_dist.push_back(diff);
    }
    sort(img_dist.begin(), img_dist.end());
    return img_dist;
}

// If single_translation is true:
//   check if there is a translational symmetry with translation
//   vector that has direction `dim' and length delta that maps 1:1 atoms with
//   `dim' coordinates in range (x0-delta, x0) to atoms in range (x0, x0+delta).
// Otherwise:
//   check if all atoms have images in translation, i.e. check
//   periodicity (in given direction) in all the system.
// slab.delta is the length of the translation vector
double check_symmetry_in_distance(CellMethod const& cm, Slab const& slab,
                                  double epsilon, bool single_translation)
{
    if (verbosity > 0) {
        printf("checking translation |t|=%g ... ", slab.delta);
        fflush(stdout);
    }
    dbr_aconf const& aconf = cm.aconf();
    double pbcd = get_pbc(aconf, slab.dim);
    StdDev avg_w;
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
        dbr_real x[3] = { xyz[0], xyz[1], xyz[2] };
        if (single_translation) { // checking one translation
            double d = dist_forward(xyz[slab.dim], slab.x0, pbcd);
            assert (d >= 0);
            assert (d <= pbcd);
            // checking only atoms in distance delta from x0 (on both sides)
            if (d > slab.delta && d <= pbcd - slab.delta)
                continue;

            if (d < slab.delta) {
                x[slab.dim] -= slab.delta;
                if (x[slab.dim] < 0)
                    x[slab.dim] += pbcd;
            }
            else {
                x[slab.dim] += slab.delta;
                if (x[slab.dim] >= pbcd)
                    x[slab.dim] -= pbcd;
            }
        }
        else { // checking if all the system is periodic
            x[slab.dim] += slab.delta;
            if (x[slab.dim] >= pbcd)
                x[slab.dim] -= pbcd;
        }

        dbr_atom const* atom = cm.get_atom_at(x, epsilon);
        if (atom == NULL || strcmp(atom->name, aconf.atoms[i].name) != 0) {
            if (verbosity > 0)
                printf("failed (atom %s)\n", atom_to_str(aconf, i));
            if (verbosity > 1) {
                int a = cm.find_nearest(x);
                if (a >= 0)
                    printf(" the nearest atom (|r|=%g) is %s\n",
                           cm.get_dist(x, aconf.atoms[a].xyz),
                           atom_to_str(aconf, a));
            }
            return 0.;
        }
        // calculate average length of translation vector in direction dim
        double exact_w = dist_forward(atom->xyz[slab.dim], xyz[slab.dim], pbcd);
        avg_w.add_x(exact_w);
    }
    if (verbosity > 0)
        printf("ok (t%c = %g +/- %g)\n", 'x'+slab.dim, avg_w.mean(),
                                                       avg_w.stddev());
    return avg_w.mean();
}

double search_for_translation(dbr_aconf const& aconf,
                              CellMethod const& cm,
                              gengetopt_args_info const& args,
                              int dim, double val)
{
    int a0 = select_atom(aconf, dim, val);
    if (verbosity > 0)
        printf("atom0: %s\n", atom_to_str(aconf, a0));
    double lo_lim = args.min_delta_given ? args.min_delta_arg : 0.;
    double hi_lim = args.max_delta_given ? args.max_delta_arg
                                         : get_pbc(aconf, dim) / 2.;
    double eps = args.epsilon_arg;
    bool single_trans = !args.periodic_given;
    // possible lengths of translation vectors, based on one atom's neighbors
    vector<double> tr = nabe_in_direction_distances(aconf, a0, dim,
                                                    lo_lim, hi_lim, eps);
    // return the first (smallest) translation vector found
    Slab slab;
    slab.dim = dim;
    slab.x0 = val;
    for (vector<double>::const_iterator d = tr.begin(); d != tr.end(); ++d) {
        slab.delta = *d;
        double t = check_symmetry_in_distance(cm, slab, eps, single_trans);
        if (t != 0.) {
            if (!single_trans && verbosity > -1)
                printf("PBC[%d] / |t| = %g\n", dim, get_pbc(aconf, dim) / t);
            return t;
        }
    }
    return 0.;
}

template <typename T>
struct index_sorter
{
    index_sorter(vector<T> const& values) : values_(values) {}
    bool operator() (int a, int b) const { return values_[a] < values_[b]; }
    vector<T> const& values_;
};

struct XYZ
{
    dbr_real x, y, z;
};

// search for translational symmetry in all directions
void find_trans_sym(dbr_aconf const& aconf, Slab const& slab,
                    gengetopt_args_info const& args)
{
    CellMethod cm(aconf, args.min_cell_arg);

    // print min. interatomic distance
    double min_dist = cm.find_min_distance();
    double eps = args.epsilon_arg;
    if (eps >= min_dist / 2.) {
        printf("WARNING: Min. interatomic distance is: %g.\n", min_dist);
        printf("WARNING: Epsilon should be smaller than %g.\n", min_dist/2);
    }
    double lo_lim = args.min_delta_given ? args.min_delta_arg : 0.;
    double hi_lim = args.max_delta_given ? args.max_delta_arg
                                         : get_pbc(aconf, slab.dim) / 2.;

    vector<double> dd = cm.get_image_distances(0, lo_lim, hi_lim, eps);

    vector<int> pp;
    for (int i = 0; i != aconf.n; ++i)
        if (dd[i] != 0.)
            pp.push_back(i);
    sort(pp.begin(), pp.end(), index_sorter<double>(dd));

    vector<XYZ> tt;

    dbr_real* xyz0 = aconf.atoms[0].xyz;
    for (vector<int>::const_iterator i = pp.begin(); i != pp.end(); ++i) {
        dbr_real r[3];
        cm.get_shortest_vec(xyz0, aconf.atoms[*i].xyz, r);

        // average the vector
        StdDev avg_r[3];
        for (int j = 0; j != aconf.n; ++j) {
            dbr_real delta[3];
            dbr_atom const* img = cm.get_image(r, j, eps, delta);
            assert(img != NULL);
            for (int k = 0; k < 3; ++k)
                avg_r[k].add_x(r[k] + delta[k]);
        }
        XYZ t = { avg_r[0].mean(), avg_r[1].mean(), avg_r[2].mean() };

        // make the first index non-negative
        if (t.x < 0) {
            t.x = -t.x;
            t.y = -t.y;
            t.z = -t.z;
        }

        // check if it's not n * another translation
        bool is_dup = false;
        for (vector<XYZ>::const_iterator j = tt.begin(); j != tt.end(); ++j) {
            // calculate average multiplier
            double m = (t.x + t.y + t.z) / (j->x + j->y + j->z);
            // check if this is a multiplier for all coordinates
            if (fabs(m * j->x - t.x) < eps
                    && fabs(m * j->y - t.y) < eps
                    && fabs(m * j->z - t.z) < eps) {
                is_dup = true;
                break;
            }
        }
        if (is_dup)
            continue;

        printf("T = (% f % f % f) +/- (%g %g %g)\n",
               t.x, t.y, t.z,
               avg_r[0].stddev(), avg_r[1].stddev(), avg_r[2].stddev());
        tt.push_back(t);
#if 0
        // remove translations that are n * base_translation (in PBC)
        while (n != 0) { // n == 0 means we got back to the original atom
            dbr_real* xyz = aconf.atoms[n].xyz;
            dbr_real x[3] = { xyz[0]+r[0], xyz[1]+r[1], xyz[2]+r[2] };
            n = cm.find_nearest(x);
            dd[n] = 0.;
        }
#endif
    }
}

// If output is false:
//   delete atoms with `dim' coord in <x0, x0+delta). 
//   If delta < 0, deleted range is <x0-|delta|, x0).
// Otherwise delete the rest of atoms.
void delete_atoms(dbr_aconf& aconf, Slab const& slab, bool del_outside)
{
    //printf("slab: %c=%g, delta: %g\n", 'x'+slab.dim, slab.x0, slab.delta);
    int n = 0;
    double pbcd = get_pbc(aconf, slab.dim);
    double delta = fabs(slab.delta);
    double x0 = slab.delta >= 0 ? slab.x0 : slab.x0 - delta;
    for (int i = 0; i != aconf.n; ++i) {
        double d = dist_forward(aconf.atoms[i].xyz[slab.dim], x0, pbcd);
        bool atom_in_slab = (d < delta);
        // if del_outside, we keep (i.e. don't delete) atoms in the slab
        // otherwise, we keep atom out of the slab
        if ((del_outside && atom_in_slab) || (!del_outside && !atom_in_slab)){
            dbr_copy_atom(aconf, i, n);
            ++n;
        }
    }
    int deleted = aconf.n - n;
    aconf.n = n;
    if (!aconf.auxiliary.empty())
        aconf.auxiliary.resize(aconf.n);
    if (verbosity > -1)
        printf("%d atoms were deleted.\n", deleted);
}

void add_vacuum(dbr_aconf& aconf, Slab const& slab)
{
    for (int i = 0; i != aconf.n; ++i)
        if (aconf.atoms[i].xyz[slab.dim] > slab.x0)
            aconf.atoms[i].xyz[slab.dim] += slab.delta;
    expand_pbc(aconf, slab.dim, slab.delta);
}

void add_copy(dbr_aconf& aconf, Slab const& slab,
              CellMethod const& cm, double epsilon, int N)
{
    // bookmark atoms that are to be copied
    vector<int> orig;
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
        double pbcd = get_pbc(aconf, slab.dim);
        double dist = dist_forward(slab.x0, xyz[slab.dim], pbcd);
        if (dist < slab.delta - epsilon)
            orig.push_back(i);
        // to avoid overlapping atoms or gaps, we must ensure that all images
        // of bookmarked atoms are moved, and atoms which images are not to be
        // moved are not copied.
        else if (dist < slab.delta + epsilon) {
            dbr_real x[3] = { xyz[0], xyz[1], xyz[2] };
            x[slab.dim] += slab.delta;
            dbr_atom const* img = cm.get_atom_at(x, epsilon);
            if (img && img->xyz[slab.dim] >= slab.x0)
                orig.push_back(i);
        }
    }

    if (!orig.empty()) {
        // resize atoms table
        dbr_atom* new_atoms = new dbr_atom[aconf.n + N * orig.size()];
        memcpy(new_atoms, aconf.atoms, sizeof(dbr_atoms) * aconf.n);
        delete [] aconf.atoms;
        aconf.atoms = new_atoms;
    }

    double delta = slab.delta;
    for (int i = 0; i < N; ++i) {
        add_vacuum(aconf, slab);

        // if delta is relative, we rescale it to keep the same absolute value
        if (aconf.reduced_coordinates)
            delta /= 1 + delta;

        // copy bookmarked atoms into the just created vacuum
        for (size_t i = 0; i != orig.size(); ++i) {
            dbr_copy_atom(aconf, orig[i], aconf.n + i);
            aconf.atoms[aconf.n + i].xyz[slab.dim] += delta;
        }
        aconf.n += orig.size();
    }
    if (verbosity > -1)
        printf("%d atoms were added.\n", int(N * orig.size()));
}

// Don't search for translation vectors, just multiply the configuration
// by integer number in all directions.
void multiply_pbc(dbr_aconf& aconf, int x, int y, int z)
{
    int new_n = aconf.n * x * y * z;
    dbr_atom* new_atoms = new dbr_atom[new_n];
    int mult[3] = {x, y, z};
    for (int ix = 0; ix < x; ++ix)
        for (int iy = 0; iy < y; ++iy)
            for (int iz = 0; iz < z; ++iz) {
                int img_start = (ix * y * z + iy * z + iz) * aconf.n;
                dbr_real r[3];
                dbr_real img[3] = {ix, iy, iz};
                dbr_vec3_mult_pbc(img, aconf.pbc, r);
                for (int i = 0; i < aconf.n; ++i) {
                    dbr_atom const& source = aconf.atoms[i];
                    dbr_atom& dest = new_atoms[i + img_start];
                    strcpy(dest.name, source.name);
                    if (aconf.reduced_coordinates) {
                        for (int i = 0; i < 3; ++i)
                            dest.xyz[i] = (source.xyz[i] + img[i]) / mult[i];
                    }
                    else {
                        for (int i = 0; i < 3; ++i)
                            dest.xyz[i] = source.xyz[i] + r[i];
                    }
                }
            }

    aconf.atoms = new_atoms;
    aconf.pbc.v00 *= x; aconf.pbc.v01 *= x; aconf.pbc.v02 *= x;
    aconf.pbc.v10 *= y; aconf.pbc.v11 *= y; aconf.pbc.v12 *= y;
    aconf.pbc.v20 *= z; aconf.pbc.v21 *= z; aconf.pbc.v22 *= z;
    string comment = "configuration multiplied by " + S(x) + "x" + S(y)
                     + "x" + S(z);
    aconf.comments.insert(aconf.comments.begin(), comment);
    aconf.auxiliary.resize(new_n);
    for (int img = 1; img < x*y*z; ++img)
        for (int i = 0; i != aconf.n; ++i)
            aconf.auxiliary[img * aconf.n + i] = aconf.auxiliary[i];
    aconf.n = new_n;
}

bool multiply_pbc(dbr_aconf& aconf, const char* mult)
{
    // parse string
    const char *start = mult;
    char *endptr = NULL;
    int x = strtol(start, &endptr, 10);
    if (endptr == start || *endptr != 'x')
        return false;
    start = endptr + 1;
    int y = strtol(start, &endptr, 10);
    if (endptr == start || *endptr != 'x')
        return false;
    start = endptr + 1;
    int z = strtol(start, &endptr, 10);
    if (endptr == start || *endptr != '\0')
        return false;
    if (x == 0 || y == 0 || z == 0)
        return false;
    multiply_pbc(aconf, x, y, z);
    return true;
}

bool shift_under_pbc(dbr_aconf& aconf, const char* arg)
{
    // parse string
    const char *start = arg;
    char *endptr = NULL;
    double d[3];
    d[0] = strtod(start, &endptr);
    if (endptr == start || *endptr != ',')
        return false;
    start = endptr + 1;
    d[1] = strtod(start, &endptr);
    if (endptr == start || *endptr != ',')
        return false;
    start = endptr + 1;
    d[2] = strtod(start, &endptr);
    if (endptr == start || *endptr != '\0')
        return false;

    // shift atoms
    for (int i = 0; i != aconf.n; ++i)
        for (int j = 0; j < 3; ++j)
            aconf.atoms[i].xyz[j] += d[j];
    wrap_to_pbc(aconf);
    return true;
}

double find_extreme_coord(dbr_aconf const& aconf, int dim, bool upper)
{
    vector<dbr_real> v(aconf.n);
    for (int i = 0; i != aconf.n; ++i)
        v[i] = aconf.atoms[i].xyz[dim];
    sort(v.begin(), v.end());
    double pbcd = get_pbc(aconf, dim);
    int pos = 0;
    double max_delta = v.front() + pbcd - v.back();
    for (int i = 1; i != aconf.n; ++i)
        if (v[i] - v[i-1] > max_delta) {
            pos = i;
            max_delta = v[i] - v[i-1];
        }
    // if the vacuum is found between x1 and x2, x1<x2, x2 is the lower bound
    // of the material slab, and x1 -- the upper bound
    if (upper)
        pos = (pos > 0 ? pos-1 : aconf.n-1);
    double margin = (upper ? 1e-9 : -1e-9);
    return v[pos] + margin;
}

string argv_as_str(int argc, char **argv)
{
    string s;
    for (int i = 0; i < argc; ++i) {
        if (i != 0)
            s += " ";
        // TODO
        //if (should be quoted)
        //    s += str('"') + argv[i] + '"';
        //else
            s += argv[i];
    }
    return s;
}

int main(int argc, char **argv)
{
    // processing cmd line args
    gengetopt_args_info args;
    if (cmdline_parser(argc, argv, &args) != 0)
        return EXIT_FAILURE;

    if (args.show_examples_given) {
        int len = sizeof(usage_examples)/sizeof(usage_examples[0]);
        for (int i = 0; i < len; ++i)
            printf("%s\n", usage_examples[i]);
        return EXIT_SUCCESS;
    }

    verbosity = args.verbose_given; // global variable

    if (args.x_given + args.y_given + args.z_given + args.bound_given> 1) {
        fprintf(stderr, "Only one of -x, -y, -z, -b options can be given "
                "(-h will show all valid options).\n");
        return EXIT_FAILURE;
    }
    if (!args.x_given && !args.y_given && !args.z_given && !args.bound_given
            && !args.multiply_given && !args.shift_given
            && !args.find_trans_given && !args.merge_given && !args.t1_given) {
        fprintf(stderr, "One of -x, -y, -z, -b, -N, -S, -F options must be "
                        "given (-h will show all valid options).\n");
        return EXIT_FAILURE;
    }
    if (args.in_place_given && args.output_given) {
        fprintf(stderr, "Only one of -i and -o options can be given "
                "(-h will show explaination).\n");
        return EXIT_FAILURE;
    }
    if (args.add_vacuum_given && args.add_copy_given) {
        fprintf(stderr, "Only one of add-copy and add_vacuum parameters "
                "can be given.\n");
        return EXIT_FAILURE;
    }

    if (args.inputs_num != 1) {
        fprintf(stderr, "Exactly one input file should be given.\n");
        return EXIT_FAILURE;
    }

    // reading configuration file
    bool reduced_coords = args.reduced_given;
    dbr_aconf aconf = read_atoms_from_file(args.inputs[0], reduced_coords);

    Slab slab;
    slab.dim = -1;
    slab.x0 = 0.;
    slab.delta = 0.;
    if (args.x_given) {
        slab.dim = 0;
        slab.x0 = args.x_arg;
    }
    else if (args.y_given) {
        slab.dim = 1;
        slab.x0 = args.y_arg;
    }
    else if (args.z_given) {
        slab.dim = 2;
        slab.x0 = args.z_arg;
    }
    else if (args.bound_given) {
        bool upper;
        switch (args.bound_arg) {
            case bound_arg_x: slab.dim = 0; upper = false; break;
            case bound_arg_y: slab.dim = 1; upper = false; break;
            case bound_arg_z: slab.dim = 2; upper = false; break;
            case bound_arg_X: slab.dim = 0; upper = true; break;
            case bound_arg_Y: slab.dim = 1; upper = true; break;
            case bound_arg_Z: slab.dim = 2; upper = true; break;
        }
        slab.x0 = find_extreme_coord(aconf, slab.dim, upper);
        if (verbosity > 0)
            printf("Plane set to %c=%g.\n", 'x'+slab.dim, slab.x0);
    }

    if (slab.dim >= 0) {
        double pbcd = get_pbc(aconf, slab.dim);
        if (verbosity > 1)
            printf("PBC[%d]=%g\n", slab.dim, pbcd);
        if (slab.x0 < 0) {
            fprintf(stderr,
             "The given point (%g) should be in PBC (non-negative)\n", slab.x0);
            return EXIT_FAILURE;
        }
        if (slab.x0 > pbcd) {
            fprintf(stderr,
             "The given point (%g) should be in PBC ( < %g)\n", slab.x0, pbcd);
            return EXIT_FAILURE;
        }

        // put atoms into cells
        CellMethod cm(aconf, args.min_cell_arg);

        // print min. interatomic distance
        double min_dist = cm.find_min_distance();
        printf("Min. interatomic distance: %g\n", min_dist);

        // if epsilon is larger than half of cell, increase cell sizes
        if (args.epsilon_arg >= cm.half_cell()) {
            fprintf(stderr, "Epsilon too large.\n");
            return EXIT_FAILURE;
        }

        if (args.width_given)
            slab.delta = args.width_arg;
        else {
            slab.delta = search_for_translation(aconf, cm, args,
                                                slab.dim, slab.x0);
            if (verbosity > -1) {
                if (slab.delta != 0.)
                    printf("translation length: %g\n", slab.delta);
                else {
                    printf("translation not found\n");
                    return EXIT_FAILURE;
                }
            }
        }
        if (args.delete_given)
            delete_atoms(aconf, slab, /*del_outside=*/false);
        if (args.add_vacuum_given)
            add_vacuum(aconf, slab);
        if (args.add_copy_given)
            add_copy(aconf, slab, cm, args.epsilon_arg,
                     args.add_copy_arg);
        if (args.cut_given) {
            delete_atoms(aconf, slab, /*del_outside=*/true);
            expand_pbc(aconf, slab.dim, slab.delta-pbcd);
        }
    }

    if (args.multiply_given)
        multiply_pbc(aconf, args.multiply_arg);

    if (args.shift_given)
        shift_under_pbc(aconf, args.shift_arg);

    if (args.find_trans_given)
        find_trans_sym(aconf, slab, args);

    if (args.t1_given)
        transform1(aconf, args);

    if (args.merge_given)
        merge_atoms(aconf, args);

    aconf.comments.insert(aconf.comments.begin(), argv_as_str(argc, argv));
    if (args.in_place_given)
        write_file_with_atoms(aconf, aconf.orig_filename);
    else if (args.output_given)
        write_file_with_atoms(aconf, args.output_arg);

    return EXIT_SUCCESS;
}

