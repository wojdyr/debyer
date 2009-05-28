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
// $Id: debyer.ggo 112 2009-04-14 00:44:25Z wojdyr $

// TODO: 
// * it doesn't work when there are two cells in given direction,
//   ./dbr_extend -r  -e0.02 -c2.99 -M0.03  -v -v  -z0.5 mydisl.cfg.gz

#include <algorithm>
#include <cstring>

#include "fileio.h"
#include "iloops.h"
#include "utils.h"
#include "extend_cmd.h"

using namespace std;

static int verbosity = 0;

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

// expand PBC, but don't move atoms, i.e. change reduced coordinates
// in case of reduced coords, delta is also reduced
void expand_pbc(dbr_aconf& aconf, int dim, double delta)
{
    double* pp = get_pbc_ptr(aconf, dim);
    double old_p = *pp;
    double new_p = aconf.reduced_coordinates ? old_p * (1. + delta)
                                             : old_p + delta;

    if (aconf.reduced_coordinates)
        for (int i = 0; i != aconf.n; ++i)
            aconf.atoms[i].xyz[dim] /= (1. + delta);

    *pp = new_p;

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

    double get_dist(dbr_real const* a, dbr_real const* b) const;

    double half_cell() const
        { return min(min(cell_size(0), cell_size(1)), cell_size(2)); }

    double const* pbc() const { return pbc_; }
    dbr_aconf const& aconf() const { return aconf_; }

//private:
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

// gets cell and atom coordinates; 
// if the cell indices are correct returns pointer to data of this cell,
// otherwise wraps the cell indices and shifts the atom.
int const* CellMethod::adjust_atom(int idx0, int idx1, int idx2,
                                   dbr_real* xyz) const
{
    int a[3] = { idx0, idx1, idx2 };
    for (int i = 0; i < 3; ++i)
        if (a[i] == -1) {
            a[i] += n_[i];
            xyz[i] += pbc_[i];
        }
        else if (a[i] == n_[i]) {
            a[i] = 0;
            xyz[i] -= pbc_[i];
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


dbr_atom const*
CellMethod::get_atom_at(dbr_real const* xyz, double epsilon) const
{
    // pre: epsilon < r/2, where r=pbc_[i]/n_[i]
    int idx[3];
    find_cell_indices(xyz, idx);

    int *p = cell_start(idx);
    while (*p != -1) {
        if (get_sq_dist(xyz, aconf_.atoms[*p].xyz) < epsilon * epsilon)
            return &aconf_.atoms[*p];
        ++p;
    }

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
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
        if (fabs(xyz[dim1] - xyz0[dim1]) > epsilon ||
                fabs(xyz[dim2] - xyz0[dim2]) > epsilon ||
                i == a0)
            continue;
        double diff = dist_forward(xyz[dim], xyz0[dim], pbcd);
        if (verbosity > 2)
            printf("distance %g to %s\n", diff, atom_to_str(aconf, i));
        if (diff > lo_lim && diff < hi_lim)
            img_dist.push_back(diff);
    }
    sort(img_dist.begin(), img_dist.end());
    return img_dist;
}

// check if there is a translational symmetry with translation
// vector that has direction `dim' and length w that maps 1:1 atoms with
// `dim' coordinates in range (x0-w, x0) to atoms in range (x0, x0+w).
bool check_symmetry_in_distance(CellMethod const& cm,
                                int dim, double x0, double w, double epsilon)
{
    if (verbosity > 0) {
        printf("checking translations |v|=%g ... ", w);
        fflush(stdout);
    }
    dbr_aconf const& aconf = cm.aconf();
    double pbcd = get_pbc(aconf, dim);
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
        double d = dist_forward(xyz[dim], x0, pbcd);
        assert (d >= 0);
        assert (d <= pbcd);
        if (d < w || d > pbcd - w) {
            dbr_real x[3] = { xyz[0], xyz[1], xyz[2] };
            if (d < w) {
                x[dim] -= w;
                if (x[dim] < 0)
                    x[dim] += pbcd;
            }
            else {
                x[dim] += w;
                if (x[dim] >= pbcd)
                    x[dim] -= pbcd;
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
                return false;
            }
        }
    }
    if (verbosity > 0)
        printf("ok\n");
    return true;
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
    double hi_lim = args.max_delta_given ? args.min_delta_arg
                                         : get_pbc(aconf, dim) / 2.;
    double eps = args.epsilon_arg;
    // possible lengths of translation vectors, based on one atom's neighbors
    vector<double> tr = nabe_in_direction_distances(aconf, a0, dim,
                                                    lo_lim, hi_lim, eps);
    // return the first (smallest) translation vector found
    for (vector<double>::const_iterator d = tr.begin(); d != tr.end(); ++d)
        if (check_symmetry_in_distance(cm, dim, val, *d, eps))
            return *d;
    return 0.;
}

// delete atoms with `dim' coord in <x0, x0+delta). 
// If delta < 0, deleted range is <x0-|delta|, x0).
void delete_atoms(dbr_aconf& aconf, int dim, double x0, double delta)
{
    int n = 0;
    double pbcd = get_pbc(aconf, dim);
    if (delta < 0) {
        x0 -= delta;
        delta = -delta;
    }
    for (int i = 0; i != aconf.n; ++i)
        if (dist_forward(aconf.atoms[i].xyz[dim], x0, pbcd) >= delta) {
            dbr_copy_atom(aconf, i, n);
            ++n;
        }
    int deleted = aconf.n - n;
    aconf.n = n;
    if (verbosity > 0)
        printf("%d atoms were deleted.\n", deleted);
}

void add_vacuum(dbr_aconf& aconf, int dim, double x0, double delta)
{
    for (int i = 0; i != aconf.n; ++i)
        if (aconf.atoms[i].xyz[dim] > x0)
            aconf.atoms[i].xyz[dim] += delta;
    expand_pbc(aconf, dim, delta);
}

void add_copy(dbr_aconf& aconf, int dim, double x0, double delta,
              CellMethod const& cm, double epsilon)
{
    // bookmark atoms that are to be copied
    double pbcd = get_pbc(aconf, dim);
    vector<int> orig;
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
        double dist = dist_forward(x0, xyz[dim], pbcd);
        if (dist < delta - epsilon)
            orig.push_back(i);
        // to avoid overlapping atoms or gaps, we must ensure that all images
        // of bookmarked atoms are moved, and atoms which images are not to be
        // moved are not copied.
        else if (dist < delta + epsilon) {
            dbr_real x[3] = { xyz[0], xyz[1], xyz[2] };
            x[dim] += delta;
            dbr_atom const* img = cm.get_atom_at(x, epsilon);
            if (img && img->xyz[dim] >= x0)
                orig.push_back(i);
        }
    }

    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
        if (xyz[dim] >= x0)
            xyz[dim] += delta;
    }
    expand_pbc(aconf, dim, delta);

    // if delta is relative, we rescale it to keep the same absolute value
    if (aconf.reduced_coordinates)
        delta /= 1 + delta;

    if (!orig.empty()) {
        // resize atoms table
        dbr_atom* new_atoms = new dbr_atom[aconf.n + orig.size()];
        memcpy(new_atoms, aconf.atoms, sizeof(dbr_atoms) * aconf.n);
        delete [] aconf.atoms;
        aconf.atoms = new_atoms;

        // copy bookmarked atoms into the just created vacuum
        for (size_t i = 0; i != orig.size(); ++i) {
            dbr_copy_atom(aconf, orig[i], aconf.n + i);
            aconf.atoms[aconf.n + i].xyz[dim] += delta;
        }
        aconf.n += orig.size();
        if (verbosity > 0)
            printf("%d atoms were added.\n", (int) orig.size());
    }
}

int main(int argc, char **argv)
{
    // processing cmd line args
    gengetopt_args_info args;
    if (cmdline_parser(argc, argv, &args) != 0)
        return -1;

    verbosity = args.verbose_given; // global variable

    if (args.x_given + args.y_given + args.z_given != 1) {
        fprintf(stderr, "One of -x, -y, -z parameters should be given "
                "(-h will show all valid options).\n");
        return -1;
    }
    if (args.add_vacuum_given && args.add_copy_given) {
        fprintf(stderr, "Only one of add-copy and add_vacuum parameters "
                "can be given.\n");
        return -1;
    }

    if (args.inputs_num != 1) {
        fprintf(stderr, "Exactly one input file should be given.\n");
        return -1;
    }

    int dim = -1;
    double val = 0;
    if (args.x_given) {
        dim = 0;
        val = args.x_arg;
    }
    else if (args.y_given) {
        dim = 1;
        val = args.y_arg;
    }
    else if (args.z_given) {
        dim = 2;
        val = args.z_arg;
    }

    // reading configuration file
    bool reduced_coords = args.reduced_given;
    dbr_aconf aconf = read_atoms_from_file(args.inputs[0], reduced_coords);

    double pbcd = get_pbc(aconf, dim);
    if (verbosity > 1)
        printf("PBC[%d]=%g\n", dim, pbcd);
    if (val < 0) {
        fprintf(stderr,
                "The given point (%g) should be in PBC (non-negative)\n", val);
        return -1;
    }
    if (val > pbcd) {
        fprintf(stderr,
                "The given point (%g) should be in PBC ( < %g)\n", val, pbcd);
        return -1;
    }

    // put atoms into cells
    CellMethod cm(aconf, args.min_cell_arg);

    // if epsilon is larger than half of cell, increase cell sizes
    if (args.epsilon_arg >= cm.half_cell()) {
        fprintf(stderr, "Epsilon too large.\n");
        return -1;
    }

    double delta = 0.;
    if (args.width_given)
        delta = args.width_arg;
    else {
        delta = search_for_translation(aconf, cm, args, dim, val);
        if (delta != 0.)
            printf("translation length: %g\n", delta);
        else
            printf("translation not found\n");
    }
    if (args.delete_given)
        delete_atoms(aconf, dim, val, delta);
    if (args.add_vacuum_given)
        add_vacuum(aconf, dim, val, delta);
    if (args.add_copy_given)
        add_copy(aconf, dim, val, delta, cm, args.epsilon_arg);

    if (args.in_place_given)
        write_file_with_atoms(aconf, aconf.orig_filename);
    else if (args.output_given)
        write_file_with_atoms(aconf, args.output_arg);

    return 0;
}

