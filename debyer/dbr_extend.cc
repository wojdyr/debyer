// dbr_extend -- a tool originally written to extend a system (configuration
// for atomistic simulation) in x, y or z direction by increasing PBC and
// filling the new space with periodic structure, a copy of the structure
// around it. Now it can also do several other things.
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

// TODO: 
// * option -cN with N>1 should change PBC only once

#include <algorithm>
#include <cstring>

#include "fileio.h"
#include "utils.h"
#include "extend_cmd.h"
#include "cells.h"

using namespace std;

static int verbosity = 0;

const char* usage_examples[] = {
"                              USAGE EXAMPLES",
"",
"dbr_extend -z41.5 -e0.2 -vvv file.cfg",
"  Try to find periodicity of the structure in the z direction starting from",
"  z0=41.5. It tries to find z1 such that each atom with z0 < z < z1 has",
"  a periodic image with x'=x, y'=y, z'=z+delta, delta=z1-z0.",
"  Comparisons of coordinates are done with epsilon=0.2.",
"",
"dbr_extend -c10 -z41.5 -e0.2 -i file.cfg",
"  Extend PBC box by in the z direction by 10 times the value reported",
"  from the command above, copy atoms to the newly created space,",
"  write the configuration back to the same file.",
"",
"dbr_extend -bz -w3 -d -i file.cfg",
"dbr_extend -bZ -w-3 -d -i file.cfg",
"  (Initially, file.cfg contained a slab with surfaces normal to z axis.)",
"  Delete surfaces (3A deep) of the slab.",
"",
"dbr_extend -v -bz -w2 -a3 -o tmp2.cfg tmp1.cfg",
"dbr_extend -v -bZ -w-2 -a3 -o tmp3.cfg tmp2.cfg",
"  (Initially, tmp1.cfg contained a slab with surfaces normal to z axis.)",
"  Extend the slab in the z direction, 3A from each surface.",
"",
"dbr_extend -S0,0,0.5 -r -i file.cfg",
"  Shift object under PBC, by half of the PBC box, in the z direction,",
"  write the configuration back to file.cfg.",
"",
"dbr_extend -N1x2x1 -o out.cfg file.cfg",
"  Duplicate the system in the y direction (create a supercell).",
"",
"dbr_extend -v -z-7 -w12. --density file.cfg",
"  Calculate (in a smart way) numeric density of the slab defined by planes",
"  z=-7 and z=5 (its in PBC, so it's continuus region that includes z=0).",
"  Designed to calculate density of a GB in bicrystal geometry.",
"",
"dbr_extend --resize=ref.cfg -o output.cfg input.cfg",
"  Resize the PBC box, make it the same as the size of the file ref.cfg.",
"  Atomic positions are scaled with the box."
};


struct Slab
{
    int dim;
    double x0;
    double delta;
};

struct DeltaLim
{
    double lo, hi;
};

struct Trans // translation vector
{
    double dist;
    double r[3];
    double stddev[3];
    bool operator<(Trans const& t) const { return dist < t.dist; }
};


bool parse_n_numbers(const char* arg, int n, double* d)
{
    const char *start = arg;
    char *endptr = NULL;

    for (int i = 0; i != n; ++i) {
        d[i] = strtod(start, &endptr);
        if (endptr == start || *endptr != (i < n-1 ? ',' : '\0'))
            return false;
        start = endptr + 1;
    }
    return true;
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


class SlabInPBC
{
public:
    SlabInPBC(dbr_aconf const& aconf, Slab const& slab)
        : dim_(slab.dim),
          pbcd_(get_pbc(aconf, slab.dim)),
          delta_(fabs(slab.delta)),
          x0_(wrapped(slab.delta >= 0 ? slab.x0 : slab.x0 - delta_)) {}

    // if reduced coordinates are used, pbcd_==1 and x-floor(x) is returned
    double wrapped(double x) const { return x - floor(x / pbcd_) * pbcd_; }

    bool has(dbr_real const* xyz) const
        { return dist_forward(xyz[dim_], x0_, pbcd_) < delta_; }
    bool ok() const { return dim_ >= 0; }
    int dim() const { return dim_; }
    double x0() const { return x0_; }
    double delta() const { return delta_; }

private:
    int dim_;
    double pbcd_, delta_, x0_;
};

int get_any_atom_in_slab(dbr_aconf const& aconf, SlabInPBC const& pslab)
{
    if (pslab.dim() == -1 && aconf.n > 0)
        return 0;
    for (int i = 0; i != aconf.n; ++i)
        if (pslab.has(aconf.atoms[i].xyz))
            return i;
    return -1;
}

// returns the number of removed atoms
int merge_atoms(dbr_aconf& aconf, gengetopt_args_info const& args, bool quiet)
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
    if (!quiet) {
        double ratio = double(aconf.n) / new_n;
        printf("merging, i.e. removing duplicates: %g%% (1/%g) of atoms left\n",
                100./ratio, ratio);
    }

    int n = 0;
    for (int i = 0; i != aconf.n; ++i) {
        if (!is_null(aconf.atoms[i])) {
            dbr_copy_atom(aconf, i, n);
            ++n;
        }
    }
    assert (n == new_n);
    if (!aconf.auxiliary.empty())
        aconf.auxiliary.resize(new_n);

    int n_removed = aconf.n - new_n;
    if (!quiet) {
        string comment = "merging atoms, " + S(aconf.n) + " -> " + S(new_n);
        aconf.comments.insert(aconf.comments.begin(), comment);
    }
    aconf.n = new_n;
    return n_removed;
}

void rotate_atoms(dbr_aconf& aconf, double angle, double axis[3])
{
    double rot[3][3];
    rodrigues(angle, axis, rot);
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* x = aconf.atoms[i].xyz;
        // x = rot . x
        dbr_real r[3];
        for (int j = 0; j < 3; ++j)
            r[j] = rot[j][0] * x[0] + rot[j][1] * x[1] + rot[j][2] * x[2];
        for (int j = 0; j < 3; ++j)
            x[j] = r[j];
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
    merge_atoms(aconf, args, false);
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
                                           DeltaLim const& lim,
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
        printf(" %g < |%c-%g| < %g\n", lim.lo, 'x'+dim, xyz0[dim], lim.hi);
    }
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
        if (fabs(xyz[dim1] - xyz0[dim1]) > epsilon ||
                fabs(xyz[dim2] - xyz0[dim2]) > epsilon ||
                i == a0)
            continue;
        double diff = dist_forward(xyz[dim], xyz0[dim], pbcd);
        bool ok = (diff > lim.lo && diff < lim.hi);
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
        bool forward = true;
        if (single_translation) { // checking one translation
            double d = dist_forward(xyz[slab.dim], slab.x0, pbcd);
            assert (d >= 0);
            assert (d <= pbcd);
            // checking only atoms in distance delta from x0 (on both sides)
            if (d > slab.delta && d <= pbcd - slab.delta)
                continue;
            if (d < slab.delta)
                forward = false;
        }

        dbr_real x[3] = { xyz[0], xyz[1], xyz[2] };
        if (forward) {
            x[slab.dim] += slab.delta;
            if (x[slab.dim] >= pbcd)
                x[slab.dim] -= pbcd;
        }
        else {
            x[slab.dim] -= slab.delta;
            if (x[slab.dim] < 0)
                x[slab.dim] += pbcd;
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
        double exact_w = forward ?
                    dist_forward(atom->xyz[slab.dim], xyz[slab.dim], pbcd)
                  : dist_forward(xyz[slab.dim], atom->xyz[slab.dim], pbcd);
        avg_w.add_x(exact_w);
    }
    if (verbosity > 0)
        printf("ok (t_%c = %g +/- %g)\n", 'x'+slab.dim, avg_w.mean(),
                                                        avg_w.stddev());
    return avg_w.mean();
}

DeltaLim parse_delta_lim(gengetopt_args_info const& args)
{
    DeltaLim lim;
    lim.lo = args.min_delta_given ? args.min_delta_arg : 0.;
    lim.hi = args.max_delta_given ? args.max_delta_arg : 1e9;
    return lim;
}

double search_for_translation(dbr_aconf const& aconf,
                              CellMethod const& cm,
                              gengetopt_args_info const& args,
                              int dim, double val)
{
    int a0 = select_atom(aconf, dim, val);
    if (verbosity > 0)
        printf("atom0: %s\n", atom_to_str(aconf, a0));
    DeltaLim lim = parse_delta_lim(args);

    double eps = args.epsilon_arg;
    bool single_trans = !args.periodic_given;
    // possible lengths of translation vectors, based on one atom's neighbors
    vector<double> tr = nabe_in_direction_distances(aconf, a0, dim, lim, eps);
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

// get distances from atom atom_nr to its images that are in distance d
// lim.lo < d < lim.hi, but only if all atoms in the slab (or in the whole
// system) have images in the same relative position.
vector<Trans> get_translations(CellMethod const& cm,
                               int atom_nr,
                               SlabInPBC const& pslab,
                               DeltaLim const& lim,
                               double epsilon)
{
    vector<Trans> ret;
    dbr_atom const& atom0 = cm.aconf().atoms[atom_nr];
    dbr_real const* xyz0 = atom0.xyz;
    for (int i = 0; i != cm.aconf().n; ++i) {
        if (i == atom_nr || strcmp(atom0.name, cm.aconf().atoms[i].name) != 0)
            continue;
        double dist = cm.get_dist(xyz0, cm.aconf().atoms[i].xyz);
        if (dist < lim.lo || dist > lim.hi)
            continue;
        dbr_real r[3];
        //dbr_diff3(cm.aconf().atoms[i].xyz, xyz0, r);
        cm.get_shortest_vec(xyz0, cm.aconf().atoms[i].xyz, r);
        // check for the translation symmetry and calculate average
        bool ok = true;
        StdDev sd[3];
        for (int j = 0; j != cm.aconf().n; ++j) {
            if (pslab.ok() && !pslab.has(cm.aconf().atoms[j].xyz))
                continue;
            dbr_real delta[3];
            dbr_atom const* img = cm.get_image(r, j, epsilon, delta);
            if (img == NULL) {
                ok = false;
                break;
            }
            for (int k = 0; k < 3; ++k)
                sd[k].add_x(r[k] + delta[k]);
        }
        if (ok) {
            Trans trans;
            for (int k = 0; k < 3; ++k) {
                trans.r[k] = sd[k].mean();
                trans.stddev[k] = sd[k].stddev();
            }
            trans.dist = dbr_len3(sd[0].mean(), sd[1].mean(), sd[2].mean());
            ret.push_back(trans);
            //printf("debug: %d r=(%f %f %f)\n", i, r[0], r[1], r[2]);
        }
    }
    return ret;
}

// search for translational symmetry in all directions
vector<Trans> find_trans_sym(CellMethod const& cm, SlabInPBC const& pslab,
                             double epsilon, DeltaLim const& lim)
{
    // print min. interatomic distance
    double min_dist = cm.find_min_distance();
    if (epsilon >= min_dist / 2.) {
        printf("WARNING: Min. interatomic distance is: %g.\n", min_dist);
        printf("WARNING: Epsilon should be smaller than %g.\n", min_dist/2);
    }

    int atom0_nr = get_any_atom_in_slab(cm.aconf(), pslab);

    vector<Trans> tt = get_translations(cm, atom0_nr, pslab, lim, epsilon);

    if (tt.empty()) {
        printf("Translation symmetry not found in the %s.\n",
                pslab.ok() ? "slab" : "system");
    }
    return tt;
}

// returns (sorted) translation vectors T that are not T=n*T'
vector<Trans> get_basic_trans(CellMethod const& cm, SlabInPBC const& pslab,
                              gengetopt_args_info const& args)
{
    vector<Trans> basic;
    vector<Trans> tt = find_trans_sym(cm, pslab, args.epsilon_arg,
                                      parse_delta_lim(args));

    //// make the first index non-negative (if x < 0: x=-x, y=-y, z=-z)
    //for (vector<Trans>::iterator i = tt.begin(); i != tt.end(); ++i)
    //    if (i->r[0] < 0)
    //        for (int k = 0; k < 3; ++k)
    //            i->r[k] = - i->r[k];

    sort(tt.begin(), tt.end());

    // flag duplicates by setting dist=-1
    double eps = args.epsilon_arg;
    for (vector<Trans>::iterator i = tt.begin(); i != tt.end(); ++i) {

        // check if it's not n * another translation
        for (vector<Trans>::const_iterator j = tt.begin(); j != i; ++j) {
            if (j->dist < 0)
                continue;
            // calculate average multiplier
            double m = i->dist / j->dist;
            // check if this is a multiplier for all coordinates
            if (fabs(m * j->r[0] - i->r[0]) < eps
                    && fabs(m * j->r[1] - i->r[1]) < eps
                    && fabs(m * j->r[2] - i->r[2]) < eps) {
                i->dist = -1; // mark duplicate
                break;
            }
        }
    }
    for (vector<Trans>::const_iterator i = tt.begin(); i != tt.end(); ++i)
        if (i->dist > 0) // i->dist is also a flag for duplicates
            basic.push_back(*i);
    return basic;
}

void print_trans_sym(dbr_aconf const& aconf, Slab const& slab,
                     gengetopt_args_info const& args)
{
    SlabInPBC pslab(aconf, slab);
    CellMethod cm(aconf, args.min_cell_arg);
    vector<Trans> tt = get_basic_trans(cm, pslab, args);

    for (vector<Trans>::const_iterator i = tt.begin(); i != tt.end(); ++i) {
        printf("T = (% f % f % f) +/- (%.2g %.2g %.2g)  |T|=%g\n",
               i->r[0], i->r[1], i->r[2],
               i->stddev[0], i->stddev[1], i->stddev[2],
               i->dist);
    }
}

// returns true if the angle is 90 +- 3 degrees
bool is_orthogonal(vector<Trans>::const_iterator u,
                   vector<Trans>::const_iterator v)
{
    double dot = u->r[0] * v->r[0] + u->r[1] * v->r[1] + u->r[2] * v->r[2];
    double cos_angle = dot / (u->dist * v->dist);
    return fabs(cos_angle) < cos(87*M_PI/180.);
}

double find_vector_for_cubic_cell(vector<Trans> const& tt, double eps,
                                  const double* a[3])
{
    for (vector<Trans>::const_iterator i = tt.begin(); i != tt.end(); ) {
        vector<Trans>::const_iterator begin = i;
        for (i = begin; i != tt.end() && fabs(i->dist - begin->dist) < eps; ++i)
            ;
        if (verbosity > 1)
            printf("checking for cell parameter a=%g +- %g (%d vectors)\n",
                   begin->dist, eps, (int) (i - begin));
        if (i - begin < 3)
            continue;
        // find 3 orthogonal vectors
        for (vector<Trans>::const_iterator u = begin; u != i; ++u)
            for (vector<Trans>::const_iterator v = u+1; v != i; ++v) {
                if (!is_orthogonal(u, v))
                    continue;
                for (vector<Trans>::const_iterator w = v+1; w != i; ++w)
                    if (is_orthogonal(u, w) && is_orthogonal(v, w)) {
                        a[0] = u->r;
                        a[1] = v->r;
                        a[2] = w->r;
                        return (u->dist + v->dist + w->dist) / 3.;
                    }
            }
    }
    return 0.;
}

void change_order_and_sign_of_vectors(const double** v, double a[][3])
{
    int ax_idx = 0;
    for (int i = 1; i != 3; ++i)
        if (fabs(v[i][0]) > fabs(v[ax_idx][0]))
            ax_idx = i;
    int ay_idx = (ax_idx + 1) % 3;
    int az_idx = (ax_idx + 2) % 3;
    if (fabs(v[az_idx][1]) > fabs(v[ay_idx][1]))
        swap(ay_idx, az_idx);
    for (int i = 0; i != 3; ++i) {
        a[0][i] = v[ax_idx][i];
        a[1][i] = v[ay_idx][i];
        a[2][i] = v[az_idx][i];
    }

    for (int i = 0; i != 3; ++i)
        if (a[i][0] < 0)
            for (int j = 0; j != 3; ++j)
                a[i][j] = -a[i][j];
}

struct AtomNewPos
{
    int atom_idx;
    double xyz[3]; // in not normalized coordinates
};

void make_cubic(dbr_aconf& aconf, Slab const& slab,
                gengetopt_args_info const& args)
{
    if (aconf.reduced_coordinates) {
        printf("Reduced coordinates cannot be used with --make-cubic option\n");
        exit(EXIT_FAILURE);
    }
    SlabInPBC pslab(aconf, slab);
    double eps = args.epsilon_arg;
    double cell_param = 0.;
    double a[3][3];

    if (args.make_cubic_arg == NULL)
    {
        CellMethod cm(aconf, args.min_cell_arg);
        vector<Trans> tt = get_basic_trans(cm, pslab, args);
        const double* v[3];
        cell_param = find_vector_for_cubic_cell(tt, eps, v);
        if (cell_param == 0.) {
            printf("Cubic cell not found.\n");
            return;
        }

        change_order_and_sign_of_vectors(v, a);
    }
    else {
        bool r = parse_n_numbers(args.make_cubic_arg, 9, (double*) a);
        if (!r) {
            printf("Wrong argument to --make-cubic option.\n");
            return;
        }
        cell_param = (dbr_len3v(a[0]) + dbr_len3v(a[1]) + dbr_len3v(a[2])) / 3.;
    }

    if (verbosity > 1) {
        printf("Making cubic cell with a=%g and vectors:\n", cell_param);
        for (int i = 0; i != 3; ++i)
            printf("\t(%.9g, %.9g, %.9g)\n", a[i][0], a[i][1], a[i][2]);
    }

    // pick a reference atom
    const dbr_real* ref_xyz = NULL;
    for (int i = 0; i != aconf.n; ++i) {
        if (!pslab.ok() || pslab.has(aconf.atoms[i].xyz)) {
            ref_xyz = aconf.atoms[i].xyz;
            break;
        }
    }
    assert (ref_xyz != NULL);

    // find sites of atoms in the unit cell
    double b[3][3];
    inverse_3x3_matrix(a, b);

    vector<AtomNewPos> new_positions;
    // shift_sd is used to calculate a correction to the reference atom
    // position, that minimizes shifts of atoms
    StdDev shift_sd[3];
    for (int i = 0; i != aconf.n; ++i) {
        if (pslab.ok() && !pslab.has(aconf.atoms[i].xyz))
            continue;
        AtomNewPos newpos;
        newpos.atom_idx = i;
        const dbr_real *ax = aconf.atoms[i].xyz;
        const dbr_real xyz[3] = { ax[0] - ref_xyz[0],
                                  ax[1] - ref_xyz[1],
                                  ax[2] - ref_xyz[2] };
        // r = xyz . b
        double r[3];
        for (int j = 0; j != 3; ++j)
            r[j] = b[0][j] * xyz[0] + b[1][j] * xyz[1] + b[2][j] * xyz[2];

        double normalized_pos[3];
        for (int j = 0; j != 3; ++j) {
            normalized_pos[j] = floor(r[j] * 4 + 0.5) / 4.;
            double diff = (normalized_pos[j] - r[j]) * cell_param;
            if (fabs(diff) > eps) {
                printf("Atom %d in unexpected site (%g %g %g) -> "
                        "(%.2f %.2f %.2f)\n",
                        i, ax[0], ax[1], ax[2], r[0], r[1], r[2]);
                return;
            }
            shift_sd[j].add_x(diff);
        }
        // newpos.xyz = ref_xyz + normalized_pos . a
        for (int j = 0; j != 3; ++j)
            newpos.xyz[j] = ref_xyz[j] + a[0][j] * normalized_pos[0] +
                a[1][j] * normalized_pos[1] + a[2][j] * normalized_pos[2];

        new_positions.push_back(newpos);
        //// it doesn't include shift_sd:
        //printf("Atom %d: (%g %g %g) -> (%.2f %.2f %.2f) -> (%g %g %g)\n",
        //        i, ax[0], ax[1], ax[2], r[0], r[1], r[2],
        //        newpos.xyz[0], newpos.xyz[1], newpos.xyz[2]);
    }
    if (verbosity > 1)
      printf("Changing %d atoms, std. dev. of distortions: %.3g, %.3g, %.3g\n",
             (int) new_positions.size(),
             shift_sd[0].stddev(), shift_sd[1].stddev(), shift_sd[2].stddev());

    // change atom coordinates
    for (vector<AtomNewPos>::const_iterator i = new_positions.begin();
                                               i != new_positions.end(); ++i)
        for (int j = 0; j != 3; ++j)
            aconf.atoms[i->atom_idx].xyz[j] = i->xyz[j] - shift_sd[j].mean();
}

// If output is false:
//   delete atoms with `dim' coord in <x0, x0+delta). 
//   If delta < 0, deleted range is <x0-|delta|, x0).
// Otherwise delete the rest of atoms.
void delete_atoms(dbr_aconf& aconf, Slab const& slab, bool del_outside)
{
    //printf("slab: %c=%g, delta: %g\n", 'x'+slab.dim, slab.x0, slab.delta);
    int n = 0;
    SlabInPBC pslab(aconf, slab);
    for (int i = 0; i != aconf.n; ++i) {
        bool atom_in_slab = pslab.has(aconf.atoms[i].xyz);
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

void add_vacuum(dbr_aconf& aconf, int dim, double x0, double delta)
{
    for (int i = 0; i != aconf.n; ++i)
        if (aconf.atoms[i].xyz[dim] > x0)
            aconf.atoms[i].xyz[dim] += delta;
    expand_pbc(aconf, dim, delta);
}

void resize_atoms_table(dbr_aconf& aconf, int new_size)
{
    dbr_atom* new_atoms = new dbr_atom[new_size];
    memcpy(new_atoms, aconf.atoms, sizeof(dbr_atoms) * min(aconf.n, new_size));
    delete [] aconf.atoms;
    aconf.atoms = new_atoms;
    if (!aconf.auxiliary.empty())
        aconf.auxiliary.resize(new_size);
}

void find_translation_vector(CellMethod const& cm,
                             SlabInPBC const& pslab,
                             double epsilon,
                             DeltaLim const& lim,
                             int tr_sign,
                             double *r)
{
    vector<Trans> tt = find_trans_sym(cm, pslab, epsilon, lim);
    if (tt.empty()) {
        exit(EXIT_FAILURE);
    }
    // we choose the vector in direction closest to slab normal
    Trans const* t = NULL;
    for (vector<Trans>::const_iterator i = tt.begin(); i != tt.end(); ++i) {
        double norm_r = fabs(i->r[pslab.dim()]);
        if (norm_r < 1e-9)
            continue;
        if (t == NULL || norm_r / i->dist > fabs(t->r[pslab.dim()]) / t->dist)
            t = &(*i);
    }
    if (t == NULL) {
        printf("Can not find translation vector.\n");
        exit(EXIT_FAILURE);
    }
    // in the case of surface, translation is possible only in toward
    // the system, and we want translation (to extend the system)
    // toward the vacuum.
    int sign = t->r[pslab.dim()] > 0 ? 1 : -1;
    for (int i = 0; i < 3; ++i)
        r[i] = (sign == tr_sign ? t->r[i] : -t->r[i]);

    if (t && verbosity > 0)
        printf("T = (% f % f % f) +/- (%g %g %g)\n",
               r[0], r[1], r[2],
               t->stddev[0], t->stddev[1], t->stddev[2]);
}

// bookmark atoms that are to be copied
vector<int> bookmark_slab_atoms(dbr_aconf const& aconf, Slab const& slab,
                                CellMethod const& cm, double epsilon)
{
    vector<int> orig;
    double pbcd = get_pbc(aconf, slab.dim);
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
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
    return orig;
}

vector<int> bookmark_one_shift_group(dbr_aconf const& aconf, double const* r,
                                     int dim, double x0,
                                     CellMethod const& cm, double epsilon)
{
    vector<int> orig;
    double pbcd = get_pbc(aconf, dim);
    x0 -= floor(x0 / pbcd) * pbcd;
    double rd = r[dim];
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real* xyz = aconf.atoms[i].xyz;
        double dist = (rd >= 0 ? dist_forward(x0, xyz[dim], pbcd)
                               : dist_forward(xyz[dim], x0, pbcd));
        if (dist < fabs(rd) - epsilon)
            orig.push_back(i);
        // we must ensure that of every pair original-image in this area,
        // exactly one atom is bookmarked
        else if (dist < fabs(rd) + epsilon) {
            dbr_real x[3] = { xyz[0]+r[0], xyz[1]+r[1], xyz[2]+r[2] };
            dbr_atom const* img = cm.get_atom_at(x, epsilon);
            if (!img)
                continue;
            double di = (rd >= 0 ? dist_forward(x0, img->xyz[dim], pbcd)
                                 : dist_forward(img->xyz[dim], x0, pbcd));
            if (di > fabs(rd))
                orig.push_back(i);
        }
    }
    return orig;
}

double calc_avg_coordinate(dbr_aconf const& aconf, vector<int> const& sel,
                           int dim)
{
    StdDev sd;
    double pbcd = get_pbc(aconf, dim);
    for (vector<int>::const_iterator i = sel.begin(); i != sel.end(); ++i) {
        double t = aconf.atoms[*i].xyz[dim];
        if (sd.n() != 0 && fabs(t - sd.mean()) > pbcd/2.) {
            if (t < sd.mean())
                t += pbcd;
            else
                t -= pbcd;
        }
        sd.add_x(t);
    }
    //printf("DEBUG: %s\n", sd.str().c_str());
    return sd.mean();
}

void print_density(dbr_aconf const& aconf, Slab const& slab,
                   CellMethod const& cm, gengetopt_args_info const& args)
{
    if (aconf.reduced_coordinates) {
        printf("Reduced coordinates can not be used with --density option\n");
        exit(EXIT_FAILURE);
    }
    int dim = slab.dim;
    SlabInPBC pslab(aconf, slab);
    double x0 = pslab.x0();
    double delta = pslab.delta();
    double pbc_area = cm.pbc()[(dim+1) % 3] * cm.pbc()[(dim+2) % 3];

    if (verbosity > -1) {
        int atom_count = 0;
        for (int i = 0; i != aconf.n; ++i)
            if (pslab.has(aconf.atoms[i].xyz))
                ++atom_count;
        double vol = pbc_area * delta;
        printf("N / V = %d / %g = %g\n", atom_count, vol, atom_count / vol);
    }

    double epsilon = args.epsilon_arg;
    double d = 1.0;
    DeltaLim lim = parse_delta_lim(args);
    if (!args.max_delta_given)
        lim.hi = 4*d;
    vector<bool> margin_atoms(aconf.n, false);
    double r[3];

    Slab right = { dim, x0 + delta, -d };
    find_translation_vector(cm, SlabInPBC(aconf, right), epsilon, lim, 1, r);
    vector<int> bright =
                bookmark_one_shift_group(aconf, r, dim, x0+delta, cm, epsilon);
    double r_avg = calc_avg_coordinate(aconf, bright, dim);
    for (vector<int>::const_iterator i = bright.begin(); i != bright.end(); ++i)
        margin_atoms[*i] = true;
    double vol_marginr = pbc_area * fabs(r[dim]);

    Slab left = { dim, x0, d };
    find_translation_vector(cm, SlabInPBC(aconf, left), epsilon, lim, -1, r);
    vector<int> bleft =
                bookmark_one_shift_group(aconf, r, dim, x0, cm, epsilon);
    double l_avg = calc_avg_coordinate(aconf, bleft, dim);
    for (vector<int>::const_iterator i = bleft.begin(); i != bleft.end(); ++i)
        margin_atoms[*i] = true;
    double vol_marginl = pbc_area * fabs(r[dim]);

    double pbcd = get_pbc(aconf, dim);
    double center_to_center = dist_forward(r_avg, l_avg, pbcd);
    if (verbosity > 0) {
        printf("-- Margins of sizes defined by translation vectors T are used"
                " to --\n");
        printf("-- make density independent from small changes of slab size."
                "     --\n");
        printf("right margin at %g has %d atoms\n", x0+delta,
                                                    (int) bright.size());
        printf("left margin at %g has %d atoms\n", x0, (int) bleft.size());
        printf("margin centers: %g, %g, delta=%g\n", r_avg, l_avg,
                                                     center_to_center);
        printf("margins density: %g, %g\n", bright.size() / vol_marginr,
                                            bleft.size() / vol_marginl);
    }

    Slab center = { dim, x0 + d/2., delta - d };
    SlabInPBC pcenter(aconf, center);
    int center_counter = 0;
    for (int i = 0; i != aconf.n; ++i)
        if (pcenter.has(aconf.atoms[i].xyz) && !margin_atoms[i])
            ++center_counter;
    double count = center_counter + (bleft.size() + bright.size()) / 2.;
    double vol = pbc_area * center_to_center;
    printf("numeric density: %g\n", count / vol);
#if 0
    double reference_density = 8 / pow(4.321035, 3); // at. / A^3
    double ref_vol = count / reference_density;
    double extra_vol = vol - ref_vol;
    double vacuum_width = extra_vol / pbc_area;
    printf("equivalent of vacuum width: %g\n", vacuum_width);
#endif
}

void add_copy(dbr_aconf& aconf, Slab const& slab,
              CellMethod const& cm, double epsilon, int N)
{
    vector<int> orig = bookmark_slab_atoms(aconf, slab, cm, epsilon);
    if (!orig.empty())
        resize_atoms_table(aconf, aconf.n + N * orig.size());

    double delta = slab.delta;
    for (int i = 0; i < N; ++i) {
        add_vacuum(aconf, slab.dim, slab.x0, slab.delta);

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

void add_by_translation(dbr_aconf& aconf, Slab const& slab,
                        CellMethod const& cm,
                        gengetopt_args_info const& args)
{
    double extra_width = args.add_arg;
    if (extra_width <= 0) {
        printf("WARNING: WIDTH in --add=WIDTH should be positive.\n");
    }
    double eps = args.epsilon_arg;
    double r[3];
    find_translation_vector(cm,
                            SlabInPBC(aconf, slab),
                            eps, parse_delta_lim(args), 0, r);

    int orig_n = aconf.n;
    // fill the added space with atoms
    double pbcd = get_pbc(aconf, slab.dim);
    // we will also fill the margins (of width=eps) around the added space
    double rd = r[slab.dim];
    // first count the new atoms
    vector<pair<int, int> > cc;// pairs (original atom index, number of copies)
    int sum = 0;
    for (int i = 0; i != aconf.n; ++i) {
        dbr_real x = aconf.atoms[i].xyz[slab.dim];
        double dist_to_x0 = (rd >= 0 ? dist_forward(slab.x0, x, pbcd)
                                     : dist_forward(x, slab.x0, pbcd));
        if (dist_to_x0 > fabs(rd) + eps)
            continue;
        int n = (int) ((dist_to_x0 + extra_width + 2*eps) / fabs(rd));
        if (n > 0) {
            cc.push_back(make_pair(i, n));
            sum += n;
        }
    }
    if (verbosity > 1)
        printf("%d total copies of %d original atoms are considered.\n",
               sum, (int) cc.size());

    // if option -b is given, the outer margin is not filled
    if (args.bound_given) {
        int n_rej = 0;
        for (vector<pair<int, int> >::iterator i = cc.begin();
                i != cc.end(); ++i) {
            double new_x = aconf.atoms[i->first].xyz[slab.dim] + i->second * rd;
            double dist = (rd >= 0 ? dist_forward(new_x, slab.x0, pbcd)
                                   : dist_forward(slab.x0, new_x, pbcd));
            if (dist > extra_width) {
                --(i->second);
                ++n_rej;
            }
        }
        sum -= n_rej;
        if (verbosity > 1)
            printf("%d copies rejected: in distance > %g from %c=%g.\n",
                   n_rej, extra_width, 'x'+slab.dim, slab.x0);
    }

    // resize PBC and shift atoms
    add_vacuum(aconf, slab.dim, slab.x0, fabs(extra_width));

    // add new atoms to aconf
    resize_atoms_table(aconf, aconf.n + sum);
    int pos = aconf.n;
    aconf.n += sum;
    for (vector<pair<int, int> >::const_iterator i = cc.begin();
            i != cc.end(); ++i)
        for (int j = 0; j < i->second; ++j) {
            dbr_copy_atom(aconf, i->first, pos);
            for (int k = 0; k < 3; ++k)
                aconf.atoms[pos].xyz[k] += (j+1) * r[k];
            ++pos;
    }
    assert (pos == aconf.n);

    // remove duplicates
    wrap_to_pbc(aconf);
    int n_dups = merge_atoms(aconf, args, true);
    if (verbosity > 1)
        printf("%d copies rejected: duplicates.\n", n_dups);

    if (verbosity > -1)
        printf("%d atoms were added (%d -> %d).\n", aconf.n - orig_n,
                                                    orig_n, aconf.n);
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
    if (!aconf.auxiliary.empty()) {
        aconf.auxiliary.resize(new_n);
        for (int img = 1; img < x*y*z; ++img)
            for (int i = 0; i != aconf.n; ++i)
                aconf.auxiliary[img * aconf.n + i] = aconf.auxiliary[i];
    }
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
    double d[3];
    bool r = parse_n_numbers(arg, 3, d);
    if (!r)
        return false;

    // shift atoms
    for (int i = 0; i != aconf.n; ++i)
        for (int j = 0; j < 3; ++j)
            aconf.atoms[i].xyz[j] += d[j];
    wrap_to_pbc(aconf);
    return true;
}

bool resize_system(dbr_aconf& aconf, const char* arg)
{
    double d[3];
    bool r = parse_n_numbers(arg, 3, d);
    // if arg is not x,y,z it must be a filename
    if (!r) {
        dbr_aconf ref = read_atoms_from_file(arg, false, "");
        for (int i = 0; i != 3; ++i)
            d[i] = *get_pbc_ptr(ref, i);
    }

    for (int i = 0; i != 3; ++i) {
        if (d[i] <= 0)
            continue;
        double* pp = get_pbc_ptr(aconf, i);
        if (!aconf.reduced_coordinates)
            for (int j = 0; j != aconf.n; ++j)
                aconf.atoms[j].xyz[i] *= (d[i] / *pp);
        *pp = d[i];
    }
    return true;
}

// finds the largest span of vacuum in direction dim and returns one of the 
// boundaries of the vacuum.
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
            && !args.find_trans_given && !args.make_cubic_given
            && !args.merge_given && !args.t1_given && !args.resize_given) {
        fprintf(stderr, "One of -x, -y, -z, -b, -N, -S, -F options must be "
                        "given (-h will show all valid options).\n");
        return EXIT_FAILURE;
    }
    if (args.in_place_given && args.output_given) {
        fprintf(stderr, "Only one of -i and -o options can be given "
                "(-h will show explaination).\n");
        return EXIT_FAILURE;
    }
    if (args.add_vacuum_given + args.add_copy_given + args.add_given > 1) {
        fprintf(stderr, "Only one of add* options can be given."
                "can be given.\n");
        return EXIT_FAILURE;
    }

    if (args.inputs_num != 1) {
        fprintf(stderr, "Exactly one input file should be given.\n");
        return EXIT_FAILURE;
    }

    // reading configuration file
    bool reduced_coords = args.reduced_given;
    dbr_aconf aconf = read_atoms_from_file(args.inputs[0], reduced_coords, "");

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
        if (slab.x0 < 0)
            slab.x0 += pbcd;
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
        if (verbosity > 1) {
            double min_dist = cm.find_min_distance();
            printf("Min. interatomic distance: %g\n", min_dist);
        }

        // if epsilon is larger than half of cell, increase cell sizes
        if (args.epsilon_arg >= cm.min_cell_size() / 2.) {
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
        if (args.density_given)
            print_density(aconf, slab, cm, args);
        if (args.delete_given)
            delete_atoms(aconf, slab, /*del_outside=*/false);
        if (args.add_vacuum_given)
            add_vacuum(aconf, slab.dim, slab.x0, slab.delta);
        if (args.add_copy_given)
            add_copy(aconf, slab, cm, args.epsilon_arg, args.add_copy_arg);
        if (args.add_given)
            add_by_translation(aconf, slab, cm, args);
        if (args.cut_given) {
            delete_atoms(aconf, slab, /*del_outside=*/true);
            expand_pbc(aconf, slab.dim, slab.delta-pbcd);
        }
    }

    if (args.multiply_given) {
        bool r = multiply_pbc(aconf, args.multiply_arg);
        if (!r) {
            fprintf(stderr,
                    "error parsing -N arguments, expected e.g. -N2x1x1\n");
            return EXIT_FAILURE;
        }
    }

    if (args.shift_given)
        shift_under_pbc(aconf, args.shift_arg);

    if (args.resize_given)
        resize_system(aconf, args.resize_arg);

    if (args.find_trans_given)
        print_trans_sym(aconf, slab, args);

    if (args.make_cubic_given)
        make_cubic(aconf, slab, args);

    if (args.t1_given)
        transform1(aconf, args);

    if (args.merge_given)
        merge_atoms(aconf, args, false);

    aconf.comments.insert(aconf.comments.begin(), argv_as_str(argc, argv));
    if (args.in_place_given)
        write_file_with_atoms(aconf, aconf.orig_filename, "");
    else if (args.output_given)
        write_file_with_atoms(aconf, args.output_arg, "");

    return EXIT_SUCCESS;
}

