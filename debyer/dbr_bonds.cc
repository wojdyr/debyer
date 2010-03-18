// dbr_bonds -- coordination number, bond length and angles between bonds

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
//
// $Id$

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cfloat>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric> // accumulate
#include "fileio.h"
#include "cells.h"
#include "utils.h"

using namespace std;


int main(int argc, char **argv)
{
    if (argc != 2 && argc != 3) {
        printf("Usage: dbr_bonds filename max-bond-length\n");
        return -1;
    }
    bool reduced_coords = false;
    dbr_aconf aconf = read_atoms_from_file(argv[1], reduced_coords);


    // check PBC
    if (aconf.pbc.v01 != 0 || aconf.pbc.v02 != 0 || aconf.pbc.v10 != 0 ||
            aconf.pbc.v12 != 0 || aconf.pbc.v20 != 0 || aconf.pbc.v21 != 0) {
        fprintf(stderr, "Non-orthorhombic systems are not handled.\n");
        return EXIT_FAILURE;
    }
    printf("system size: %g x %g x %g\n", aconf.pbc.v00, aconf.pbc.v11,
                                                            aconf.pbc.v22);
    // atom statistics
    printf("atom count: %d total", aconf.n);
    fflush(stdout);
    vector<int> species(aconf.n, -1);
    vector<string> symbols;
    for (int i = 0; i != aconf.n; ++i) {
        string name(aconf.atoms[i].name);
        vector<string>::const_iterator it
            = find(symbols.begin(), symbols.end(), name);
        species[i] = it - symbols.begin();
        if (it == symbols.end())
            symbols.push_back(name);
    }
    int n_species = symbols.size();
    for (int i = 0; i != n_species; ++i) {
        int n = count(species.begin(), species.end(), i);;
        printf(", %d %s", n, symbols[i].c_str());
    }
    printf("\n");
    fflush(stdout);

    if (argc == 2) {
        printf("To show info about bonds, specify max. bond length.\n");
        return 0;
    }

    char *endptr;
    double rcut = strtod(argv[2], &endptr);
    if (*endptr != 0 || rcut <= 0) {
        fprintf(stderr, "Wrong value for max_bondlength\n");
        return EXIT_FAILURE;
    }

    // put atoms into cells and find bonds
    CellMethod cm(aconf, rcut); //args.min_cell_arg);
    vector<Bond> bonds;
    cm.find_all_bonds(rcut, bonds);
    if (bonds.empty()) {
        fprintf(stderr, "No bonds were found.\n");
        return EXIT_FAILURE;
    }

    // coordination number statistics (total and per species)
    vector<int> cn(aconf.n, 0);
    for (vector<Bond>::const_iterator i = bonds.begin(); i != bonds.end(); ++i){
        ++cn[i->a];
        ++cn[i->b];
    }
    int max_cn = *max_element(cn.begin(), cn.end());
    vector<vector<int> > cn_stats(n_species, vector<int>(max_cn+1, 0));
    vector<int> cn_tot_stats(max_cn+1, 0);
    printf("\ncoordination number\n");
    printf("  total:");
    for (int i = 0; i != aconf.n; ++i) {
        ++cn_stats[species[i]][cn[i]];
        ++cn_tot_stats[cn[i]];
    }
    for (int i = 0; i != max_cn+1; ++i)
        if (cn_tot_stats[i] != 0)
            printf(" %6d CN=%d", cn_tot_stats[i], i);
    for (int sp = 0; sp != n_species; ++sp) {
        printf("\n  %-5s:", symbols[sp].c_str());
        for (int i = 0; i != max_cn+1; ++i)
            if (cn_tot_stats[i] != 0)
                printf(" %6d CN=%d", cn_stats[sp][i], i);
    }
    printf("\n");
    fflush(stdout);

    // bond statistics
    double min_bond = DBL_MAX, max_bond = 0;
    for (vector<Bond>::const_iterator i = bonds.begin(); i != bonds.end(); ++i){
        if (i->length < min_bond)
            min_bond = i->length;
        if (i->length > max_bond)
            max_bond = i->length;
    }
    printf("\nmin. bond length: %g", min_bond);
    printf("\nmax. bond length: %g (%g was the limit)", max_bond, rcut);
    printf("\nbond count: %d total", (int) bonds.size());
    fflush(stdout);

    printf("\nbond statistics: # bonds, avg bond length, std. dev.\n");
    vector<vector<StdDev> > bond_stats(n_species, vector<StdDev>(n_species));
    for (vector<Bond>::const_iterator i = bonds.begin(); i != bonds.end(); ++i){
        int sp1 = species[i->a];
        int sp2 = species[i->b];
        if (sp1 > sp2)
            swap(sp1, sp2);
        bond_stats[sp1][sp2].add_x(i->length);
    }

    for (int i = 0; i != n_species; ++i)
        for (int j = i; j != n_species; ++j) {
            const StdDev &sd = bond_stats[i][j];
            printf("b.s. %2s-%-2s %8d  %.3f  %.3g\n", symbols[i].c_str(),
                           symbols[j].c_str(), sd.n(), sd.mean(), sd.stddev());
        }
    // make the matrix symmetric, just for easier handling
    for (int i = 0; i != n_species; ++i)
        for (int j = 0; j != i; ++j)
            bond_stats[i][j] = bond_stats[j][i];

#if 0
    // GB statistics
    double gb_area = aconf.pbc.v00 * aconf.pbc.v11;
    printf("\nGB area: %g", gb_area);
    int si = find(symbols.begin(), symbols.end(), "Si") - symbols.begin();
    int c = find(symbols.begin(), symbols.end(), "C") - symbols.begin();
    printf("\nGB: Si-Si:%g  C-C:%g  Si<4:%g  C<4:%g  Si>4:%g  C>4:%g\n",
           bond_stats[si][si].n() / gb_area,
           bond_stats[c][c].n() / gb_area,
           accumulate(cn_stats[si].begin(), cn_stats[si].begin()+4, 0) /gb_area,
           accumulate(cn_stats[c].begin(), cn_stats[c].begin()+4, 0) / gb_area,
           accumulate(cn_stats[si].begin()+5, cn_stats[si].end(), 0) / gb_area,
           accumulate(cn_stats[c].begin()+5, cn_stats[c].end(), 0) / gb_area);
#endif

    // bond angle statistics
    //  list neighbours for each atom
    vector<vector<int> > nabes(aconf.n);
    for (vector<vector<int> >::iterator i = nabes.begin(); i!=nabes.end(); ++i)
        i->reserve(4);
    for (vector<Bond>::const_iterator i = bonds.begin(); i != bonds.end(); ++i){
        nabes[i->a].push_back(i->b);
        nabes[i->b].push_back(i->a);
    }
    //  calculate angles
    StdDev cn4_angles;
    //vector<int> cn4_histogram(180, 0);
    StdDev all_angles;
    for (size_t i = 0; i != nabes.size(); ++i) {
        const vector<int>& nn = nabes[i];
        const dbr_real* xyz = aconf.atoms[i].xyz;
        string name(aconf.atoms[i].name);
        if (name != "Si" && name != "C")
            continue;
        for (size_t ileft = 0; ileft != nn.size(); ++ileft) {
            const dbr_real* left = aconf.atoms[nn[ileft]].xyz;
            for (size_t iright = ileft+1; iright != nn.size(); ++iright) {
                const dbr_real* right = aconf.atoms[nn[iright]].xyz;
                double angle = cm.get_angle(left, xyz, right) * 180. / M_PI;
                assert (angle >= 0 && angle <= 180.);
                all_angles.add_x(angle);
                if (nn.size() == 4) {
                    cn4_angles.add_x(angle);
                    //++cn4_histogram[angle];
                }
            }
        }
    }
    printf("\nbond angle statistics: # angles, avg angle, std. dev.\n");
    printf("b.a.s. for cn=4  %8d  %.2f   %.3g\n",
            cn4_angles.n(), cn4_angles.mean(), cn4_angles.stddev());
    printf("b.a.s. for all   %8d  %.2f   %.3g\n",
            all_angles.n(), all_angles.mean(), all_angles.stddev());
    //for (size_t i = 0; i != cn4_histogram.size(); ++i)
    //    printf("%g: %d\n", 1.0 * i, cn4_histogram[i]);
}


