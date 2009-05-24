// Tool to calculate diffusion and make .aux files with displacement
// for a sequance of .cfg files
// .cfg and .aux files are used by AtomEye visualization program

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <cstring>

#include "fileio.h"
#include "utils.h"

using namespace std;

// define SPLIT_DIFFUSION to calculate diffusion separately for multiple
// categories of atoms.
// Category for atoms is taken from file atoms.cat. Each line in this file
// corresponds to one atom, and should contain integer number from 0 to
// (SPLIT_DIFFUSION - 1).

#define SPLIT_DIFFUSION 2

void write_displacement_aux(ofstream& f, int n,
                            xyz_name const* c1, xyz_name const* c2,
                            dbr_pbc const& pbc,
                            vector<string> const& old_lines)
{
    StdDev diffus;
#ifdef SPLIT_DIFFUSION
    vector<StdDev> cat_diffus(SPLIT_DIFFUSION);
    vector<int> cat;
    const char* cat_fn = "atoms.cat";
    ifstream split_file(cat_fn);
    if (split_file) {
        string s;
        while (getline(split_file, s)) {
            int c = strtol(s.c_str(), 0, 10);
            cat.push_back(c);
        }
        assert(cat.size() == (size_t) n);
    }
    else {
        cerr << "WARNING: Can not open file: " << cat_fn << endl;
    }
#endif
    for (int i = 0; i < n; ++i)  {
        dbr_xyz d;
        for (int j = 0; j < 3; ++j) {
            d[j] = c1[i].xyz[j] - c2[i].xyz[j];
            if (d[j] > 0.5)
                d[j] -= 1.;
            else if (d[j] <= -0.5)
                d[j] += 1;
        }
        dbr_xyz r;
        dbr_vec3_mult_pbc(d, pbc, r);
        dbr_real displ2 = dbr_dot3(r, r);
        diffus.add_x(displ2);
#ifdef SPLIT_DIFFUSION
        if (!cat.empty()) {
            int c = cat[i];
            if (c < SPLIT_DIFFUSION)
                cat_diffus[c].add_x(displ2);
        }
#endif
        if (old_lines.size() == (size_t) n)
            f << old_lines[i] << " ";
        f <<  displ2 << " " << sqrt(displ2) << " "
            << r[0] << " " << r[1] << " " << r[2] << endl;
    }
    cerr << "Diffusion: " << diffus.str() << endl;
#ifdef SPLIT_DIFFUSION
    if (!cat.empty()) {
        for (int i = 0; i < SPLIT_DIFFUSION; ++i)
            if (cat_diffus[i].get_n() > 100)
                cerr << "Diffusion of cat. " << i << ": "
                     << cat_diffus[i].str() << endl;
    }
#endif
}

void read_lines(string const& filename, vector<string>& old_lines)
{
    old_lines.clear();
    ifstream f(filename.c_str());
    if (!f)
        return;
    string s;
    while (getline(f, s))
        old_lines.push_back(s);
}

int main(int argc, char **argv)
{
    xyz_name *coords = 0, *coords2 = 0;
    dbr_pbc pbc;
    int old_n;
    int f_start = 1;
    bool append = false;
    vector<string> old_lines;
    while (f_start < argc && argv[f_start][0] == '-') {
        if (argv[f_start][1] == 'a')  // append to file
            append = true;
    }
    for (int i = f_start; i < argc; ++i) {
        int n = open_atoms_file(argv[i], &coords, pbc);
        if (i != 1) {
            assert(n == old_n);
            string argvi = argv[i];
            string fn_base(argvi, 0, strlen(argv[i])-4);
            string output_fn = fn_base + ".aux";
            if (append)
                read_lines(output_fn, old_lines);
            ofstream aux(output_fn.c_str());
            write_displacement_aux(aux, n, coords, coords2, pbc, old_lines);
            delete [] coords2;
        }
        old_n = n;
        coords2 = coords;
    }

    return 0;
}


