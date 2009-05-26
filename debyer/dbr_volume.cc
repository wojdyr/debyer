// Tool to calculate volume and a number of voids

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <cstring>

#include "fileio.h"
#include "utils.h"

using namespace std;

static const char* usage =
"Usage:\n"
"       dbr_volume [-m N] FILENAME A [AMAX ASTEP]\n"
"   N = min. number of empty cells that make void\n"
"   A = size of cubic cell\n"
"   if AMAX and ASTEP are given, a = A, A+ASTEP, ... AMAX.\n";

int main(int argc, char **argv)
{
    int c = 3;
    int first_param = 1;
    while (first_param + 2 < argc) {
        if (strcmp(argv[first_param], "-m") == 0) {
            c = strtol(argv[2], 0, 10);
            first_param += 2;
        }
        else
            break;
    }

    if (first_param + 2 != argc && first_param + 4 != argc) {
        cerr << usage;
        return -1;
    }
    double a0 = strtod(argv[first_param + 1], 0);
    double astep = 0,
           amax = 0;
    if (first_param + 4 != argc) {
        amax = strtod(argv[first_param + 2], 0);
        astep = strtod(argv[first_param + 3], 0);
    }

    dbr_atom *coords = 0;
    dbr_pbc pbc;
    int n = open_atoms_file(argv[first_param], &coords, pbc);
    int has_pbc = (pbc.v00 != 0.);
    assert (has_pbc);

    double a = a0;
    int nx = iround(pbc.v00 / a);
    int ny = iround(pbc.v11 / a);
    int nz = iround(pbc.v22 / a);
    vector<vector<vector<int> > >
        boxes(nx, vector<vector<int> >(ny, vector<int>(nz, 0)));
    for (int i = 0; i < n; ++i) {
       int jx = int(coords[i].xyz[0] * nx) % nx;
       int jy = int(coords[i].xyz[1] * ny) % ny;
       int jz = int(coords[i].xyz[2] * nz) % nz;
       boxes[jx][jy][jz] += 1;
    }
    int all = nx * ny * nz;
    int empty = 0;
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                if (boxes[nx][ny][nz] == 0)
                    empty++;
    int voids = 0;
    cout << "# dx  dy  dz  nboxes  empty  V[%]  voids\n";
    cout << "# " << pbc.v00 / nx << " " << pbc.v11 / ny << " "
         << pbc.v22 / nz  << " " << all << " " << empty << " "
         << 100. * empty / all << " " << voids << endl;

    delete [] coords;
    return 0;
}


