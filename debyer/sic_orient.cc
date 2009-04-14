
// Tool to find crystal (zinc-blende structure) orientation changes 
// as a function of position and/or time

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>

#include "debyer.h"
#include "fileio.h"
#include "iloops.h"
#include "utils.h"

using namespace std;

double max_disp = 5; // max. displacement
double max_ctr_dist = 10; // max. displacement

class CutSlice
{
public:
    int axis; // 0=x, 1=y, 2=z
    double from, to;

    bool check(xyz_name const& atom, dbr_pbc const& pbc) const;
};

bool CutSlice::check(xyz_name const& atom, dbr_pbc const& /*pbc*/) const
{
    //TODO use pbc
    assert (axis >= 0 && axis < 3);
    double v = atom.xyz[axis];
    if (from < to) 
        return from < v && v < to;
    else
        return v < from || v > to;
}


vector<int> find_atoms(vector<CutSlice> const& cut,
                       vector<int> const& grain_indices, int idx,
                       int n, xyz_name const* coords, dbr_pbc const& pbc)
{
    // select atoms 
    vector<int> slice_atoms;
    for (int i = 0; i < n; ++i) {
        bool ok = true;
        for (vector<CutSlice>::const_iterator j = cut.begin(); 
                                                      j != cut.end(); ++j) {
            if (!j->check(coords[i], pbc)) {
                ok = false;
                break;
            }
        }
        if (!grain_indices.empty()) 
            if (grain_indices[i] != idx)
                ok = false;
        if (ok)
            slice_atoms.push_back(i);
    }

    return slice_atoms;
}


xyz_name find_center(vector<int> const& slice_atoms, 
                     xyz_name const* coords)
{
    xyz_name r;
    strcpy(r.name, "center");
    dbr_real const* first = coords[*slice_atoms.begin()].xyz;
    dbr_real cur;
    for (int i = 0; i < 3; ++i)
        r.xyz[i] = 0;
    for (vector<int>::const_iterator a = slice_atoms.begin(); 
                                                  a != slice_atoms.end(); ++a)
        for (int i = 0; i < 3; ++i) {
            cur = coords[*a].xyz[i];
            if (cur > first[i] + 0.5)
                cur -= 1.;
            else if (cur < first[i] - 0.5)
                cur += 1.;
            r.xyz[i] += cur;
        }
    for (int i = 0; i < 3; ++i)
        r.xyz[i] /= slice_atoms.size();
    return r;
}


void print_usage_and_exit()
{
    cerr << "Usage: ..." << endl;
    dbr_abort(EXIT_FAILURE);
}

// returns angle between two vectors a and b.
double calc_angle(const dbr_real* a, const dbr_real* b)
{
    double arg = dbr_dot3(a,b) / sqrt(dbr_dot3(a,a) * dbr_dot3(b,b));
    if (arg >= 1)
        return 0;
    //cross-product
    double c[3];
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    int sign = (c[0]+c[1]+c[2] > 0) ? 1 : -1;
    return sign * acos(arg) * 180 / M_PI;
}

// input: two points a and b in reduced coordinates and PBC
// output: vector ab in real (not reduced) coordinates
void get_diff_in_pbc(const dbr_real* a, const dbr_real* b, dbr_pbc const& pbc,
                     dbr_real* abr)
{
    dbr_xyz ab;
    dbr_diff3(a, b, ab);
    dbr_vec3_mult_pbc(ab, pbc, abr);
}

// input: two points in reduced coordinates and PBC
// output: distance^2 in real (not reduced) coordinates
double get_sq_dist_in_pbc(const dbr_real* a, const dbr_real* b, 
                          dbr_pbc const& pbc)
{
    dbr_xyz dreal;
    get_diff_in_pbc(a, b, pbc, dreal);
    return dbr_dot3(dreal, dreal);
}

double calc_angle_around_axis(const dbr_real* r1, const dbr_real* r2, int axis)
{
    dbr_real r2prim[] = { r2[0], r2[1], r2[2] };
    r2prim[axis] = r1[axis];
    return calc_angle(r1, r2prim);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        print_usage_and_exit();
    }
    const char *infn = argv[1];

    vector<CutSlice> cut;
    vector<int> grain_indices;
    int idx = -1;
    int p = 2;
    while (p < argc - 2) {
        if (strlen(argv[p]) == 1) {
            CutSlice c;
            if (argv[p][0] == 'x')
                c.axis = 0; 
            else if (argv[p][0] == 'y')
                c.axis = 1; 
            else if (argv[p][0] == 'z')
                c.axis = 2; 
            else
                break;
           c.from = strtod(argv[p+1], 0);
           c.to = strtod(argv[p+2], 0);
           cut.push_back(c);
           p += 3;
        }
        else if (argv[p] == string("idx")) {
           ifstream f(argv[p+1]);
           string s;
           grain_indices.clear();
           while (getline(f, s))
               grain_indices.push_back(strtol(s.c_str(), 0, 10));
           idx = strtol(argv[p+2], 0, 10);
           p += 3;
        }
        else if (argv[p] == string("-m")) { 
            max_disp = strtod(argv[p+1], 0);
            p += 2;
        }
        else if (argv[p] == string("-r")) { 
            max_ctr_dist = strtod(argv[p+1], 0);
            p += 2;
        }
        else {
            cerr << "Unexpected argument " << p << ": " << argv[p] << endl;
            break;
        }
    }
    if (p != argc - 1) 
        print_usage_and_exit();
    const char* infn2 = argv[p];

    xyz_name *coords = 0;
    dbr_pbc pbc;
    int n = open_atoms_file(infn, &coords, pbc);
    vector<int> slice_atoms = find_atoms(cut, grain_indices, idx,
                                         n, coords, pbc);

    cerr << slice_atoms.size() << " atoms selected.\n";
    xyz_name center = find_center(slice_atoms, coords);
    dbr_xyz ctr_angstr;
    dbr_vec3_mult_pbc(center.xyz, pbc, ctr_angstr);
    cerr << "Center: (" << ctr_angstr[0] << ", " << ctr_angstr[1] 
         << ", " << ctr_angstr[2] << ")\n";

    // read atoms
    dbr_pbc pbc2;
    xyz_name *coords2 = 0;
    int n2 = open_atoms_file(infn2, &coords2, pbc2);
    assert (n2 == n);
    xyz_name center2 = find_center(slice_atoms, coords2);

    cerr << "Center move: " 
         << sqrt(get_sq_dist_in_pbc(center.xyz, center2.xyz, pbc)) << endl;

    StdDev x_rot_sd, y_rot_sd, z_rot_sd;
    vector<double> colors(n, -1.);
    for (vector<int>::const_iterator a = slice_atoms.begin(); 
                                            a != slice_atoms.end(); ++a) {
        dbr_real* at1 = coords[*a].xyz;
        dbr_real* at2 = coords2[*a].xyz;
        double ctr_sq_dist1 = get_sq_dist_in_pbc(center.xyz, at1, pbc);
        if (get_sq_dist_in_pbc(at1, at2, pbc) < max_disp * max_disp 
                && ctr_sq_dist1 < max_ctr_dist * max_ctr_dist
                && ctr_sq_dist1 > 1 * 1) {
            dbr_xyz r1; 
            get_diff_in_pbc(at1, center.xyz, pbc, r1);
            dbr_xyz r2; 
            get_diff_in_pbc(at2, center2.xyz, pbc, r2);

            double x_rot = calc_angle_around_axis(r1, r2, 0);
            x_rot_sd.add_x(x_rot);

            colors[*a] = x_rot;

            double y_rot = calc_angle_around_axis(r1, r2, 1);
            y_rot_sd.add_x(y_rot);

            double z_rot = calc_angle_around_axis(r1, r2, 2);
            z_rot_sd.add_x(z_rot);
        }
    }
    cerr << x_rot_sd.get_n() << " atoms considered..." << endl;
    cerr << "x rot [deg]: " << x_rot_sd.str() << endl;
    cerr << "y rot [deg]: " << y_rot_sd.str() << endl;
    cerr << "z rot [deg]: " << z_rot_sd.str() << endl;

    string fn_base(string(infn), 0, strlen(infn) - 4);
    ofstream clr((fn_base + ".clr").c_str());
    for (int i = 0; i < n; ++i) {
        double v = colors[i];
        if (v == -1)
            clr << "-1 0 0\n";
        else
            clr << v/180. << " " << v/180. << " " << v/180. << endl;
    }

    return 0;
}



