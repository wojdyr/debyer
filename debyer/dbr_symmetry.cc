// Tool to check translational symmetry of the system in given directions.

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <cfloat>
#include <map>

#include "fileio.h"
#include "utils.h"

using namespace std;

//------------------- C O N F I G U R A T I O N ---------------------------
// there can't be two atoms in cube cell of this size
const dbr_real cs_cell_size = 0.5; //[A]

// min. interactomic distance in the final structure
const dbr_real min_dist = 1.0; // [A]
//-------------------------------------------------------------------------

void print_stoichiometry_info(map<string,int> const& st)
{
    cout << "Stoichiometry:";
    int sum = 0;
    for (map<string,int>::const_iterator i = st.begin(); i != st.end(); ++i)
        sum += i->second;
    for (map<string,int>::const_iterator i = st.begin(); i != st.end(); ++i)
        cout << "   " << i->first << ":" << i->second
             << " (" << 100. * i->second / sum << "%)";
    cout << endl;
}

map<string,int> get_stoichiometry(dbr_aconf const& aconf)
{
    map<string,int> stoich;
    for (int i = 0; i < aconf.n; ++i)
        stoich[aconf.atoms[i].name] += 1;
    return stoich;
}

map<string,int> get_stoichiometry_of(dbr_aconf const& aconf,
                                     vector<int> const& sel)
{
    map<string,int> stoich;
    for (vector<int>::const_iterator i = sel.begin(); i != sel.end(); ++i)
        stoich[aconf.atoms[*i].name] += 1;
    return stoich;
}

int find_atom_near_0(dbr_aconf const& aconf)
{
    cerr << "Searching for atom near (0,0,0)...  ";
    int a0 = 0;
    dbr_real *xyz = aconf.atoms[0].xyz;
    dbr_real min_d = xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
    for (int i = 1; i < aconf.n; ++i) {
        xyz = aconf.atoms[i].xyz;
        dbr_real d = xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
        if (d < min_d)
            a0 = i;
    }
    cerr << "atom " << a0 << " (" << aconf.atoms[a0].name << ")" << endl;
    return a0;
}

bool contains_element(std::vector<dbr_real> const& vec, dbr_real t, dbr_real epsilon)
{
    for (vector<dbr_real>::const_iterator i = vec.begin(); i != vec.end(); ++i)
        if (fabs(*i - t) < epsilon)
            return true;
    return false;
}


dbr_real wrapped_diff(dbr_real x2, dbr_real x1)
{
    double d = x2 - x1;
    if (d < 0)
        d += 1.;
    assert (d >= 0.);
    assert (d < 1.);
    return d;
}

dbr_real wrapped_dist(dbr_real x1, dbr_real x2)
{
    dbr_real d = fabs(x2 - x1);
    assert (d <= 1.);
    return d < 0.5 ? d : 1. - d;
}

dbr_real wrapped_dist3(dbr_real *x1, dbr_real *x2)
{
    dbr_real d[3];
    for (int k = 0; k < 3; ++k)
        d[k] = wrapped_dist(x1[k], x2[k]);
    return sqrt(dbr_dot3(d, d));
}

dbr_real wrapped_dist3(dbr_real *x1, dbr_real *x2, dbr_pbc const& pbc)
{
    dbr_real d_reduced[3];
    for (int k = 0; k < 3; ++k)
        d_reduced[k] = wrapped_dist(x1[k], x2[k]);
    dbr_real d[3];
    dbr_vec3_mult_pbc(d_reduced, pbc, d);
    return sqrt(dbr_dot3(d, d));
}

dbr_real wrapped_sum(dbr_real x1, dbr_real d)
{
    assert(d >= 0.);
    dbr_real x2 = x1 + d;
    if (x2 >= 1)
        x2 -= 1.;
    assert (x2 >= 0.);
    assert (x2 < 1.);
    return x2;
}


int mod(int a, int n)
{
    // pre: n > 0
    int r = a % n;
    return r >= 0 ? r : r + n;
}


// 3D table with indices of atoms; each cell contains 0 or 1 atom.
// Atomic coordinates are assumed to be relative to PBC (i.e. in [0,1) range)
//  or PBC box size is assumed to be 1.0 x 1.0 x 1.0
class Cells
{
public:
    Cells(int n0, int n1, int n2) : coords(NULL) { init(n0, n1, n2); }

    Cells(dbr_real a0, dbr_real a1, dbr_real a2, dbr_real max_cell_size,
          xyz_name const* coords_, int coords_size=0)
    {
        init( (int) ceil(a0 / max_cell_size),
              (int) ceil(a1 / max_cell_size),
              (int) ceil(a2 / max_cell_size) );
        coords = coords_;
        for (int i = 0; i < coords_size; ++i) {
            int r = put(i);
            assert(r == -1);
        }
    }

    void init(int n0, int n1, int n2)
    {
        ncells[0] = n0;
        ncells[1] = n1;
        ncells[2] = n2;
        data.resize(n0*n1*n2, -1);
        set_scale(1, 1, 1);
        cout << "[Cells] " << n0 << " x " << n1 << " x " << n2 << endl;
    }

    void set_coords(xyz_name const* coords_) { coords = coords_; }

    void clear_data() { fill(data.begin(), data.end(), -1); }

    int get_atom_at(dbr_real const* x, dbr_real *pbc, dbr_real epsilon) const
    {
        int k[3];
        get_indices(x, k);
        for (int i0 = k[0] - 1; i0 <= k[0] + 1; ++i0)
            for (int i1 = k[1] - 1; i1 <= k[1] + 1; ++i1)
                for (int i2 = k[2] - 1; i2 <= k[2] + 1; ++i2) {
                    int a = get(i0, i1, i2);
                    //if (a != -1 && same_xyz(x, a, epsilon))
                    if (a != -1 && dist_lt(x, coords[a].xyz, pbc, epsilon))
                        return a;
                }
        return -1;
    }

    dbr_real dist_lt(dbr_real const* x1, dbr_real const* x2,
                     dbr_real const* pbc,
                     dbr_real epsilon) const
    {
        dbr_real d[3];
        for (int k = 0; k < 3; ++k) {
            d[k] = wrapped_dist(x1[k], x2[k]);
            if (d[k] * scaling[k] > 0.5)
                d[k] = 1./scaling[k] - d[k];
            assert(d[k] * scaling[k] >= 0);
            assert(d[k] * scaling[k] <= 0.5);
            d[k] *= pbc[k];
        }
        return dbr_dot3(d, d) < epsilon * epsilon;
    }

    bool same_xyz(dbr_real const* x, int a, dbr_real epsilon) const
    {
        for (int i = 0; i < 3; ++i) {
            dbr_real wd = wrapped_dist(coords[a].xyz[i], x[i]);
            dbr_real d = wd * scaling[i] < 0.5 ? wd : 1./scaling[i] - wd;
            if (d > epsilon / ncells[i])
                return false;
        }
        return true;
    }

    void set_scale(int s0, int s1, int s2)
    {
        scaling[0] = s0;
        scaling[1] = s1;
        scaling[2] = s2;
    }

    int put(int atom)
    {
        int p = pos_c(coords[atom].xyz);
        int r = data[p];
        data[p] = atom;
        return r;
    }

    int get(dbr_real const *x) const
    {
        return data[pos_c(x)];
    }

protected:
    int ncells[3];
    vector<int> data;
    xyz_name const *coords;
    // scaling is used only to cut rectangular from the system and make
    // it periodic
    int scaling[3];

    int pos(int k0, int k1, int k2) const
    {
        return mod(k0, ncells[0]) * ncells[1] * ncells[2]
               + mod(k1, ncells[1]) * ncells[2]
               + mod(k2, ncells[2]);
    }

    void get_indices(dbr_real const *x, int *k) const
    {
        for (int i = 0; i < 3; ++i)
            k[i] = mod(int(x[i] * scaling[i] * ncells[i]), ncells[i]);
    }

    int pos_c(dbr_real const *x) const
    {
        int k[3];
        get_indices(x, k);
        return pos(k[0], k[1], k[2]);
    }

    int get(int k0, int k1, int k2) const
    {
        return data[pos(k0, k1, k2)];
    }
};


void normalize_coordinates(int n, xyz_name *coords)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < 3; ++j) {
            if (coords[i].xyz[j] >= 1. || coords[i].xyz[j] < 0.)
                coords[i].xyz[j] -= floor(coords[i].xyz[j]);
        }
}

void assert_coordinates_normalized(dbr_aconf const& aconf)
{
    for (int i = 0; i < aconf.n; ++i)
        for (int j = 0; j < 3; ++j) {
            dbr_real c = aconf.atoms[i].xyz[j];
            if (c >= 1. || c < 0.) {
                cerr << "error: coordinates not normalized (as expected)\n";
                cerr << "error: atom " << i << " xyz["<<j<<"]=" << c << endl;
                dbr_abort(EXIT_FAILURE);
            }
        }
}

// one line per atom
vector<dbr_real> read_energy_file(const char* energy_file)
{
    vector<dbr_real> result;
    ifstream efile(energy_file);
    if (!efile) {
        cerr << "energy file not found\n";
        dbr_abort(EXIT_FAILURE);
    }
    string line;
    while (getline(efile, line)) {
        dbr_real d = strtod(line.c_str(), NULL);
        result.push_back(d);
    }
    return result;
}

dbr_real get_pbc_length(dbr_pbc const& pbc, int dir)
{
    if (dir == 0)
        return sqrt(pbc.v00 * pbc.v00 + pbc.v01 * pbc.v01 + pbc.v02 * pbc.v02);
    else if (dir == 1)
        return sqrt(pbc.v10 * pbc.v10 + pbc.v11 * pbc.v11 + pbc.v12 * pbc.v12);
    else if (dir == 2)
        return sqrt(pbc.v20 * pbc.v20 + pbc.v21 * pbc.v21 + pbc.v22 * pbc.v22);
    assert(0);
    return 0;
}

void check_symmetry(const char *filename)
{
    // used to compare coordinates of real and perfect atom images
    const dbr_real epsilon = 1e-2;
    dbr_aconf aconf = read_atoms_from_file(filename, true);
    //normalize_coordinates(aconf.n, aconf.atoms);
    assert_coordinates_normalized(aconf);
    xyz_name *coords = aconf.atoms;

    dbr_real pbc[] = { get_pbc_length(aconf.pbc, 0),
                       get_pbc_length(aconf.pbc, 1),
                       get_pbc_length(aconf.pbc, 2) };
    Cells cells(pbc[0], pbc[1], pbc[2], cs_cell_size, aconf.atoms, aconf.n);

    // pick one atom
    int a0 = find_atom_near_0(aconf);
    cerr << "Checking translation symmetry..."  << endl;
    dbr_real *xyz0 = coords[a0].xyz;

    // find number of images and translation vectors
    vector<dbr_real> imgs[3];
    for (int i = 0; i < aconf.n; ++i) {
        if (i == a0)
            continue;
        dbr_real diff[3];
        for (int j = 0; j < 3; ++j)
            diff[j] = wrapped_diff(coords[i].xyz[j], xyz0[j]);
        bool same[3];
        for (int j = 0; j < 3; ++j)
            same[j] = fabs(diff[j]) < epsilon / pbc[j];
        if ((int)same[0] + (int)same[1] + (int)same[2] == 2) {
            int not_same;
            if (!same[0])
                not_same = 0;
            else if (!same[1])
                not_same = 1;
            else /*!same[2]*/
                not_same = 2;
            imgs[not_same].push_back(diff[not_same]);
        }
    }
    int mults[3];
    for (int i = 0; i < 3; ++i) {
        vector<dbr_real>& img = imgs[i];
        sort(img.begin(), img.end());
        int m = img.size() + 1;
        int mult = 1;
        for (int j = m; j > 1 && mult == 1; --j) {
            if (m % j != 0)
                continue;
            mult = j;
            // check images of atom a0
            for (int k = 1; k < j; ++k)
                if (!contains_element(img, dbr_real(k)/j, epsilon / pbc[i])) {
                    mult = 1;
                    break;
                }
            if (mult == 1)
                continue;
            // check all atoms
            for (int a = 0; a < aconf.n; ++a) {
                dbr_real x[3];
                for (int ii = 0; ii < 3; ++ii)
                    x[ii] = coords[a].xyz[ii];
                x[i] = wrapped_sum(x[i], 1. / j);
                int b = cells.get_atom_at(x, pbc, epsilon);
                if (b == -1 || strcmp(coords[a].name, coords[b].name) != 0) {
                    mult = 1;
                    break;
                }
            }
        }
        mults[i] = mult;
    }
    cerr << "number of cells: "
         << mults[0] << " x " << mults[1] << " x " << mults[2] << endl;
    cerr << "primitive cell size: " << aconf.pbc.v00 / mults[0] << " x "
         << aconf.pbc.v11 / mults[1] << " x " << aconf.pbc.v22 / mults[2]
         << endl;
}

struct atom_idx_sorter {
    atom_idx_sorter(xyz_name *coords_, int t_) : coords(coords_), t(t_) {}
    bool operator() (int i, int j)
        { return coords[i].xyz[t] < coords[j].xyz[t]; }

    xyz_name *coords;
    int t;
};

class CyclicVectorInt : public vector<int>
{
public:
    CyclicVectorInt() {}
    CyclicVectorInt(size_type n) : vector<int>(n) {}
    template <class InputIterator>
    CyclicVectorInt(InputIterator a, InputIterator b) : vector<int>(a,b) {}

    int operator[] (int n) const
    {
        int r = n % (int) size();
        return vector<int>::operator[] (r >= 0 ? r : r + size());
    }
    reference at(size_type n)
    {
        assert(n < size());
        return vector<int>::operator[] (n);
    }
};

dbr_real find_largest_gap(vector<dbr_real> const& zz,
                          dbr_real dmin, dbr_real dmax)
{
    vector<dbr_real>::const_iterator a = lower_bound(zz.begin(),zz.end(), dmin);
    vector<dbr_real>::const_iterator b = lower_bound(zz.begin(),zz.end(), dmax);
    assert (a != zz.begin() && a != zz.end());
    assert (b != zz.begin() && b != zz.end());
    dbr_real max_delta = 0.;
    dbr_real r = -1.;
    for (vector<dbr_real>::const_iterator i = a; i <= b; ++i) {
        dbr_real delta = *i - *(i-1);
        if (delta > max_delta) {
            max_delta = delta;
            r = (*i + *(i-1)) / 2.;
        }
    }
    //cerr << "max_delta " << max_delta << "  " << r << endl;
    return r;
}

CyclicVectorInt select_interesting_atoms(dbr_aconf const& aconf)
{
    // filter atoms
    CyclicVectorInt sel; // positions of atoms taken into account
    const dbr_real dmin = 0.15;
    const dbr_real dmax = 0.20;
    vector<dbr_real> zz(aconf.n);
    for (int i = 0; i < aconf.n; ++i)
        zz[i] = aconf.atoms[i].xyz[2];
    sort(zz.begin(), zz.end());
    dbr_real dlower = find_largest_gap(zz, dmin, dmax);
    dbr_real dupper = find_largest_gap(zz, 1 - dmax, 1 - dmin);

    for (int i = 0; i < aconf.n; ++i) {
        dbr_real z = aconf.atoms[i].xyz[2];
        if (z < dlower || z > dupper)
            sel.push_back(i);
    }
    return sel;
}

class CyclicVectorIntSlice
{
public:
    CyclicVectorIntSlice(CyclicVectorInt const& ref_, int start_, int length_)
        : ref(ref_), start(start_), length(length_) {}
    int size() const { return length; }
    int operator[] (size_t n) const { return ref[start+n]; }
private:
    CyclicVectorInt const& ref;
    int start, length;
};

// points A B C D; returns true if |BC| < r < |AD| in 1D; the points are
// coordinates in direction `dir' of given atoms
bool is_distance_proper(dbr_aconf const& aconf, int dir, dbr_real r,
                        int a, int b, int c, int d)
{
    dbr_real xa = aconf.atoms[a].xyz[dir];
    dbr_real xb = aconf.atoms[b].xyz[dir];
    dbr_real xc = aconf.atoms[c].xyz[dir];
    dbr_real xd = aconf.atoms[d].xyz[dir];
    dbr_real minr = wrapped_diff(xc, xb);
    dbr_real maxr = wrapped_diff(xd, xa);
    return minr < r && r < maxr;
}

/*
bool can_be_safely_pbced(dbr_aconf const& aconf, int dir,
                         CyclicVectorIntSlice const& primitive, dbr_real r,
                         dbr_real x0)
{
    static vector<int> left;
    static vector<int> right;
    dbr_real min_dist_reduced = min_dist / get_pbc_length(aconf.pbc, dir);
    //dbr_real x0 = aconf.atoms[primitive[0]].xyz[dir];
    left.clear();
    right.clear();
    for (int i = 0; i != primitive.size(); ++i) {
        int idx = primitive[i];
        dbr_real d0 = wrapped_diff(aconf.atoms[idx].xyz[dir], x0);
        assert (d0 >= 0 && d0 < r);
        if (d0 < min_dist_reduced)
            left.push_back(idx);
        else if (d0 > r - min_dist_reduced)
            right.push_back(idx);
    }
    for (vector<int>::const_iterator i = left.begin(); i != left.end(); ++i) {
        dbr_real const* xyz = aconf.atoms[*i].xyz;
        dbr_real img[3] = { xyz[0], xyz[1], xyz[2] };
        img[dir] = wrapped_sum(xyz[dir], r);
        for (vector<int>::const_iterator j = right.begin();
                                                       j != right.end(); ++j) {
            dbr_real const* xyz2 = aconf.atoms[*j].xyz;
            dbr_real d_reduced[3];
            for (int k = 0; k < 3; ++k)
                d_reduced[k] = wrapped_dist(img[k], xyz2[k]);
            dbr_real d[3];
            dbr_vec3_mult_pbc(d_reduced, aconf.pbc, d);
            if (dbr_dot3(d, d) < min_dist*min_dist)
                return false;
        }
    }
    return true;
}
*/

// for debugging
bool check_for_dups(dbr_aconf const& aconf, dbr_real epsilon)
{
    int dup_counter = 0;
    for (int i = 0; i < aconf.n - 1; ++i) {
        dbr_real *a = aconf.atoms[i].xyz;
        for (int j = i + 1; j < aconf.n; ++j) {
            dbr_real *b = aconf.atoms[j].xyz;
            dbr_real d = wrapped_dist3(a, b, aconf.pbc);
            if (d < epsilon) {
                cerr << "Dupl/ " << d << " " << i << " " << j
                    << " (" << a[0] << ", " << a[1] << ", " << a[2] << ")\n";
                dup_counter++;
            }
        }
    }
    if (dup_counter)
        cerr << dup_counter << " duplicates.\n";
    if (dup_counter == 0)
        cerr << "DEBUG: No duplicates.\n";
    return dup_counter != 0;
}

/*
bool check_distances(dbr_aconf const&aconf, int dir, Cells const* cells,
                     CyclicVectorInt const& sel, int nparts, int i)
{
    // if the difference in coordinates in analyzed direction of any two atoms
    // is less than min_delta, these atoms shall not be split.
    //const dbr_real min_delta = 1e-6; // reduced unit (0,1)

    // min. interactomic distance in the final structure
    const dbr_real min_dist = 1.0; // [A]

    dbr_real dx = 1.0 / nparts;
    int ac = sel.size() / nparts;

    dbr_real x1m = aconf.atoms[sel[i]     ].xyz[dir];
    dbr_real x1M = aconf.atoms[sel[i+1]   ].xyz[dir];
    dbr_real x2m = aconf.atoms[sel[ac+i]  ].xyz[dir];
    dbr_real x2M = aconf.atoms[sel[ac+i+1]].xyz[dir];
    dbr_real delta_x1 = wrapped_diff(x1M, x1m);
    dbr_real delta_x2 = wrapped_diff(x2M, x2m);

    //if (delta_x1 < min_delta || delta_x2 < min_delta)
    //    return false;

    dbr_real dxm = wrapped_diff(x2m, x1M);
    dbr_real dxM = wrapped_diff(x2M, x1m);
    // cout << "n" << name0 << ": " << an0 << "  E= " << slice_e
    //      << "  range: " << x1M << " - " << x2m
    //      << "  delta: " << (dxm < dx && dx < dxM ? "ok" : "wrong") << endl;
    if (!(dxm < dx && dx < dxM))
        return false;

    // check minimal interatomic distance (atoms in cell and its image)
    dbr_real min_dist_reduced = min_dist / get_pbc_length(aconf.pbc, dir);
    for (int t = i+1; t != ac+i+1; ++t) {
        dbr_real const* xyz = aconf.atoms[sel[t]].xyz;
        if (wrapped_diff(xyz[dir], x1M) > min_dist_reduced)
            break;
        dbr_real img[3] = { xyz[0], xyz[1], xyz[2] };
        img[dir] = wrapped_sum(xyz[dir], dx);
        int r = cells->get_atom_at(img, min_dist_reduced);
        if (r != -1)
            return false;
    }
    return true;
}
*/

void copy_atom(xyz_name const& source, xyz_name& dest)
{
    strcpy(dest.name, source.name);
    for (int k = 0; k < 3; ++k)
        dest.xyz[k] = source.xyz[k];
}

void copy_atom(xyz_name *coords, int source_pos, int dest_pos)
{
    copy_atom(coords[source_pos], coords[dest_pos]);
}

// Atoms from the slice with the lowest energy are duplicated to other slices.
// Only atoms that were selected (i.e. that are in `sel') are changed.
void multiply_selected_atoms(xyz_name *coords, int dir,
                             CyclicVectorInt const& sel,
                             int min_pos, int nparts)
{
    assert (sel.size() % nparts == 0);
    int ac = sel.size() / nparts;
    for (int j = 0; j < ac; ++j) {
        int orig_pos = sel[min_pos+j];
        for (int i = 1; i < nparts; ++i) {
            int copy_pos = sel[i * ac + min_pos+j];
            copy_atom(coords, orig_pos, copy_pos);
            coords[copy_pos].xyz[dir]
                  = wrapped_sum(coords[orig_pos].xyz[dir], double(i) / nparts);
        }
    }
}

int axis2dir(char axis)
{
    int dir = -1;
    if (axis == 'x')
        dir = 0;
    else if (axis == 'y')
        dir = 1;
    //else if (axis == 'z')
    //    dir = 2;
    else
        assert(0);
    return dir;
}
char dir2axis(int dir)
{
    if (dir == 0)
        return 'x';
    else if (dir == 1)
        return 'y';
    else if (dir == 2)
        return 'z';
    assert(0);
    return 0;
}

int count_atoms(vector<int>::const_iterator begin,
                vector<int>::const_iterator end,
                xyz_name const* coords, const char *name0)
{
    int n = 0;
    for (vector<int>::const_iterator i = begin; i != end; ++i)
        if (strcmp(coords[*i].name, name0) == 0)
            ++n;
    return n;
}

void get_primitive(xyz_name const* coords, CyclicVectorInt const& sel, int pos,
                   int dir, dbr_real dmax1, int dir2, dbr_real dmax2,
                   vector<int> &primitive)
{
    primitive.clear();
    dbr_real const* a = coords[sel[pos]].xyz;
    for (int j = pos; ; ++j) {
        dbr_real const* b = coords[sel[j]].xyz;
        if (wrapped_diff(b[dir], a[dir]) >= dmax1)
            break;
        //dbr_real d1 = wrapped_diff(b[dir], a[dir]);
        if (dir2 >= 0 && wrapped_diff(b[dir2], a[dir2]) >= dmax2)
            continue;
        primitive.push_back(sel[j]);
    }
}

void find_lowest_energy(char axis, int nparts,
                        char axis2, int nparts2,
                        const char* input_file, const char *energy_aux_file,
                        const char* output_file)
{
    int dir = axis2dir(axis);
    int dir2 = axis2 ? axis2dir(axis2) : -1;

    dbr_aconf aconf = read_atoms_from_file(input_file, true);
    assert_coordinates_normalized(aconf);
    xyz_name *coords = aconf.atoms;
    //normalize_coordinates(n, coords);

    // energies from aux file
    vector<dbr_real> energies = read_energy_file(energy_aux_file);
    assert((size_t) aconf.n == energies.size());

    dbr_real min_e = FLT_MAX;
    int min_pos = -1;

    CyclicVectorInt sel = select_interesting_atoms(aconf);
    assert (sel.size() % (nparts * nparts2) == 0);

    // calculate average energy in selected atoms
    double avg_e = 0.;
    for (vector<int>::const_iterator i = sel.begin(); i != sel.end(); ++i)
        avg_e += energies[*i];
    avg_e /= sel.size();

    // Count atoms of one of the species. This is enough to check
    // for stoichiometry in binary system.
    const char *name0 = coords[0].name;
    int an0 = count_atoms(sel.begin(), sel.end(), coords, name0);
    assert (an0 % (nparts * nparts2) == 0);

    cerr << sel.size() << " of " << aconf.n << " atoms considered, including "
         << an0 << " " << name0 << " atoms ("
         << 100. * an0 / sel.size() << "%)." << endl;

    sort(sel.begin(), sel.end(), atom_idx_sorter(coords, dir));
    int ac = sel.size() / nparts2 / nparts;
    int pc_tried_cnt = 0; // counter of all primitive cells we try

    // check if the primitive cell can be PBCed -- part 1 (prepare)
    int scale[3] = { 1, 1, 1 };
    scale[dir] = nparts;
    scale[dir2] = nparts2;
    dbr_real pbc[3];
    dbr_real prim_box[3];
    for (int i = 0; i < 3; ++i) {
        pbc[i] = get_pbc_length(aconf.pbc, i);
        prim_box[i] = pbc[i] / scale[i];
    }
    Cells prim_cells(prim_box[0], prim_box[1], prim_box[2], min_dist,
                     aconf.atoms);
    prim_cells.set_scale(scale[0], scale[1], scale[2]);

    vector<int> pr(ac);
    for (size_t i = 0; i != sel.size(); ++i) {
        get_primitive(coords, sel, i,
                      dir, 1./nparts, dir2, 1./nparts2,
                      pr);

        if ((int) pr.size() != ac)
            continue;

        dbr_real slice_e = 0.;
        for (vector<int>::const_iterator j = pr.begin(); j != pr.end(); ++j)
            slice_e += energies[*j];
        // we are interested only in minimal energy
        if (slice_e >= min_e)
            continue;

        // stoichiometry in slice should be the same as total stoichiometry
        int slice_n0 = count_atoms(pr.begin(), pr.end(), coords, name0);
        if (slice_n0 * nparts * nparts2 != an0)
            continue;

        pc_tried_cnt++;

        // check if the primitive cell can be PBCed -- part 2 (check)
        prim_cells.clear_data();
        bool ok = true;
        for (size_t j = 0; j != pr.size(); ++j) {
            int a = pr[j];
            if (prim_cells.get_atom_at(aconf.atoms[a].xyz, pbc, min_dist)
                    != -1) {
                ok = false;
                break;
            }
            else {
                int r = prim_cells.put(a);
                if (r != -1) {
                    ok = false;
                    break;
                }
            }
        }
        if (!ok)
            continue;

        min_e = slice_e;
        min_pos = i;
    }

    cerr << pc_tried_cnt << " primitive cells were investigated." << endl;

    if (min_e == FLT_MAX) {
        cerr << "primitive cell NOT FOUND" << endl;
        return;
    }

    cerr << "E min/avg: " << min_e/ac << "  " << avg_e << endl;

    get_primitive(coords, sel, min_pos,
                  dir, 1./nparts, dir2, 1./nparts2,
                  pr);
    cerr << "min slice starts at "
         << dir2axis(dir) << "=" << coords[sel[min_pos]].xyz[dir];
    if (dir2 >= 0)
        cerr << ", " << dir2axis(dir2) << "=" << coords[sel[min_pos]].xyz[dir2];
    cerr << " (#" << sel[min_pos] << ")" << endl;

    if (output_file) {
        // store atoms from primitive cell with the lowest E in new array
        vector<xyz_name> aprim(ac);
        assert((int) pr.size() == ac);
        for (int i = 0; i < ac; ++i)
            aprim[i] = aconf.atoms[pr[i]];

        // build coords using primitive cell aprim
        int counter = 0;
        for (int i = 0; i < nparts; ++i)
            for (int j = 0; j < nparts2; ++j)
                for (int k = 0; k < ac; ++k) {
                    assert (counter == (i * nparts2 + j) * ac + k);
                    xyz_name &copy = coords[sel[(i * nparts2 + j) * ac + k]];
                    copy_atom(aprim[k], copy);
                    copy.xyz[dir]
                        = wrapped_sum(copy.xyz[dir], double(i) / nparts);
                    copy.xyz[dir2]
                        = wrapped_sum(copy.xyz[dir2], double(j) / nparts2);
                    counter++;
                }
        assert (counter == (int) sel.size());

        //check_for_dups(aconf, min_dist);
        print_stoichiometry_info(get_stoichiometry(aconf));
        //print_stoichiometry_info(get_stoichiometry_of(aconf, sel2));
        write_file_with_atoms(aconf, output_file);
    }
}

int main(int argc, char **argv)
{
    if (argc == 3 && strcmp(argv[1], "check") == 0)
        check_symmetry(argv[2]);
    else if (argc >= 5 && strlen(argv[1]) == 1) {
        char axis = argv[1][0];
        int nparts = strtol(argv[2], NULL, 10);
        char axis2 = 0;
        int nparts2 = 1;
        if (argc >= 7) {
            axis2 = argv[3][0];
            nparts2 = strtol(argv[4], NULL, 10);
        }
        int file0 = (argc >= 7 ? 5 : 3);
        const char *output_file = (argc > file0 + 2 ? argv[file0 + 2] : NULL);
        find_lowest_energy(axis, nparts, axis2, nparts2,
                           argv[file0], argv[file0+1], output_file);
    }
    else {
        cerr << "Usage: \n"
                "    dbr_symmetry check input.cfg\n"
                "    dbr_symmetry D N [D N] input.cfg energy.aux [output.cfg]\n"
                "D is x or y or z \n"
                "N is integer\n";
        return -1;
    }
    return 0;
}

