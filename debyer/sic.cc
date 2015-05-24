
/// dbr_sic - utility to calculate coordination numbers of atoms,
/// so-called ring distribution and other features of zinc-blende structure.
/// Named 'sic' because it was used to study SiC structure.

#include <vector>
#include <utility>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>
#include <climits>
#include <algorithm>

#include "debyer.h"
#include "fileio.h"
#include "utils.h"

using namespace std;

static
dbr_real get_sq_dist(const dbr_real *xyz1, const dbr_real *xyz2)
{
    dbr_real dx = xyz1[0] - xyz2[0];
    dbr_real dy = xyz1[1] - xyz2[1];
    dbr_real dz = xyz1[2] - xyz2[2];
    return dx*dx + dy*dy + dz*dz;
}

static
string do_format_stats(vector<int> const& v)
{
    ostringstream s;
    for (size_t i = 0; i < v.size(); ++i)
        if (v[i])
            s << setw(6) << i << "  |";
    s << endl;
    int sum = 0;
    for (int i = 0; i < 256; ++i)
        if (v[i]) {
            s << setw(8) << v[i] << "|";
            sum += v[i];
        }
    s << endl;
    for (int i = 0; i < 256; ++i)
        if (v[i])
            s << setw(6) << round(10000. * v[i] / sum) / 100. << "% |";
    s << endl;
    return s.str();
}

/// stats about coordination numbers of atoms and 2-nd neighbours,
/// should be coupled with struct dbr_atoms, which contains atom coordinates
class NeighbourStats
{
public:
    vector<vector<double> > angles;

    NeighbourStats(dbr_cells *cells_)
        : cells(cells_), find_neigh_other(NULL), ring2_count(0), coord_array(0)
        { initialize(); }

    void find_neigbours(NeighbourStats &other, dbr_real rcut,
                        vector<int> *angle_distrib, bool set_n2=true);
    void add_self_to_n1(dbr_real rcut);

    void find_rings();

    int size() const { return n1_count.size(); }

    int get_atom_number(int cell_number, int at_in_cell_number) const
        { return cell_start[cell_number] + at_in_cell_number; }

    void set_n1(int a, int n1) { n1_count[a] = n1; }
    int get_n1(int a) const { return n1_count[a]; }

    string get_n1_stats() const
    {
        return get_name() + " coordination number statistics:\n"
            + do_format_stats(get_n1_histogram());
    }

    string get_n2_stats() const;
    string get_ring_stats(bool full) const;
    string get_name() const { return cells->name; }
    dbr_cells const* get_cells() const { return cells; }
    dbr_real const* get_coord_array() const { return coord_array; }
    vector<int> get_n1_histogram() const;
    int get_wrong_bonds(double rcut);
    int get_ring3(int n) const { return ring3_count[n]; }
    void store_indices();
    int get_index(int n) const { return indices[n]; }
    vector<pair<int,int> > const& get_neigh(int n) const { return n2[n]; }

private:
    dbr_cells *cells;
    NeighbourStats const* find_neigh_other;
    vector<unsigned char> n1_count;
    vector<unsigned char> ring3_count;
    vector<vector<pair<int,int> > > n2;
    vector<int> cell_start;
    int ring2_count;
    dbr_real *coord_array;
    vector<int> indices;

    void initialize();
};

vector<int> NeighbourStats::get_n1_histogram() const
{
    vector<int> v(256, 0);
    for (vector<unsigned char>::const_iterator i = n1_count.begin();
            i != n1_count.end(); ++i)
        ++v[*i];
    return v;
}

string print_n1_stats(NeighbourStats const& a, NeighbourStats const& b)
{
    vector<int> va = a.get_n1_histogram();
    vector<int> vb = b.get_n1_histogram();
    int suma = accumulate(va.begin(), va.end(), 0);
    int sumb = accumulate(vb.begin(), vb.end(), 0);
    ostringstream s;
    assert(va.size() == vb.size());
    assert(suma != 0 && sumb != 0);
    for (size_t i = 0; i < va.size(); ++i)
        if (va[i] || vb[i]) {
            // CN#2 | Si 2382 1.35% | C 5682 3.2% | all 8064 2.1%
            s << "CN#" << i
                << " | " << a.get_name() << " " << setw(8) << va[i] << "  "
                << setprecision(4) << va[i] * 100. / suma << "%"
                << " | " << b.get_name() << " " << setw(8) << vb[i] << "  "
                << setprecision(4) << vb[i] * 100. / sumb << "%"
                << endl;
        }
    s << endl;
    return s.str();
}


void NeighbourStats::initialize()
{
    n1_count.resize(cells->atom_count);
    n2.resize(cells->atom_count);
    ring3_count.resize(cells->atom_count);
    for (size_t i = 0; i < n1_count.size(); ++i) {
        n1_count[i] = 0;
        n2[i].reserve(12);
    }
    cell_start.resize(cells->count);
    int sum = 0;
    for (int i = 0; i < cells->count; ++i) {
        if (cells->data[i].real) {
            cell_start[i] = sum;
            sum += cells->data[i].count;
        }
        else
            cell_start[i] = INT_MAX;
    }
}

struct cell_ordering: public binary_function<int, int, bool>
{
    cell_ordering(dbr_cells const *cells) : data(cells->data) {}
    bool operator() (int a, int b) const
      { return dbr_cell_original(&data[a]) < dbr_cell_original(&data[b]); }
    dbr_cell const* const data;
};


void NeighbourStats::find_neigbours(NeighbourStats &other, dbr_real rcut,
                                    vector<int> *angle_distrib,
                                    bool set_n2/*=true*/)
{
    assert (rcut > 0);
    dbr_real rcut2 = rcut * rcut;
    vector<int> neigh;
    vector<dbr_real*> neigh_a;
    vector<int> c1_nabes(27);

    int nbins = angle_distrib ? angle_distrib->size() : 0;
    if (angle_distrib)
        other.angles.resize(other.size());
    for (int i = 0; i < other.cells->count; ++i) {
        const dbr_cell *c1 = &other.cells->data[i];
        if (!c1->real)
            continue;
        for (int j = 0; j < 27; ++j)
            c1_nabes[j] = c1->neighbours[j];
        sort(c1_nabes.begin(), c1_nabes.end(), cell_ordering(cells));
        for (int k = 0; k < c1->count; k++) {
            int on = other.get_atom_number(i, k);
            dbr_real const* first_at = c1->atoms[k];
            neigh.clear();
            neigh_a.clear();

            for (int j = 0; j < 27; ++j) {
                int c2_n = c1_nabes[j];
                const dbr_cell *c2 = &cells->data[c2_n];
                for (int m = 0; m != c2->count; ++m) {
                    dbr_real d2 = get_sq_dist(first_at, c2->atoms[m]);
                    if (d2 <= rcut2) {
                        int at = get_atom_number(c2->real ?c2_n : c2->original,
                                                 m);
                        if (set_n2) {
                            for (size_t it = 0; it != neigh.size(); ++it) {
                                n2[neigh[it]].push_back(make_pair(on, at));
                                if (angle_distrib) {
                                    dbr_real angle = dbr_get_angle(neigh_a[it],
                                                       first_at, c2->atoms[m]);
                                    int bin = int((angle / M_PI) * nbins + 0.5);
                                    if (bin < nbins) //n'th bin -> n*180/nbins.
                                        ++(*angle_distrib)[bin];
                                    other.angles[on].push_back(angle);
                                }
                            }
                            neigh_a.push_back(c2->atoms[m]);
                        }
                        neigh.push_back(at);
                    }
                }
            }

            other.set_n1(on, neigh.size());
        }
    }
    find_neigh_other = &other;
}

// TODO: optimize, don't loop over all atoms in all neighbour cell
void NeighbourStats::add_self_to_n1(dbr_real rcut)
{
    assert (rcut > 0);
    dbr_real rcut2 = rcut * rcut;
    vector<int> c1_nabes(27);

    for (int i = 0; i < this->cells->count; ++i) {
        const dbr_cell *c1 = &this->cells->data[i];
        if (!c1->real)
            continue;
        for (int j = 0; j < 27; ++j)
            c1_nabes[j] = c1->neighbours[j];
        sort(c1_nabes.begin(), c1_nabes.end(), cell_ordering(cells));
        for (int k = 0; k < c1->count; k++) {
            int on = this->get_atom_number(i, k);
            dbr_real const* first_at = c1->atoms[k];

            for (int j = 0; j < 27; ++j) {
                int c2_n = c1_nabes[j];
                const dbr_cell *c2 = &cells->data[c2_n];
                for (int m = 0; m != c2->count; ++m) {
                    dbr_real d2 = get_sq_dist(first_at, c2->atoms[m]);
                    if (d2 <= rcut2 && (c1 != c2 || k != m)) {
                        n1_count[on] ++;
                    }
                }
            }
        }
    }
}


void NeighbourStats::find_rings()
{
    ring2_count = 0;
    for (size_t i = 0; i != n2.size(); ++i) { // 1st atom is i
        vector<pair<int, int> > const& t = n2[i];
        for (size_t j = 0; j != t.size(); ++j) { // 2nd atom is t[j].second
            vector<pair<int, int> > const& jn = n2[t[j].second];
            for (size_t k = 0; k != t.size(); ++k) { // 3rd is t[k].second
                if (t[j].first == t[k].first)
                    continue;
                if (j < k && t[j].second == t[k].second)
                    ring2_count++;
                for (vector<pair<int, int> >::const_iterator m = jn.begin();
                        m != jn.end(); ++m) {
                    if (m->second == t[k].second) {
                        ring3_count[i]++;
                        ring3_count[t[j].second]++;
                        ring3_count[t[k].second]++;
                    }
                }
            }
        }
    }
}

string NeighbourStats::get_n2_stats() const
{
    int sum = 0;
    for (size_t i = 0; i < n2.size(); ++i)
        sum += n2[i].size();
    return "Avg number of " + get_name() + "-...-" + get_name()
        + " bonds (x2) per atoms is: " + S(2. * sum / n2.size()) + "\n";
}

string NeighbourStats::get_ring_stats(bool full) const
{
    vector<int> v(256, 0);
    int ring3_sum = 0;
    for (vector<unsigned char>::const_iterator i = ring3_count.begin();
            i != ring3_count.end(); ++i) {
        ring3_sum += *i;
        ++v[*i];
    }
    string s = get_name() + ": " + S(ring2_count) + " 2-fold rings, "
                + S(ring3_sum / 3) + " 3-fold rings";
    if (full)
        s += "\n" + get_name() + " 3-fold ring statistics:\n"
            + do_format_stats(v);
    else
        s += "  " + get_name() + " w/ 12 rings: " + S(v[12])
            + " (" + S(round(10000. * v[12] / n2.size()) / 100.) + "%)\n";
    return s;
}

int NeighbourStats::get_wrong_bonds(double rcut)
{
    find_neigbours(*this, rcut, NULL, false);
    int n = n1_count.size(); // atom number
    // the stats contain also n bugus "self" bonds (atom "bonded" with itself)
    return (accumulate(n1_count.begin(), n1_count.end(), 0) - n) / 2;
}

void NeighbourStats::store_indices()
{
    indices = vector<int>(n1_count.size(), -1);
    for (int j = 0; j < cells->count; ++j) {
        dbr_cell const* c = &cells->data[j];
        if (!c->real)
            continue;
        for (int k = 0; k < c->count; k++) {
            int idx = c->indices[k];
            int nr_in_ns = get_atom_number(j, k);
            assert(indices[nr_in_ns] == -1);
            indices[nr_in_ns] = idx;
        }
    }
}

bool mark_neigbours(NeighbourStats const& ns, NeighbourStats const& ns2,
                    vector<int>& id,
                    int k, int grain)
{
    int t = ns.get_index(k);
    if (id[t] == -2) {
        if (ns.get_ring3(k) == 12) { // ordered
            id[t] = grain;
            vector<pair<int,int> > const& neigh = ns.get_neigh(k);
            for (vector<pair<int,int> >::const_iterator i = neigh.begin();
                                                       i != neigh.end(); ++i) {
                if (ns2.get_ring3(i->first) == 12) {
                    id[ns2.get_index(k)] = grain;
                    mark_neigbours(ns, ns2, id, i->second, grain);
                }
                else
                    id[ns2.get_index(k)] = -1;
            }
            return true;
        }
        else // amorphous
            id[t] = -1;
    }
    return false;
}

// the algorithm doesn't work, i.e. grains are not separated
void generate_grain_ids(vector<NeighbourStats> const&nstats,
                        string const& grain_id_file)
{
    assert (nstats.size() == 2);
    cerr << "Autogenerated grain IDs -> " << grain_id_file << endl;
    vector<int> id(nstats[0].size() + nstats[1].size(), -2);
    int grain = 0;
    for (int k = 0; k < nstats[0].size(); ++k) {
        if (mark_neigbours(nstats[0], nstats[1], id, k, grain)) {
            cerr << "grains #" << grain << endl;
            ++grain;
        }
    }
    cerr << "Number of grains: " << grain << endl;
    ofstream f(grain_id_file.c_str());
    for (size_t i = 0; i < id.size(); ++i)
        f << id[i] << endl;
    f.close();
}


int main(int argc, char **argv)
{
    dbr_init(&argc, &argv);
    dbr_verbosity = 1;
    bool full_ring_stats = false;
    string ring_dump;
    string angle_distrib;
    string grain_id_file;
    int option_idx = 3;
    while (argc - option_idx >= 2) {
        if (argv[option_idx] == string("-r")) {
            ring_dump = argv[option_idx + 1];
            option_idx += 2;
        }
        else if (argv[option_idx] == string("-a")) {
            angle_distrib = argv[option_idx + 1];
            option_idx += 2;
        }
        else if (argv[option_idx] == string("-g")) {
            grain_id_file = argv[option_idx + 1];
            option_idx += 2;
        }
        else
            break;
    }

    if (option_idx != argc) {
        cerr << "Usage:\n" << argv[0] << " filename max_bondlength [OPTIONS]"
            "\n options: \n    [-a angle_distrib_file]"
                        "\n    [-r ring_dump_file]"
                        "\n    [-g grain_id_file]"
                << endl;
        dbr_abort(EXIT_SUCCESS);
    }
    char *endptr;
    const dbr_real rcut = strtod(argv[2], &endptr);
    if (*endptr != 0 || rcut <= 0) {
        cerr << "Wrong value for max_bondlength\n";
        dbr_abort(EXIT_SUCCESS);
    }
    char *infn = argv[1];
    LineInput in;
    if (!in.init(infn)) {
        cerr << in.get_error();
        dbr_abort(EXIT_FAILURE);
    }

    // read atoms
    dbr_aconf aconf = read_atoms_from_file(in, false, "");
    dbr_pbc_prop pbc_prop = get_pbc_properties(aconf.pbc);
    cerr << "Numeric density: " << aconf.n / pbc_prop.volume << endl;

    // separete different atom types (xyz_name[] -> dbr_atoms[])
    dbr_atoms *xa = 0;
    int tc = dbr_get_atoms(aconf.n, aconf.atoms, &xa, 1);
    delete [] aconf.atoms;
    if (dbr_verbosity > 0)
        cerr << "Elapsed " << dbr_get_elapsed() <<" s. Atoms were sorted."
            << endl;
    if (dbr_verbosity >= 0) {
        cerr << "Atom statistics:";
        for (int i = 0; i < tc; ++i)
            cerr << " " << xa[i].name << ":" << xa[i].count;
        cerr << endl;
    }

    //put atoms into cells and free dbr_atoms[] xa
    dbr_cells *cells = prepare_cells_all(aconf.pbc, 10, xa, tc);

    vector<NeighbourStats> nstats;
    for (int i = 0; i < tc; ++i)
        nstats.push_back(NeighbourStats(&cells[i]));
    if (dbr_verbosity > 0)
        cerr << "Elapsed " << dbr_get_elapsed()
            <<" s. Memory for neighbours reserved." << endl;
    if (dbr_verbosity > 0)
        cerr << "max. bondlength is " << rcut << endl;
    assert (tc >= 2);

    const int nbins = 360;
    vector<int> ang1(nbins, 0);
    vector<int> ang2(nbins, 0);

    cerr << "Looking for second neighbours in..." << flush;
    cerr << "  " << nstats[0].get_name() << flush;
    nstats[0].find_neigbours(nstats[1], rcut, &ang1);
    cerr << "  " << nstats[1].get_name() << flush;
    nstats[1].find_neigbours(nstats[0], rcut, &ang2);
    cerr << " ...Done." << endl;

    //for (int i = 0; i < 2; ++i)
    //    cerr << nstats[i].get_n1_stats() << endl;
    cerr << print_n1_stats(nstats[0], nstats[1]) << endl;
    for (int i = 0; i < 2; ++i)
        cerr << nstats[i].get_n2_stats() << endl;

    cerr << "Looking for rings in...";
    for (int i = 0; i < 2; ++i)  {
        cerr << "  " << nstats[i].get_name() << flush;
        nstats[i].find_rings();
    }
    cerr << " ...Done." << endl;

    for (int i = 0; i < 2; ++i)
        cerr << nstats[i].get_ring_stats(full_ring_stats) << endl;

    if (!ring_dump.empty()) {
        cerr << "Number of rings, neighbours and angle deviation for each atom"
            "\n -> " << ring_dump << endl;
        vector<int> r3(aconf.n, -1);
        vector<int> nabes(aconf.n, -1);
        vector<StdDev> aa;
        if (!angle_distrib.empty())
            aa.resize(aconf.n);
        const double angle0 = acos(-1./3);
        for (vector<NeighbourStats>::iterator i = nstats.begin();
                                                    i != nstats.end(); ++i) {
            i->store_indices();
            for (int j = 0; j < i->size(); ++j) {
                int idx = i->get_index(j);
                r3[idx] = i->get_ring3(j);
                nabes[idx] = i->get_n1(j);
                if (!angle_distrib.empty()) {
                    vector<double> const& aj = i->angles[j];
                    for (vector<double>::const_iterator k = aj.begin();
                                                          k != aj.end(); ++k)
                        aa[idx].add_x(*k - angle0);
                }
            }
        }
        ofstream f(ring_dump.c_str());
        for (int i = 0; i < aconf.n; ++i) {
            f << r3[i] << " " << (r3[i] == 12 ? 1 : 0) << " " << nabes[i];
            if (!angle_distrib.empty())
                f << " " << aa[i].mean() * 180. / M_PI;
            f << endl;
        }
    }

    if (!angle_distrib.empty()) {
        cerr << "Angle distrib (" << nbins << " bins) -> "
            << angle_distrib << endl;
        ofstream f(angle_distrib.c_str());
        f << "#angle\t" << nstats[0].get_name() << "\t" << nstats[1].get_name()
            << "\tsum" << endl;
        for (int j = 0; j != nbins; ++j)
            f << j * 180./nbins << "\t" << ang1[j] << "\t" << ang2[j] << "\t"
                << ang1[j] + ang2[j] << endl;
        f.close();
    }

    if (!grain_id_file.empty()) {
        nstats[0].store_indices();
        nstats[1].store_indices();
        generate_grain_ids(nstats, grain_id_file);
    }

    cerr << "Adding bonds to atoms of the same species..." << endl;
    for (int i = 0; i < 2; ++i) {
        nstats[i].add_self_to_n1(rcut);
        cerr << nstats[i].get_n1_stats() << endl;
    }

    cerr << "Boundary only ??..." << endl;
    for (int i = 0; i < 2; ++i) {
        NeighbourStats &st = nstats[i];
        vector<int> v(256, 0);
        for (int j = 0; j < st.size(); ++j)
            if (st.get_ring3(j) != 12)
                ++v[st.get_n1(j)];
        cerr << st.get_name() + ": CN for boundary:\n"
            << do_format_stats(v) << endl;
    }

    // Wrong bonds: C-C:41 Si-Si:471
    cerr << "Wrong bonds: ";
    for (int i = 0; i < 2; ++i)
        cerr << nstats[i].get_name() << "-" << nstats[i].get_name() << flush
            << ":" << nstats[i].get_wrong_bonds(rcut) << " ";
    cerr << endl;

    free_cells_all(cells, tc);
    dbr_finalize();
    return 0;
}


