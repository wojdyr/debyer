//  debyer -- program for calculation of diffration patterns
//  Copyright (C) 2006-2007 Marcin Wojdyr
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

#include "fileio.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cctype>
#include <cassert>
#include <map>
#include <vector>
#include <algorithm>

#include "atomtables.h"

using namespace std;


#define mcerr \
    if (dbr_nid == 0) \
        cerr

string argv_as_str(int argc, char **argv)
{
    string s;
    for (int i = 0; i < argc; ++i) {
        if (i != 0)
            s += " ";
        bool need_quote = false;
        for (const char* c = argv[i]; *c != '\0'; ++c)
            if (strchr("|&;<>()$`\\' \t\n*?[#~%", *c) != NULL) {
                need_quote = true;
                break;
            }
        if (need_quote)
            s += string("\"") + argv[i] + "\"";
        else
            s += argv[i];
    }
    return s;
}

bool is_xyz_format(char const* buffer)
{
    // find new lines
    const char *e[5];
    for (int i = 0; i < 5; ++i) {
        char const* start = i > 0 ? e[i-1]+1 : buffer;
        e[i] = strchr(start, '\n');
        if (!e[i])
            return false;
    }
    char* endptr = 0;
    strtol(buffer, &endptr, 10);
    while (isspace(*endptr))
        ++endptr;
    if (endptr < e[0])
        return false;
    for (int i = 2; i < 5; ++i) {
        const char *ptr = e[i-1];
        while (isspace(*ptr))
            ++ptr;
        if (!isalpha(*ptr)) //should it start with letter?
            return false;
        while (isalnum(*ptr)) //atom name
            ++ptr;
        strtod(ptr, &endptr); //x
        ptr = endptr;
        strtod(ptr, &endptr); //y
        ptr = endptr;
        strtod(ptr, &endptr); //z
        ptr = endptr;
        if (ptr > e[i])
            return false;
        while (isspace(*ptr))
            ++ptr;
        if (ptr < e[i])
            return false;
    }
    return true;
}

bool is_plain_format(char const* /*buffer*/)
{
    return true;  //TODO
}

void read_xyz(LineInput& in, dbr_aconf *aconf)
{
    const char* line = in.get_line(); //first line
    aconf->n = strtol(line, 0, 10);

    line = in.get_line(); //2nd line
    if (dbr_verbosity >= 0)
        mcerr << aconf->n << " atoms, title: " << line << endl;

    aconf->atoms = new dbr_atom[aconf->n];
    for (int i = 0; i < aconf->n; ++i) {
        dbr_atom& atom = aconf->atoms[i];
        line = in.get_line();
        if (!line) {
            mcerr << "Error: Reading line " << i+3 << " failed." << endl;
            mcerr << in.get_error() << endl;
            dbr_abort(EXIT_FAILURE);
        }
        if (sscanf(line, "%7s "DBR_F" "DBR_F" "DBR_F, atom.name, &atom.xyz[0],
                    &atom.xyz[1], &atom.xyz[2]) != 4) {
            mcerr << "Format error in line " << i+3 << ":" << endl
                  << line << endl;
            dbr_abort(EXIT_FAILURE);
        }
    }
}


#define ASSERT_FORMAT(condition) \
    if (!(condition)) { \
        mcerr << "Error in line: " << in.get_line_number() << endl; \
        dbr_abort(EXIT_FAILURE); \
    }


static
void resize_atoms(size_t& atoms_size, dbr_atom** atoms)
{
    dbr_atom *new_atoms = new dbr_atom[atoms_size*2];
    for (size_t i = 0; i < atoms_size; ++i)
        new_atoms[i] = (*atoms)[i];
    atoms_size *= 2;
    delete [] *atoms;
    *atoms = new_atoms;
}

void read_plain(LineInput& in, dbr_aconf *aconf)
{
    size_t atoms_size = 1024;

    size_t counter = 0;
    aconf->atoms = new dbr_atom[atoms_size];
    const char *line;
    while ((line = in.get_line())) {
        if (counter == atoms_size)
            resize_atoms(atoms_size, &aconf->atoms);
        dbr_atom& atom = aconf->atoms[counter];
        ++counter;
        const char *ptr = line;
        while (isspace(*ptr))
            ++ptr;
        string name;
        int r = 0;
        if (isalpha(*ptr)) { //name, x, y, z
            r = sscanf(line, "%7s "DBR_F" "DBR_F" "DBR_F,
                        atom.name, &atom.xyz[0], &atom.xyz[1], &atom.xyz[2]);
        }
        else if (*ptr) {//x, y, z, name
            r = sscanf(line, DBR_F" "DBR_F" "DBR_F" %7s",
                        &atom.xyz[0], &atom.xyz[1], &atom.xyz[2], atom.name);
        }
        if (!r) {
            mcerr << "Error in line " << counter << endl;
            dbr_abort(EXIT_FAILURE);
        }
    }
    aconf->n = counter;
}

void write_atoms_to_xyz_file(dbr_aconf const& aconf, string const& filename)
{
    if (dbr_nid != 0)
        return;
    if (dbr_verbosity >= 0)
        mcerr << "Writing atoms to XMOL xyz file: " << filename << endl;
    ofstream f(filename.c_str(), ios::out|ios::binary);
    if (!f) {
        mcerr << "Error. Can't open file " << filename << " for writing.\n";
        return;
    }
    f << aconf.n << "\nconverted by debyer\n";
    for (int i = 0; i < aconf.n; ++i) {
        dbr_atom const& atom = aconf.atoms[i];
        f << atom.name << " "
          << atom.xyz[0] << " " << atom.xyz[1] << " " << atom.xyz[2] << endl;
    }
    f.close();
}

static bool endswith(string const& s, char const* e)
{
    int slen = s.size();
    int elen = strlen(e);
    return slen >= elen && strcmp(s.c_str() + slen - elen, e) == 0;
}

string detect_input_format(LineInput &in)
{
    char const* cfg_magic = "Number of particles =";
    if (!strncmp(in.get_buffer(), cfg_magic, strlen(cfg_magic))) {
        if (dbr_verbosity > 0)
            mcerr << "Detected AtomEye CFG format" << endl;
        return "atomeye";
    }
    else if (endswith(in.get_orig_filename(), "CONFIG")
             || endswith(in.get_orig_filename(), "REVCON")) {
        if (dbr_verbosity > 0)
            mcerr << "DL_POLY CONFIG format" << endl;
        return "dlpoly";
    }
    else if (endswith(in.get_orig_filename(), ".lammps")
             || endswith(in.get_orig_filename(), ".lmps")) {
        if (dbr_verbosity > 0)
            mcerr << "LAMMPS data input format" << endl;
        return "lammps";
    }
    else if (is_xyz_format(in.get_buffer())) {
        if (dbr_verbosity > 0)
            mcerr << "Detected XMol XYZ format" << endl;
        return "xyz";
    }
    else if (is_plain_format(in.get_buffer())) {
        if (dbr_verbosity > 0)
            mcerr << "Detected plain format with atom coordinates" << endl;
        return "xyza";
    }
    else {
        mcerr << "Error. Unknown format of input file" << endl;
        dbr_abort(EXIT_FAILURE);
        return ""; // suppress warning
    }
}

dbr_aconf read_atoms_from_file(LineInput &in, bool reduced_coords,
                               string const& format)
{
    dbr_aconf aconf;
    // initialize dbr_aconf
    aconf.n = 0;
    aconf.reduced_coordinates = false;
    aconf.atoms = NULL;
    aconf.pbc.v00 = aconf.pbc.v01 = aconf.pbc.v02 = 0.;
    aconf.pbc.v10 = aconf.pbc.v11 = aconf.pbc.v12 = 0.;
    aconf.pbc.v20 = aconf.pbc.v21 = aconf.pbc.v22 = 0.;

    string fmt = format.empty() ? detect_input_format(in) : format;

    if (fmt == "atomeye")
        read_atomeye(in, &aconf, reduced_coords);
    else if (fmt == "dlpoly")
        read_dlpoly_config(in, &aconf);
    else if (fmt == "lammps")
        read_lammps_data(in, &aconf, reduced_coords);
    else if (fmt == "xyz")
        read_xyz(in, &aconf);
    else if (fmt == "xyza")
        read_plain(in, &aconf);

    if (reduced_coords && !aconf.reduced_coordinates && aconf.pbc.v00 != 0) {
        double H_1[3][3]; // inverse of pbc matrix
        dbr_inverse_3x3_matrix(aconf.pbc, H_1);
        dbr_xyz reduced;
        for (int i = 0; i < aconf.n; ++i) {
            dbr_vec3_mult_mat3x3(aconf.atoms[i].xyz, H_1, reduced);
            for (int j = 0; j < 3; ++j)
                // wrap to <0,1)
                aconf.atoms[i].xyz[j] = (reduced[j] - floor(reduced[j]));
        }
        aconf.reduced_coordinates = true;
    }

    if (dbr_verbosity > 0)
        mcerr << "Elapsed " << dbr_get_elapsed() << " s. Atoms were read."
            << endl;
    aconf.orig_filename = in.get_filename();
    return aconf;
}

string get_format_from_filename(string const& filename)
{
    if (endswith(filename, ".cfg"))
        return "atomeye";
    else if (endswith(filename, ".lammps") || endswith(filename, ".lmps"))
        return "lammps";
    else if (endswith(filename, ".xyz"))
        return "xyz";
    else if (endswith(filename, ".xyza"))
        return "xyza";
    else if (endswith(filename, ".pdb"))
        return "pdb";
    else if (endswith(filename, "CONFIG") || endswith(filename, "REVCON"))
        return "dlpoly";
    else {
        mcerr << "Can't guess filetype of atoms output file from filename"
            << endl;
        dbr_abort(EXIT_FAILURE);
        return ""; // suppress warning
    }
}

void write_file_with_atoms(dbr_aconf const& aconf, string const& filename,
                           string const& format)
{
    string fmt = format.empty() ? get_format_from_filename(filename) : format;

    if (fmt == "xyz")
        write_atoms_to_xyz_file(aconf, filename);
    else if (fmt == "atomeye")
        write_atoms_to_atomeye_file(aconf, filename);
    else if (fmt == "dlpoly")
        write_dlpoly_file(aconf, filename);
    else if (fmt == "lammps")
        write_lammps_data(aconf, filename);
    else if (fmt == "pdb")
        write_pdb(aconf, filename);
    else if (fmt == "xyza")
        write_xyza(aconf, filename);
    else
        assert(!"Unknown format");
}

//=============   AtomEye file reading and writing    ==================

void read_atomeye(LineInput& in, dbr_aconf *aconf, bool reduced_coords)
{
    // only "A = 1.0 [Angstrom]" files are supported
    char const* cfg_magic = "Number of particles =";
    char const* extended_mark = "entry_count";
    int const magic_len = strlen(cfg_magic);
    double H[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    const char* line = in.get_line();
    aconf->reduced_coordinates = reduced_coords;
    assert(!strncmp(line, cfg_magic, magic_len));
    sscanf(line + magic_len, "%i", &aconf->n);
    if (dbr_verbosity > 0)
        mcerr << aconf->n << " atoms." << endl;
    aconf->atoms = new dbr_atom[aconf->n];
    int counter = 0;
    bool extended = false;
    int entry_count = 0;
    char ext_name[8] = {0, 0, 0, 0, 0, 0, 0, 0};

    while ((line = in.get_line())) {
        char const* nonblank = line;
        while (isspace(*nonblank))
            ++nonblank;

        // PBC line, H0(...) = ...
        if (nonblank[0] == 'H' && nonblank[1] == '0') {
            char const* t = strchr(line, '(');
            char *end;
            int a1 = strtol(++t, &end, 10);
            assert (a1 >= 1 && a1 <= 3);

            t = strchr(end, ',');
            int a2 = strtol(++t, &end, 10);
            assert (a2 >= 1 && a2 <= 3);
            t = strchr(end, '=');
            H[a1-1][a2-1] = strtod(++t, 0);
        }

        // entry_count line
        else if (!strncmp(nonblank, extended_mark, strlen(extended_mark))) {
            const char *eq = strchr(nonblank, '=');
            ASSERT_FORMAT(*eq);
            int r = sscanf(eq+1, "%i", &entry_count);
            ASSERT_FORMAT(r == 1);
            if (dbr_verbosity > 0)
                mcerr << "Extended CFG format with " << entry_count
                    << " entries per atom." << endl;
            extended = true;
            if (entry_count > 3)
                aconf->auxiliary.resize(aconf->n);
        }

        else if (!strncmp(nonblank, "auxiliary", strlen("auxiliary"))) {
            aconf->auxiliary_header += nonblank + string("\n");
        }

        // either atom mass or coordinates line
        else if (nonblank[0] == '-' || isdigit(nonblank[0])) {
            dbr_atom& atom = aconf->atoms[counter];
            dbr_real x=0, y=0, z=0;
            if (extended) {
                int char_count;
                int r = sscanf(nonblank, DBR_F " " DBR_F " " DBR_F "%n",
                               &x, &y, &z, &char_count);
                // it is not known if "%n" increased r above or not
                if (r < 3) {
                    line = in.get_line();
                    sscanf(line, "%7s", ext_name);
                    continue;
                }
                // if we are here, it's a line with coordinates
                assert(r == 3 || r == 4);
                assert(strlen(ext_name) > 0);
                //{
                //    mcerr << "Error in line " << i+2 << endl;
                //    dbr_abort(EXIT_FAILURE);
                //}
                strcpy(atom.name, ext_name);
                if (!aconf->auxiliary.empty())
                    aconf->auxiliary[counter] = nonblank + char_count;
            }
            else {
                dbr_real mass;
                int r = sscanf(nonblank, DBR_F" %7s "DBR_F" "DBR_F" "DBR_F,
                                         &mass, atom.name, &x, &y, &z);
                assert(r == 5);
            }
            // wrap to <0,1)
            x -= floor(x);
            y -= floor(y);
            z -= floor(z);
            if (reduced_coords) {
                atom.xyz[0] = x;
                atom.xyz[1] = y;
                atom.xyz[2] = z;
            }
            else {
                atom.xyz[0] = x * H[0][0]  +  y * H[1][0]  +  z * H[2][0];
                atom.xyz[1] = x * H[0][1]  +  y * H[1][1]  +  z * H[2][1];
                atom.xyz[2] = x * H[0][2]  +  y * H[1][2]  +  z * H[2][2];
            }
            ++counter;
            if (counter == aconf->n)
                break;
        }

        // comment line
        else if (nonblank[0] == '#') {
            nonblank++;
            if (nonblank[0] == ' ')
                nonblank++;
            aconf->comments.push_back(nonblank);
        }

        else {
            //mcerr << "This line is not handled:\n" << line << endl;
        }
    }

    if (counter != aconf->n) {
        mcerr << "Error. Only " << counter << " atoms found in cfg file.\n";
        mcerr << in.get_error();
        dbr_abort(EXIT_FAILURE);
    }

    aconf->pbc.v00 = H[0][0];
    aconf->pbc.v01 = H[0][1];
    aconf->pbc.v02 = H[0][2];
    aconf->pbc.v10 = H[1][0];
    aconf->pbc.v11 = H[1][1];
    aconf->pbc.v12 = H[1][2];
    aconf->pbc.v20 = H[2][0];
    aconf->pbc.v21 = H[2][1];
    aconf->pbc.v22 = H[2][2];
}

void write_comments_with_hashes(dbr_aconf const& aconf, ostream &f)
{
    for (vector<string>::const_iterator i = aconf.comments.begin();
            i != aconf.comments.end(); ++i)
        f << "# " << *i << "\n";
}

int count_words(const char* s)
{
    int n = 0;
    while (isspace(*s))
        ++s;
    while (*s != '\0')
    {
        while (!isspace (*s) && *s != '\0')
            ++s;
        while (isspace(*s))
            ++s;
        n++;
    }
    return n;
}

void write_atoms_to_atomeye_file(dbr_aconf const& aconf, string const& filename)
{
    if (dbr_nid != 0)
        return;
    if (dbr_verbosity >= 0)
        mcerr << "Writing " << aconf.n << " atoms to AtomEye file: "
              << filename << endl;
    ofstream f(filename.c_str(), ios::out|ios::binary);
    if (!f) {
        mcerr << "Error. Can't open file " << filename << " for writing.\n";
        return;
    }
    f << "Number of particles =" << aconf.n << endl;
    write_comments_with_hashes(aconf, f);
    f << "\n";
    f << "A = 1.0 Angstrom (basic length-scale)\n";
    f << "H0(1,1) = " << setprecision(12) << aconf.pbc.v00 << " A\n";
    f << "H0(1,2) = " << setprecision(12) << aconf.pbc.v01 << " A\n";
    f << "H0(1,3) = " << setprecision(12) << aconf.pbc.v02 << " A\n";
    f << "H0(2,1) = " << setprecision(12) << aconf.pbc.v10 << " A\n";
    f << "H0(2,2) = " << setprecision(12) << aconf.pbc.v11 << " A\n";
    f << "H0(2,3) = " << setprecision(12) << aconf.pbc.v12 << " A\n";
    f << "H0(3,1) = " << setprecision(12) << aconf.pbc.v20 << " A\n";
    f << "H0(3,2) = " << setprecision(12) << aconf.pbc.v21 << " A\n";
    f << "H0(3,3) = " << setprecision(12) << aconf.pbc.v22 << " A\n";

    f << ".NO_VELOCITY.\n";
    int entry_count = 3;
    if (!aconf.auxiliary.empty())
        entry_count += count_words(aconf.auxiliary[0].c_str());
    f << "entry_count = " << entry_count << "\n";
    if (!aconf.auxiliary_header.empty())
        f << aconf.auxiliary_header << endl;
    double H_1[3][3]; //inverse of pbc matrix
    dbr_inverse_3x3_matrix(aconf.pbc, H_1);
    for (int i = 0; i < aconf.n; ++i) {
        dbr_atom const& atom = aconf.atoms[i];
        if (i == 0 || strcmp(atom.name, aconf.atoms[i-1].name) != 0) {
            const t_pse *pse = find_in_pse(atom.name);
            f << (pse ? pse->mass : 1.0) << "\n";
            f << atom.name << "\n";
        }
        if (aconf.reduced_coordinates) {
            f << atom.xyz[0] << " " << atom.xyz[1]
                << " " << atom.xyz[2];
        }
        else {
            dbr_xyz r;
            dbr_vec3_mult_mat3x3(atom.xyz, H_1, r);
            f << r[0] << " " << r[1] << " " << r[2];
        }
        if (!aconf.auxiliary.empty() && !aconf.auxiliary[i].empty())
            f << (isspace(aconf.auxiliary[i][0]) ? "" : " ")
              << aconf.auxiliary[i];
        f << "\n";
    }
    f.close();
}

//=============    DL_POLY file reading and writing    ==================


/* checks if "atom" is a shell in core-shell model */
static bool is_shell(const char *name)
{
    static const char* postfixes[] = { "-shell", "_shell", "-shel", "_shel",
                                        "-shl", "_shl", "-sh", "_sh" };
    const int n = sizeof(postfixes) / sizeof(postfixes[0]);
    for (int i = 0; i < n; ++i) {
        int diff = strlen(name) - strlen(postfixes[i]);
        if (diff > 0 && strcmp(postfixes[i], name+diff) == 0)
            return true;
    }
    return false;
}


void read_dlpoly_config(LineInput& in, dbr_aconf *aconf)
{
    const char* line = in.get_line(); // first line is a title
    line = in.get_line(); // second line - levcfg and imcon
    int levcfg, imcon;
    int r = sscanf(line, " %d %d", &levcfg, &imcon);
    ASSERT_FORMAT(r == 2);
    //  imcon -- periodic boundary key:
    //      imcon       meaning
    //      0   no periodic boundaries
    //      1   cubic boundary conditions
    //      2   orthorhombic boundary conditions
    //      3   parallelepiped boundary conditions
    //      4   truncated octahedral boundary conditions
    //      5   rhombic dodecahedral boundary conditions
    //      6   x-y parallelogram boundary conditions with
    //          no periodicity in the z direction
    //      7   hexagonal prism boundary conditions
    //
    // Possible values of levcfg: 0 - only coordinates in file
    //                            1 - coordinates and velocities
    //                            2 - coordinates, velocities and forces

    // read boundary conditions
    if (imcon == 0)
        ; //no PBC
    else if (imcon == 1 || imcon == 2 || imcon == 3) {
        dbr_pbc &p = aconf->pbc;

        int r = sscanf(in.get_line(), " %lf %lf %lf", &p.v00, &p.v01, &p.v02);
        ASSERT_FORMAT(r == 3);

        r = sscanf(in.get_line(), " %lf %lf %lf", &p.v10, &p.v11, &p.v12);
        ASSERT_FORMAT(r == 3);

        r = sscanf(in.get_line(), " %lf %lf %lf", &p.v20, &p.v21, &p.v22);
        ASSERT_FORMAT(r == 3);
    }
    else {
        mcerr << "Error. imcon=" << imcon << " is not supported." << endl;
        dbr_abort(EXIT_FAILURE);
    }

    // read atoms
    size_t counter = 0;
    size_t atoms_size = 1024;
    aconf->atoms = new dbr_atom[atoms_size];
    if (levcfg > 2)
        levcfg = 2;

    while (1) {
        if (!(line = in.get_line()))
            break;
        if (counter + 1 == atoms_size)
            resize_atoms(atoms_size, &aconf->atoms);
        dbr_atom& atom = aconf->atoms[counter];

        const char* nonblank = line;
        while (isspace(*nonblank))
            ++nonblank;

        char atom_name[8];
        atom_name[7] = 0;
        strncpy(atom_name, nonblank, 7);

        char* blank = atom_name;
        while (*blank && !isspace(*blank))
            ++blank;
        *blank = 0;

        if (is_shell(atom_name)) {
            for (int i = 0; i < levcfg+1; ++i)
                line = in.get_line();
            continue;
        }

        strncpy(atom.name, atom_name, 7);
        line = in.get_line();
        ASSERT_FORMAT(line);
        int r = sscanf(line, " " DBR_F " " DBR_F " " DBR_F,
                             &atom.xyz[0], &atom.xyz[1], &atom.xyz[2]);
        ASSERT_FORMAT(r == 3);
        for (int i = 1; i < levcfg+1; ++i)
            line = in.get_line();

        ++counter;
    }
    aconf->n = counter;
}


map<string,int> get_symbol_map(dbr_aconf const& aconf)
{
    map<string,int> m;
    for (int i = 0; i < aconf.n; ++i) {
        map<string,int>::iterator t = m.find(aconf.atoms[i].name);
        if (t == m.end())
            m[aconf.atoms[i].name] = 1;
        else
            t->second++;
    }
    return m;
}

void write_dlpoly_file(dbr_aconf const& aconf, string const& filename)
{
    if (dbr_nid != 0)
        return;
    if (dbr_verbosity >= 0)
        mcerr << "Writing atoms to DL_POLY CONFIG cfg file: "
              << filename << endl;
    // printf formatting is more useful here
    FILE *f = fopen(filename.c_str(), "wb");
    if (!f) {
        mcerr << "Error. Can't open file " << filename << " for writing.\n";
        return;
    }
    fprintf(f, "# converted by debyer\n");
    int imcon = 0;
    dbr_pbc const& pbc = aconf.pbc;
    if (pbc.v01 != 0 || pbc.v02 != 0 || pbc.v10 != 0 || pbc.v12 != 0
            || pbc.v20 != 0 || pbc.v21 != 0)
        imcon = 3;
    else if (pbc.v00 != pbc.v11 || pbc.v00 != pbc.v22)
        imcon = 2;
    else if (pbc.v00 != 0)
        imcon = 1;
    int levcfg = 0; //only coordinates in file
    fprintf(f, "%10d%10d\n", levcfg, imcon);
    if (imcon > 0) {
        fprintf(f, "%20.10f%20.10f%20.10f\n", pbc.v00, pbc.v01, pbc.v02);
        fprintf(f, "%20.10f%20.10f%20.10f\n", pbc.v10, pbc.v11, pbc.v12);
        fprintf(f, "%20.10f%20.10f%20.10f\n", pbc.v20, pbc.v21, pbc.v22);
    }
    for (int i = 0; i < aconf.n; ++i) {
        dbr_atom const& a = aconf.atoms[i];
        fprintf(f, "%s\n", a.name);
        fprintf(f, "%20.8f%20.8f%20.8f\n", a.xyz[0], a.xyz[1], a.xyz[2]);
    }
    fclose(f);
}

/// splits string into tokens, separated by white characters
vector<string> split_string(const char* s) {
    vector<std::string> v;
    const char* p = s;
    while (*p) {
        while (isspace(*p))
            ++p;
        const char* start_pos = p;
        while (*p && !isspace(*p))
            ++p;
        v.push_back(string(start_pos, p - start_pos));
    }
    return v;
}

// This file format doesn't contain atom names, only numbers of atom types.
// The names corresponding to the numbers are in separate lammps file.
// Here we assume that names of the types follow the line "n atom types"
// as comments, e.g.: "2 atom types # C Si", or use fake names.
void read_lammps_data(LineInput& in, dbr_aconf* aconf, bool reduced_coords)
{
    aconf->reduced_coordinates = reduced_coords;
    // read header
    vector<string> symbols;
    char* line = NULL;
    int counter = 0;
    while ((line = in.get_line())) {
        ++counter;
        char* nonblank = line;
        while (isspace(*nonblank))
            ++nonblank;

        if (*nonblank == '\0') // blank line
            continue;

        // comment line
        if (nonblank[0] == '#') {
            nonblank++;
            if (nonblank[0] == ' ')
                nonblank++;
            aconf->comments.push_back(nonblank);
            continue;
        }

        char *trailing_comment = strchr(nonblank, '#');
        if (trailing_comment) {
            *trailing_comment = '\0';
        }

        if (strstr(nonblank, "atoms")) {
            aconf->n = strtol(nonblank, NULL, 10);
        }
        else if (strstr(nonblank, "atom types")) {
            int ntypes = strtol(nonblank, NULL, 10);
            if (trailing_comment) {
                symbols = split_string(trailing_comment+1);
                if ((int) symbols.size() != ntypes) {
                    mcerr << "Trailing comment in `atom types' line should "
                        " contain " << ntypes << " symbols rather than " <<
                        symbols.size() << ".\n" << trailing_comment+1 << endl;
                        dbr_abort(EXIT_FAILURE);
                }
            }
            else {
                mcerr << "Warning: No atomic symbols. "
                         "You may add the symbols as a comment after "
                         "the `atom types' line (e.g. '# Si C')" << endl;
                for (int i = 0; i < ntypes; ++i) {
                    char sym[20];
                    snprintf(sym, 20, "A%d", i);
                    symbols.push_back(sym);
                }
            }
        }
        else if (strstr(nonblank, "xlo xhi")) {
            double hi, lo;
            if (sscanf(nonblank, "%lf %lf", &lo, &hi) != 2) {
                mcerr << "Error in line:\n" << line << endl;
                dbr_abort(EXIT_FAILURE);
            }
            aconf->pbc.v00 = hi - lo;
        }
        else if (strstr(nonblank, "ylo yhi")) {
            double hi, lo;
            if (sscanf(nonblank, "%lf %lf", &lo, &hi) != 2) {
                mcerr << "Error in line:\n" << line << endl;
                dbr_abort(EXIT_FAILURE);
            }
            aconf->pbc.v11 = hi - lo;
        }
        else if (strstr(nonblank, "zlo zhi")) {
            double hi, lo;
            if (sscanf(nonblank, "%lf %lf", &lo, &hi) != 2) {
                mcerr << "Error in line:\n" << line << endl;
                dbr_abort(EXIT_FAILURE);
            }
            aconf->pbc.v22 = hi - lo;
        }
        else if (strstr(nonblank, "Atoms") == nonblank) {
            break;
        }
        else {
            mcerr << "Warning: ignoring line " << counter << ":\n" << line
                << endl;
        }
    }

    if (aconf->n == 0) {
        mcerr << "Error: line with number of atoms not found.\n";
        dbr_abort(EXIT_FAILURE);
    }
    // read data
    aconf->atoms = new dbr_atom[aconf->n];
    for (int i = 0; i < aconf->n; ++i)
        aconf->atoms[i].name[0] = 0;

    int atom_counter = 0;
    while ((line = in.get_line())) {
        const char *nonblank = line;
        while (isspace(*nonblank))
            ++nonblank;
        if (*nonblank == '\0') // blank line
            continue;
        int number, symbol_nr;
        double ax, ay, az;
        if (sscanf(line, "%d %d %lf %lf %lf",
                         &number, &symbol_nr, &ax, &ay, &az) != 5) {
            mcerr << "Format error in line:\n" << line << endl;
            dbr_abort(EXIT_FAILURE);
        }
        if (number < 1 || number > aconf->n) {
            mcerr << "Wrong atom number in line:\n" << line << endl;
            dbr_abort(EXIT_FAILURE);
        }
        if (symbol_nr < 1 || symbol_nr > (int) symbols.size()) {
            mcerr << "Wrong type number in line:\n" << line << endl;
            dbr_abort(EXIT_FAILURE);
        }
        dbr_atom *a = &aconf->atoms[number-1];
        if (a->name[0] != 0) {
            mcerr << "Atom number " << number << " was duplicated.\n";
            dbr_abort(EXIT_FAILURE);
        }
        strncpy(a->name, symbols[symbol_nr-1].c_str(), 7);
        if (reduced_coords) {
            ax /= aconf->pbc.v00;
            ay /= aconf->pbc.v11;
            az /= aconf->pbc.v22;
            // wrap to <0,1)
            ax -= floor(ax);
            ay -= floor(ay);
            az -= floor(az);
        }
        a->xyz[0] = ax;
        a->xyz[1] = ay;
        a->xyz[2] = az;
        ++atom_counter;
        if (atom_counter == aconf->n)
            return;
    }
    mcerr << "File was terminated before all the atoms were read\n";
    dbr_abort(EXIT_FAILURE);
    return;
}

void write_lammps_data(dbr_aconf const& aconf, string const& filename)
{
    if (dbr_nid != 0)
        return;
    if (aconf.pbc.v01 != 0 || aconf.pbc.v02 != 0
          || aconf.pbc.v10 != 0 || aconf.pbc.v12 != 0
          || aconf.pbc.v20 != 0 || aconf.pbc.v21 != 0) {
        mcerr << "Sorry, LAMMPS files with non-orthorhombic PBC "
                 "are not supported\n";
        return;
    }
    if (dbr_verbosity >= 0)
        mcerr << "Writing atoms to LAMMPS data file: " << filename << endl;
    ofstream f(filename.c_str(), ios::out|ios::binary);
    if (!f) {
        mcerr << "Error. Can't open file " << filename << " for writing.\n";
        return;
    }
    write_comments_with_hashes(aconf, f);
    f << "\n";
    f << "#units           real #or metal\n";
    f << "#boundary        p p p\n";
    f << "#atom_style      atomic\n";
    f << "#read_data       " << filename << "\n";
    map<string,int> m = get_symbol_map(aconf);
    vector<string> keys;
    for (map<string,int>::const_iterator i = m.begin(); i != m.end(); ++i)
        keys.push_back(i->first);
    sort(keys.begin(), keys.end());

    for (size_t i = 0; i != keys.size(); ++i) {
        const t_pse *pse = find_in_pse(keys[i].c_str());
        f << "#mass            " << i+1 << " ";
        if (pse)
            f << pse->mass;
        else
            f << "?";
        f << "  #" << keys[i] << endl;
    }
    f << "\n" << aconf.n << "\tatoms\n";
    f << keys.size() << " atom types #";
    for (size_t i = 0; i != keys.size(); ++i) {
        f << " " << keys[i];
        m[keys[i]] = i+1;
    }
    f << "\n";
    f << 0 << " " << setprecision(12) << aconf.pbc.v00 << " xlo xhi\n";
    f << 0 << " " << setprecision(12) << aconf.pbc.v11 << " ylo yhi\n";
    f << 0 << " " << setprecision(12) << aconf.pbc.v22 << " zlo zhi\n";
    f << "\nAtoms\n\n";
    for (int i = 0; i < aconf.n; ++i) {
        dbr_atom const& a = aconf.atoms[i];
        f << i+1 << "\t" << m[a.name] << "\t"
          << a.xyz[0] << "\t" << a.xyz[1] << "\t" << a.xyz[2] << "\n";
    }
    f.close();
}

void write_pdb(dbr_aconf const& aconf, string const& filename)
{
    if (dbr_nid != 0)
        return;
    if (dbr_verbosity >= 0)
        mcerr << "Writing atoms to PBC data file: " << filename << endl;
    FILE *f = fopen(filename.c_str(), "wb");
    if (f == NULL) {
        mcerr << "Error. Can't open file " << filename << " for writing.\n";
        return;
    }

    fprintf(f,"HEADER    converted by debyer; %10d atoms"
              "                                 \n", aconf.n);

    // "The CRYST1 record presents the unit cell parameters, space group,
    // and Z value. If the structure was not determined by crystallographic
    // means, CRYST1 simply defines a unit cube."
    //
    // COLUMNS       DATA TYPE      FIELD         DEFINITION
    // -------------------------------------------------------------
    //   1 -  6       Record name    "CRYST1"
    //   7 - 15       Real(9.3)      a             a (Angstroms).
    //   16 - 24      Real(9.3)      b             b (Angstroms).
    //   25 - 33      Real(9.3)      c             c (Angstroms).
    //   34 - 40      Real(7.2)      alpha         alpha (degrees).
    //   41 - 47      Real(7.2)      beta          beta (degrees).
    //   48 - 54      Real(7.2)      gamma         gamma (degrees).
    //   56 - 66      LString        sGroup        Space group.
    //   67 - 70      Integer        z             Z value.
    //
    //   here CRYST1 is used for PBC, space group is set as P1
    dbr_pbc_prop p = get_pbc_properties(aconf.pbc);
    fprintf(f, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f"
               " P 1          1           \n",
               p.lengths[0], p.lengths[1], p.lengths[2],
               acos(p.cosines[0]) * 180 / M_PI,
               acos(p.cosines[1]) * 180 / M_PI,
               acos(p.cosines[2]) * 180 / M_PI);

    // "The ATOM records present the atomic coordinates for standard residues."
    // COLUMNS        DATA TYPE       FIELD         DEFINITION
    // --------------------------------------------------------------
    // 1 -6     Record name    "ATOM  "
    // 7 -11    Integer        serial      Atom serial number.
    // 13-16    Atom           name        Atom name.
    // 17       Character      altLoc      Alternate location indicator.
    // 18-20    Residue name   resName     Residue name.
    // 22       Character      chainID     Chain identifier.
    // 23-26    Integer        resSeq      Residue sequence number.
    // 27       AChar          iCode       Code for insertion of residues.
    // 31-38    Real(8.3)      x           Orthogonal coordinates for X in A
    // 39-46    Real(8.3)      y           Orthogonal coordinates for Y in A
    // 47-54    Real(8.3)      z           Orthogonal coordinates for Z in A
    // 55-60    Real(6.2)      occupancy   Occupancy.
    // 61-66    Real(6.2)      tempFactor  Temperature factor.
    // 73-76    LString(4)     segID       Segment identifier, left-justified.
    // 77-78    LString(2)     element     Element symbol, right-justified.
    // 79-80    LString(2)     charge      Charge on the atom.
    for (int i = 0; i < aconf.n; ++i) {
        const dbr_atom& atom = aconf.atoms[i];
        fprintf(f, "ATOM  %5d %-3s          1    %8.3f%8.3f%8.3f"
                   "                      %2s  \n",
                i+1, atom.name, atom.xyz[0], atom.xyz[1], atom.xyz[2],
                atom.name);
    }
    fprintf(f, "END                                                        "
           "                     \n");
    fclose(f);
}

void write_xyza(dbr_aconf const& aconf, string const& filename)
{
    if (dbr_nid != 0)
        return;
    if (dbr_verbosity >= 0)
        mcerr << "Writing atoms to XYZA file: " << filename << endl;
    ofstream f(filename.c_str(), ios::out|ios::binary);
    if (!f) {
        mcerr << "Error. Can't open file " << filename << " for writing.\n";
        return;
    }
    if (aconf.pbc.v00 != 0)
        mcerr << "Warning: PBC is not stored in the file.\n"
                  "PBC: " << aconf.pbc.v00 << " x " << aconf.pbc.v11
                  << " x " << aconf.pbc.v22 << endl;

    for (int i = 0; i < aconf.n; ++i) {
        dbr_atom const& atom = aconf.atoms[i];
        f << atom.xyz[0] << "\t" << atom.xyz[1] << "\t" << atom.xyz[2] << "\t"
            << atom.name << "\t" << endl;
    }
    f.close();
}

// for using in C programs
int read_atoms_c(const char* filename, dbr_atom **atoms, dbr_pbc *pbc,
                 int reduced)
{
    LineInput in;
    bool r = in.init(filename);
    if (!r) {
        mcerr << in.get_error();
        dbr_abort(EXIT_FAILURE);
    }
    dbr_aconf aconf = read_atoms_from_file(in, (bool)reduced, "");
    if (atoms != NULL)
        *atoms = aconf.atoms;
    if (pbc != NULL)
        *pbc = aconf.pbc;
    return aconf.n;
}

void free_atoms(dbr_atom *atoms)
{
    delete[] atoms;
}

SimpleLineInput::SimpleLineInput(const char* filename)
{
    bool ok = init(filename);
    if (!ok) {
        mcerr << get_error();
        dbr_abort(EXIT_FAILURE);
    }
}


