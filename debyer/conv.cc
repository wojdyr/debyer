//  dbr_conv -- convert between file formats with atomistic models
//  Copyright 2010 Marcin Wojdyr
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

#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>

#include "fileio.h"

using namespace std;

namespace {

const char* usage[] = {
"Usage: dbr_conv [OPTIONS]... INPUT_FILE OUTPUT_FILE",
"       dbr_conv [OPTIONS...] -t TO_FORMAT -m INPUT_FILE1 [INPUT_FILE2]...",
"",
"  -h          Print help and exit",
"  -V          Print version and exit",
"  -q          Silent mode",
"  -v          Increase verbosity level (can be used twice)",
"  -s          Sort atoms by atomic symbols.",
"  -f FORMAT   Convert file from format.",
"  -t FORMAT   Convert file to format.",
"  -m          Convert multiple files. Output files have only file extension",
"              changed.",
"",
"Supported format names: atomeye, dlpoly, lammps, pdb, xyz, xyza.",
"Compressed files (.gz, .bz2) can be read, but not written.",
"PDB format can only be written."
};

void print_usage()
{
    int len = sizeof(usage)/sizeof(usage[0]);
    for (int i = 0; i < len; ++i)
        printf("%s\n", usage[i]);
}

void print_version()
{
#ifndef VERSION
# define VERSION "unknown"
#endif
    printf("dbr_conv version %s\n", VERSION);
}

string get_extension_for_format(const string& format)
{
    if (format == "atomeye")
        return "cfg";
    else if (format == "dlpoly")
        return "CONFIG";
    else if (format == "lammps")
        return "lammps";
    else if (format == "xyz")
        return "xyz";
    else if (format == "xyza")
        return "xyza";
    else
        return "";
}

string change_extension(const string& filename, const string& new_ext)
{
    size_t p = filename.rfind('.');
    if (p == 0 || p == string::npos)
        return filename + "." + new_ext;
    string old_ext = filename.substr(p);
    if (old_ext == ".gz" || old_ext == ".bz2") {
        size_t prev = filename.rfind('.', p-1);
        if (prev != 0 && prev != string::npos)
            p = prev;
    }
    return filename.substr(0, p+1) + new_ext;
}

struct atom_compare: public binary_function<int, int, bool>
{
    atom_compare(const dbr_atom* atoms) : atoms_(atoms) {}
    bool operator() (int a, int b) const
            { return strcmp(atoms_[a].name, atoms_[b].name) < 0; }
    const dbr_atom* atoms_;
};

void sort_atoms(dbr_aconf& aconf)
{
    vector<int> indices(aconf.n);
    for (int i = 0; i != aconf.n; ++i)
        indices[i] = i;
    sort(indices.begin(), indices.end(), atom_compare(aconf.atoms));
    dbr_atom* orig_atoms = new dbr_atom[aconf.n];
    memcpy(orig_atoms, aconf.atoms, aconf.n * sizeof(dbr_atom));
    vector<string> orig_auxiliary = aconf.auxiliary;
    for (int i = 0; i != aconf.n; ++i) {
        int new_pos = indices[i];
        aconf.atoms[new_pos] = orig_atoms[i];
        aconf.auxiliary[new_pos] = orig_auxiliary[i];
    }
    delete [] orig_atoms;
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    dbr_init(&argc, &argv); // set global variables

    // parse command line
    if (argc < 2) {
        print_usage();
        return EXIT_SUCCESS;
    }
    bool option_sort = false;
    bool option_multi = false;
    string option_from_format;
    string option_to_format;
    vector<string> files;
    for (int i = 1; i < argc; ++i) {
        const char* arg = argv[i];
        int arg_len = strlen(arg);
        if (arg_len > 1 && arg[0] == '-') {
            switch (arg[1]) {
                case 'h':
                    print_usage();
                    return EXIT_SUCCESS;
                case 'V':
                    print_version();
                    return EXIT_SUCCESS;
                case 'v':
                    ++dbr_verbosity;
                    break;
                case 'q':
                    --dbr_verbosity;
                    break;
                case 's':
                    option_sort = true;
                    break;
                case 'm':
                    option_multi = true;
                    break;
                case 'f':
                    if (arg_len > 2)
                        option_from_format = arg+2;
                    else {
                        ++i;
                        option_from_format = argv[i];
                    }
                    break;
                case 't':
                    if (arg_len > 2)
                        option_to_format = arg+2;
                    else {
                        ++i;
                        option_to_format = argv[i];
                    }
                    break;
                default:
                    fprintf(stderr, "Error: unexpected option '-%c'", arg[1]);
                    fprintf(stderr, "Run dbr_conv -h to print help.\n");
                    return EXIT_FAILURE;
            }
        }
        else
            files.push_back(argv[i]);
    }

    // convert file(s)
    bool reduced_coords = false;
    if (option_multi) {
        if (option_to_format.empty()) {
            fprintf(stderr, "Error: Option -m without -t.\n");
            fprintf(stderr, "Run dbr_conv -h to print help.\n");
            return EXIT_FAILURE;
        }
        for (size_t i = 0; i != files.size(); ++i) {
            string input_file = files[i];
            dbr_aconf aconf = read_atoms_from_file(input_file.c_str(),
                                                   reduced_coords,
                                                   option_from_format);
            if (option_sort)
                sort_atoms(aconf);
            string ext = get_extension_for_format(option_to_format);
            string output_file = change_extension(input_file, ext);
            if (input_file == output_file) {
                fprintf(stderr, "Refusing to overwrite input file.\n");
                return EXIT_FAILURE;
            }
            write_file_with_atoms(aconf, output_file, option_to_format);
            delete [] aconf.atoms;
        }
    }
    else {
        int n = files.size();
        if (n != 2) {
            fprintf(stderr, "Error: expected 2 filenames (or -m), got %d\n", n);
            fprintf(stderr, "Run dbr_conv -h to print help.\n");
            return EXIT_FAILURE;
        }
        dbr_aconf aconf = read_atoms_from_file(files[0].c_str(), reduced_coords,
                                               option_from_format);
        if (option_sort)
            sort_atoms(aconf);
        write_file_with_atoms(aconf, files[1], option_to_format);
    }
    return EXIT_SUCCESS;
}


