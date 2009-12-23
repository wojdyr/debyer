//  debyer -- program for calculation of diffration patterns
//  Copyright (C) 2006 Marcin Wojdyr
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
//
//  $Id$
//
//  Main loop of debyer program. Parses command-line arguments (using
//  gengetopt and debyer.ggo) and uses libdebyer for actual computations.

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <cstring>
#include "cmdline.h"
#include "debyer.h"
#include "fileio.h"
#include "lineio.h"

using namespace std;

#define mcerr \
    if (dbr_nid == 0) \
        cerr

output_kind get_output_from_args (gengetopt_args_info const& args)
{
    if (args.RDF_given)
        return output_rdf;
    else if (args.PDF_given)
        return output_pdf;
    else if (args.rPDF_given)
        return output_rpdf;
    else if (args.neutron_given)
        return output_neutron;
    else if (args.xray_given)
        return output_xray;
    else if (args.sf_given)
        return output_sf;
    else
        return output_none;
}

string default_fn(gengetopt_args_info const& args, const char* ext)
{
    char const* infn = args.inputs[0];
    char const* dot = strrchr(infn, '.');

    if (dot && dot > infn
            && (strcmp(dot, ".gz") == 0 || strcmp(dot, ".bz2") == 0)) {
        --dot;
        while (*dot != '.') {
            if (dot == infn) {
                dot = NULL;
                break;
            }
            --dot;
        }
    }

    if (strcmp(dot+1, ext) == 0)
        return infn + string(".new");
    string base = dot ? string(infn, dot) : string(infn);
    return base + "." + ext;
}


void print_irdfs_statistics(irdfs rdfs)
{
    // show i-rdf statistics
#ifdef __STRICT_ANSI__
    typedef long t_verylong; //may be insufficient
#else
    typedef long long t_verylong;
#endif //__STRICT_ANSI__
    t_verylong sum = 0;
    for (int i = 0; i < rdfs.pair_count; ++i) {
        for (int j = 0; j < rdfs.rdf_bins; ++j)
            sum += rdfs.data[i].nn[j];
    }
    t_verylong ac = 0;
    for (int i = 0; i < rdfs.symbol_count; ++i)
        ac += rdfs.atom_counts[i];
    if (dbr_verbosity >= 0)
        mcerr << sum << " of " << ac * (ac-1) / 2
            << " pairs included in internal PDF.\n";
}

void print_min_dist_info(irdfs const& rdfs)
{
    mcerr << "Min. interatomic dist. (+/- "  << rdfs.step / 2. << "):";
    for (int i = 0; i < rdfs.pair_count; ++i) {
        irdf const& rdf = rdfs.data[i];
        mcerr << " " << rdf.at1 << "-" << rdf.at2 << ":";
        for (int j = 0; j < rdfs.rdf_bins; ++j) {
            if (rdf.nn[j] != 0) {
                mcerr << (j+0.5) * rdfs.step;
                break;
            }
        }
    }
    mcerr << endl;
}

// get name of file to write ID in
const char* get_id_out_fn(gengetopt_args_info const& args)
{
    if (args.id_file_given) {
        if (args.id_file_arg)
            return args.id_file_arg;
        else {
            string default_id_file = default_fn(args, "id");
            char *allocated = new char[default_id_file.size() + 1];
            strcpy(allocated, default_id_file.c_str());
            return allocated;
        }
    }
    else
        return 0;
}

dbr_picker prepare_picker(gengetopt_args_info const& args, int n_atoms)
{
    dbr_picker picker;
    picker.probab = args.sample_given ? (dbr_real)args.sample_arg / n_atoms
                                      : 0.;
    if (args.x_gt_given || args.x_lt_given
            || args.y_gt_given || args.y_lt_given
            || args.z_gt_given || args.z_lt_given)
        picker.cut = 1;
    else
        picker.cut = 0;
    picker.x_min = args.x_gt_given ? args.x_gt_arg : -DBL_MAX;
    picker.x_max = args.x_lt_given ? args.x_lt_arg : +DBL_MAX;
    picker.y_min = args.y_gt_given ? args.y_gt_arg : -DBL_MAX;
    picker.y_max = args.y_lt_given ? args.y_lt_arg : +DBL_MAX;
    picker.z_min = args.z_gt_given ? args.z_gt_arg : -DBL_MAX;
    picker.z_max = args.z_lt_given ? args.z_lt_arg : +DBL_MAX;
    picker.all = (!picker.cut && picker.probab == 0.);
    return picker;
}

void write_atoms_to_file(dbr_aconf const& aconf,
                         gengetopt_args_info const& args)
{
    // format conversion
    if (args.write_xyz_given)
        write_atoms_to_xyz_file(aconf, args.write_xyz_arg
                                       ? args.write_xyz_arg
                                       : default_fn(args, "xyz"));
    if (args.write_cfg_given)
        write_atoms_to_atomeye_file(aconf, args.write_cfg_arg
                                           ? args.write_cfg_arg
                                           : default_fn(args, "cfg"));
    if (args.write_dlpoly_given)
        write_dlpoly_file(aconf, args.write_dlpoly_arg
                                 ? args.write_dlpoly_arg
                                 : default_fn(args, "dlpoly"),
                          false);
    if (args.write_dlpoly_s_given)
        write_dlpoly_file(aconf, args.write_dlpoly_s_arg
                                 ? args.write_dlpoly_s_arg
                                 : default_fn(args, "dlpoly"),
                          true);
    if (args.write_lammps_data_given)
        write_lammps_data(aconf, args.write_lammps_data_arg
                                 ? args.write_lammps_data_arg
                                 : default_fn(args, "lammps"));
    if (args.write_pdb_given)
        write_pdb(aconf, args.write_pdb_arg
                         ? args.write_pdb_arg
                         : default_fn(args, "pdb"));
    if (args.write_xyza_given)
        write_xyza(aconf, args.write_xyza_arg
                         ? args.write_xyza_arg
                         : default_fn(args, "xyza"));

    // is that all to do?
    if (args.output_group_counter == 0 && !args.id_file_given)
        dbr_abort(EXIT_SUCCESS);
}

// atoms array in the dbr_aconf needs to be delete'd [] later
dbr_aconf prepare_aconf(LineInput &in, gengetopt_args_info const& args)
{
    dbr_aconf aconf = read_atoms_from_file(in, false);
    if (args.pbc_a_given)
        aconf.pbc.v00 = args.pbc_a_arg;
    if (args.pbc_b_given)
        aconf.pbc.v11 = args.pbc_b_arg;
    if (args.pbc_c_given)
        aconf.pbc.v22 = args.pbc_c_arg;
    return aconf;
}


irdfs calculate_id_from_datafile(dbr_aconf &aconf, dbr_picker const& picker,
                                 gengetopt_args_info const& args)
{
    // read atoms
    irdfs rdfs;
    // separete different atom types (dbr_atom[] -> dbr_atoms[])
    dbr_atoms *xa = 0;
    int tc = dbr_get_atoms(aconf.n, aconf.atoms, &xa, 0);
    delete [] aconf.atoms;

    if (dbr_verbosity > 0)
        mcerr << "Elapsed " << dbr_get_elapsed() <<" s. Atoms were sorted."
            << endl;
    if (dbr_verbosity >= 0) {
        mcerr << "Atom statistics:";
        for (int i = 0; i < tc; ++i)
            mcerr << " " << xa[i].name << ":" << xa[i].count;
        mcerr << endl;
    }

    // if we calculate RDF or (r)PDF, quanta and cutoff don't need to be given
    // explicitely
    bool is_direct = dbr_is_direct(get_output_from_args(args));

    dbr_real quanta = args.quanta_arg;
    if (!args.id_file_given && is_direct && args.step_given) {
        quanta = args.step_arg;
        if (dbr_verbosity > 0)
            mcerr << "Option `quanta' set to " << quanta << endl;
    }

    dbr_real cutoff = args.cutoff_given ? args.cutoff_arg : 0.;
    if (!args.cutoff_given && is_direct && args.to_given) {
        cutoff = args.to_arg;
        if (dbr_verbosity > 0)
            mcerr << "Option `cutoff' set to " << cutoff << endl;
    }

    //calculate internal PDFs - most of computer time is spent here
    rdfs = calculate_irdfs(tc, xa, cutoff,
                           quanta, aconf.pbc, &picker, get_id_out_fn(args));
    print_irdfs_statistics(rdfs);
    return rdfs;
}


int do_benchmark1(int count)
{
    srand(12345);
    if (dbr_verbosity >= 0)
        mcerr << "Performing ID calculation benchmark with " << count
            << " atoms." << endl;
    dbr_atom* coords = new dbr_atom[count];
    dbr_real a = pow(count / 0.095, 1./3.);
    for (int i = 0; i < count; ++i) {
        strcpy(coords[i].name, i % 2 ? "Si" : "C");
        for (int j = 0; j < 3; ++j)
        coords[i].xyz[j] = a * rand() / RAND_MAX;
    }
    dbr_atoms *xa = 0;
    int tc = dbr_get_atoms(count, coords, &xa, 0);
    delete [] coords;
    if (dbr_verbosity > 0)
        mcerr << "SiC pseudo-grain (cube a=" << a << ") was generated." << endl;
    dbr_pbc pbc;
    pbc.v00 = pbc.v11 = pbc.v22 = 0.;
    dbr_picker picker;
    picker.all = 1;
    picker.cut = 0;
    picker.probab = 0.;
    irdfs rdfs = calculate_irdfs(tc, xa, 0., 0.001, pbc, &picker, 0);
    print_irdfs_statistics(rdfs);
    if (dbr_verbosity > 0)
        mcerr << "Benchmark finished." << endl;
    return 0;
}

string get_output_filename(gengetopt_args_info const& args)
{
    if (args.output_given)
        return strcmp(args.output_arg, "stdout") != 0 ? args.output_arg : "";
    output_kind kind = get_output_from_args(args);
    if (kind == output_none)
        return "";

    const char *ext = 0;
    if (kind == output_rdf)
        ext = "rdf";
    else if (kind == output_pdf)
        ext = "g";
    else if (kind == output_rpdf)
        ext = "rpdf";
    else if (kind == output_neutron)
        ext = args.lambda_given ? "n" : "nq";
    else if (kind == output_xray)
        ext = args.lambda_given ? "x" : "xq";
    else if (kind == output_sf)
        ext = args.lambda_given ? "s" : "sq";
    else
        assert(0);
    return default_fn(args, ext);
}

int main(int argc, char **argv)
{
    int r;
    gengetopt_args_info args;
    dbr_init(&argc, &argv);
    if (cmdline_parser(argc, argv, &args) != 0)
        dbr_abort(EXIT_FAILURE);
    if (args.version_given) {
        dbr_print_version();
        dbr_abort(EXIT_SUCCESS);
    }
    dbr_verbosity = args.verbose_given - args.quiet_given;
    if (dbr_verbosity > 0)
        mcerr << "Verbosity level: " << dbr_verbosity << endl;
    if (args.ben1_given) {
        return do_benchmark1(args.ben1_arg);
    }
    //if (args.inputs_num > 1) {
    //    mcerr << "Only one input file can be given.";
    if (args.inputs_num != 1) {
        mcerr << "Exactly one input file should be given." << endl;
        dbr_abort(EXIT_FAILURE);
    }
    bool output_atoms = (args.write_xyz_given || args.write_cfg_given
            || args.write_dlpoly_given || args.write_dlpoly_s_given
            || args.write_lammps_data_given || args.write_pdb_given
            || args.write_xyza_given);
    if (args.output_group_counter == 0 && !args.id_file_given && !output_atoms){
        mcerr << "No output format specified, nothing to do." << endl;
        dbr_abort(EXIT_SUCCESS);
    }
    if (args.inputs_num == 0) {
        mcerr << "Reading from stdin .... TODO" << endl;
        dbr_abort(EXIT_FAILURE);
    }

    char *infn = args.inputs[0];
    if (dbr_verbosity >= 0)
        mcerr << "examining file: " << infn << endl;

    SimpleLineInput in(infn);

    output_kind kind = get_output_from_args(args);
    string ofname = get_output_filename(args);

    const char *id_magic = "debyer-id ";
    bool input_with_id = !strncmp(in.get_buffer(), id_magic, strlen(id_magic));
    irdfs rdfs;
    if (input_with_id) {
        if (dbr_verbosity >= 0)
            mcerr << "ID file detected, skipping ID calculation phase." << endl;
        if (args.pbc_a_given || args.pbc_b_given || args.pbc_c_given) {
            mcerr << "Error. PBC options can't be applied to ID file." << endl;
            dbr_abort(EXIT_FAILURE);
        }
        if (args.quanta_given) {
            mcerr << "Error. --quanta option can't be applied to ID input file."
                << endl;
            dbr_abort(EXIT_FAILURE);
        }
        if (output_atoms) {
            mcerr << "Error. Can't deduce atom coordinates from ID file."
                << endl;
            dbr_abort(EXIT_FAILURE);
        }
        rdfs = read_irdfs_from_file(infn);
        if (!rdfs.data)
            dbr_abort(EXIT_FAILURE);
        if (args.id_file_given) {
            mcerr << "Refusing to write ID file when input is also ID file. "
                     "Just copy the file." << endl;
            //write_irdfs_to_file(rdfs, get_id_out_fn(args));
        }
    }
    else { //need to calculate id
        dbr_aconf aconf = prepare_aconf(in, args);
        dbr_picker picker = prepare_picker(args, aconf.n);

        if (output_atoms) {
            aconf.comments.insert(aconf.comments.begin(),
                                  argv_as_str(argc, argv));
            if (picker.cut) {
                // write to file only selected atoms
                dbr_atom *all = aconf.atoms;
                dbr_atom *sel = new dbr_atom[aconf.n];
                int n_all = aconf.n;
                int j = 0;
                for (int i = 0; i < n_all; ++i)
                    if (dbr_is_atom_in_sector(aconf.atoms[i].xyz, &picker)) {
                        memcpy(&sel[j], &all[i], sizeof(dbr_atom));
                        ++j;
                    }
                aconf.atoms = sel;
                aconf.n = j;
                write_atoms_to_file(aconf, args);
                aconf.atoms = all;
                aconf.n = n_all;
                delete [] sel;
            }
            else
                write_atoms_to_file(aconf, args);
        }

        rdfs = calculate_id_from_datafile(aconf, picker, args);
    }
    if (dbr_verbosity > 0)
        mcerr << "Elapsed " << dbr_get_elapsed() <<" s. ID is ready." << endl;
    if (dbr_verbosity > 0)
        print_min_dist_info(rdfs);

    const char *type_desc = 0;
    if (kind == output_rdf)
        type_desc = "RDF";
    else if (kind == output_pdf)
        type_desc = "PDF";
    else if (kind == output_rpdf)
        type_desc = "reduced PDF";
    else if (kind == output_neutron)
        type_desc = "neutron diffraction pattern";
    else if (kind == output_xray)
        type_desc = "x-ray diffraction pattern";
    else if (kind == output_sf)
        type_desc = "structure function";
    else { // kind == output_none
        free_irdfs(&rdfs);
        dbr_finalize();
        return 0;
    }

    if (dbr_is_inverse(kind) && args.cutoff_given
            && !(args.ro_given || rdfs.density > 0))
        if (dbr_verbosity >= 0)
            mcerr << "No number density - no cut-off correction.\n";
    if (dbr_verbosity >= 0)
        mcerr << "Writing " << type_desc << " to file " << ofname << " ... ";
    if (dbr_is_direct(kind))
        r = write_pdfkind_to_file(kind,
                                  rdfs,
                                  args.from_given ? args.from_arg : 0.,
                                  args.to_given ? args.to_arg : 0.,
                                  args.step_given ? args.step_arg : 0.,
                                  args.ro_given ? args.ro_arg : 0.,
                                  args.weight_arg[0],
                                  ofname.c_str());
    else {
        r = write_diffraction_to_file(kind,
                                      rdfs,
                                      args.from_given ? args.from_arg : 0.,
                                      args.to_given ? args.to_arg : 0.,
                                      args.step_given ? args.step_arg : 0.,
                                      args.lambda_given ? args.lambda_arg : 0.,
                                      args.ro_given ? args.ro_arg : 0.,
                                      args.cutoff_given ? args.cutoff_arg : 0.,
                                      ofname.c_str());
    }
    free_irdfs(&rdfs);
    if (dbr_verbosity >= 0)
        mcerr << (r == 0 ? "OK" : "Failed") << endl;
    dbr_finalize();
    cmdline_parser_free(&args);
    if (dbr_verbosity > 0)
        mcerr << "Elapsed " << dbr_get_elapsed() <<" s. Finished." << endl;
    return 0;
}

