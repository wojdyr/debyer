/*  debyer -- program for calculation of diffration patterns
 *  Copyright (C) 2006-2007 Marcin Wojdyr
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  $Id: fileio.h 107 2009-04-12 22:39:39Z wojdyr $
 *
 *  Functions for reading and writing atomistic configurations in several
 *  file formats.
 */

#ifndef DEBYER_FILEIO_H_
#define DEBYER_FILEIO_H_

#include "debyer.h" /* for dbr_pbc and dbr_atom */


#ifdef __cplusplus
extern "C" {
#endif

/* ANSI C API */
int read_atoms_c(const char* filename, dbr_atom **atoms, dbr_pbc *pbc,
                 int reduced);
void free_atoms(dbr_atom *atoms);


/**************** the rest of this file is C++ only ******************/
#ifdef __cplusplus
}

#include <fstream>
#include <string>
#include <vector>
#include "lineio.h"

// atomistic configuration
struct dbr_aconf
{
    dbr_atom *atoms; // coordinates of atoms
    int n; // number of atoms

    // reduced coordinates (i.e. values in <0,1) range) are stored 
    // Some functions may not work well with reduced_coordinates.
    bool reduced_coordinates;

    dbr_pbc pbc; // PBC, zeroed if not applicable
    std::string orig_filename;
    std::vector<std::string> comments;  // no newlines here
};


bool is_xyz_format(char const* buffer);
bool is_plain_format(char const* /*buffer*/);

void read_xyz(LineInput& in, dbr_aconf *aconf);
void read_atomeye(LineInput& in, dbr_aconf *aconf, bool reduced_coords);
void read_dlpoly_config(LineInput& in, dbr_aconf *aconf);
void read_plain(LineInput& in, dbr_aconf *aconf);
void read_lammps_data(LineInput& in, dbr_aconf *aconf, bool reduced_coords);

/// write atomistic configuration aconf to file in XMOL XYZ format
void write_atoms_to_xyz_file(dbr_aconf const& aconf,
                             std::string const& filename);
/// write atomistic configuration aconf to file in AtomEye CFG format
void write_atoms_to_atomeye_file(dbr_aconf const& aconf,
                                 std::string const& filename);
/// write atomistic configuration aconf to file in DL_POLY format
/// if sorted is true, atoms will be sorted using symbol strings
void write_dlpoly_file(dbr_aconf const& aconf,
                       std::string const& filename, bool sorted);
/// write atomistic configuration aconf to file in LAMMPS input format
void write_lammps_data(dbr_aconf const& aconf, std::string const& filename);

/// write atomistic configuration aconf to file in PDB format
void write_pdb(dbr_aconf const& aconf, std::string const& filename);

/// write atomistic configuration aconf to file in plain TSV format:
/// x y z atom-symbol
void write_xyza(dbr_aconf const& aconf, std::string const& filename);

void write_file_with_atoms(dbr_aconf const& aconf, std::string const& filename);

// if reduced_coords is true, store reduced coordinates in <0,1) range 
// (relative to PBC box). Not all reading/writing functions respect it.
dbr_aconf read_atoms_from_file(LineInput &in, bool reduced_coords);

class SimpleLineInput : public LineInput
{
public:
    SimpleLineInput(const char* filename);
};

inline
dbr_aconf read_atoms_from_file(const char* filename, bool reduced_coords)
{
    SimpleLineInput lineinp(filename);
    return read_atoms_from_file(lineinp, reduced_coords);
}

inline
int open_atoms_file(const char* fn, dbr_atom **coords, dbr_pbc& pbc,
                    bool reduced_coords = true)
{
    return read_atoms_c(fn, coords, &pbc, reduced_coords);
}


// r = xyz * m
inline
void dbr_vec3_mult_pbc(dbr_xyz const xyz, dbr_pbc const& m, dbr_xyz r)
{
    //for (int j = 0; j < 3; ++j)
    //    r[j] = xyz[0] * m[0][j] + xyz[1] * m[1][j] + xyz[2] * m[2][j];
    r[0] = xyz[0] * m.v00 + xyz[1] * m.v10 + xyz[2] * m.v20;
    r[1] = xyz[0] * m.v01 + xyz[1] * m.v11 + xyz[2] * m.v21;
    r[2] = xyz[0] * m.v02 + xyz[1] * m.v12 + xyz[2] * m.v22;
}

inline
void dbr_vec3_mult_mat3x3(dbr_xyz const xyz, double const m[3][3], dbr_xyz r)
{
    for (int j = 0; j < 3; ++j)
        r[j] = xyz[0] * m[0][j] + xyz[1] * m[1][j] + xyz[2] * m[2][j];
}


#endif /* __cplusplus */
#endif /* DEBYER_FILEIO_H_ */

