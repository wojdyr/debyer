/*  debyer -- program for calculation of diffration patterns
 *  Copyright (C) 2006 Marcin Wojdyr
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
 */

#include "debyer.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>

#include "atomtables.h"
#include "iloops.h"


#ifndef VERSION
#   define VERSION "unknown"
#endif

#ifndef M_PI
# define M_PI    3.1415926535897932384626433832795029
#endif

int dbr_nid; /* rank of process (0 if serial) */
int dbr_noprocs; /* number of processes (1 if serial) */
time_t dbr_starttime;
int dbr_verbosity;

static void* xmalloc (size_t size)
{
    void *value = malloc(size);
    if (value == 0) {
        dbr_mesg("Error: virtual memory exhausted (%i bytes requested).\n",
                 size);
        dbr_abort(2);
    }
    return value;
}

static void* xrealloc (void *ptr, size_t size)
{
    void *value = realloc(ptr, size);
    if (value == 0) {
        dbr_mesg("Error: virtual memory exhausted  (%i bytes requested).\n",
                 size);
        dbr_abort(2);
    }
    return value;
}

/* printf-like function which outputs message to stderr (only from node 0) */
void dbr_mesg(char *fmt, ...)
{
    va_list args;

    if (dbr_nid != 0)
        return;
    fflush(stdout);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fflush(stderr);
}

/* get time (in seconds) from dbr_init() call */
int dbr_get_elapsed()
{
    time_t now;
    time(&now);
    return (int) difftime(now, dbr_starttime);
}

/* The same as Python % operator. In ANSI C a%b for a<0 has
 * implementation-defined sign, it is negative in popular implementations */
static int mod(int a, int b)
{
    int r = a % b;
    return r >= 0 ? r : r + b;
}

int dbr_is_direct(OutputKind k)
{
    return k == output_rdf || k == output_pdf || k == output_rpdf;
}

int dbr_is_inverse(OutputKind k)
{
    return k == output_xray || k == output_neutron || k == output_sf;
}

dbr_real get_xray_scattering_factor(const char* at, dbr_real q)
{
    dbr_real stol = q / (4*M_PI); /* Q -> sin(theta)/lambda */
    const t_it92_coeff *caa = find_in_it92(at);
    if (caa == 0) {
        dbr_mesg("Error: Scattering factor not found for atom: %s\n", at);
        dbr_abort(1);
    }
    return calculate_it92_factor(caa, stol*stol);
}

dbr_real get_neutron_scattering_factor(const char* at)
{
    const t_nn92_record *nbsl = find_in_nn92(at);
    if (nbsl == 0) {
        dbr_mesg("Error: Scattering factor not found for atom: %s\n", at);
        dbr_abort(1);
    }
    return nbsl->bond_coh_scatt_length;
}

void dbr_init(int *argc, char ***argv)
{
#ifdef USE_MPI
    int rc = MPI_Init(argc, argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &dbr_noprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &dbr_nid);
#else /* serial version */
    dbr_nid = 0;
    dbr_noprocs = 1;
    // suppress warning about unused vars
    assert(argv || argc);
#endif /*USE_MPI*/
    dbr_verbosity = 0;
    time(&dbr_starttime);
}

void dbr_finalize()
{
#ifdef USE_MPI
    MPI_Finalize();
#endif /*USE_MPI*/
}

void dbr_abort(int err_code)
{
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, err_code);
#else /* serial version */
    exit(err_code);
#endif /*USE_MPI*/
}


/* group atoms by name and store in `result',
 * return number of different names */
int dbr_get_atoms(int n, dbr_atom* coords, dbr_atoms** result,
                  int store_indices)
{
    int i, j;
    int type_count = 0;
    char *name = 0;
    dbr_atoms* atoms = *result;
    for (i = 0; i < n; ++i) {
        name = coords[i].name;
        for (j = 0; j < type_count; ++j) {
            dbr_atoms* a = &atoms[j];
            if (!strncmp(name, a->name, 8)) {
                if (a->count == a->asize) {
                    a->asize *= 2;
                    a->xyz = (dbr_xyz*) xrealloc(a->xyz,
                                                 a->asize * sizeof(dbr_xyz));
                    if (store_indices)
                        a->indices = (int*) xrealloc(a->indices,
                                                 a->asize * sizeof(int));
                }
                memcpy(a->xyz[a->count], coords[i].xyz, sizeof(dbr_xyz));
                if (store_indices)
                    a->indices[a->count] = i;
                a->count++;
                break;
            }
        }
        if (j == type_count) { /* new atom type */
            dbr_atoms *a = NULL;
            ++type_count;
            atoms = *result = (dbr_atoms*)
                              xrealloc(*result, type_count*sizeof(dbr_atoms));
            a = &atoms[type_count-1];
            strncpy(a->name, name, 8);
            a->asize = 8192;
            a->xyz = (dbr_xyz*) xmalloc(8192 * sizeof(dbr_xyz));
            a->indices = store_indices ? (int*) xmalloc(8192 * sizeof(int))
                                       : NULL;
            a->count = 1;
            memcpy(a->xyz[0], coords[i].xyz, sizeof(dbr_xyz));
            if (store_indices)
                a->indices[0] = i;
        }
    }
    /* shrink allocated space */
    i = 0;
    for (j = 0; j < type_count; ++j) {
        atoms[j].asize = atoms[j].count;
        atoms[j].xyz = (dbr_xyz*) xrealloc(atoms[j].xyz,
                                            atoms[j].asize * sizeof(dbr_xyz));
        i += atoms[j].asize * sizeof(dbr_xyz);
        if (store_indices) {
            atoms[j].indices = (int*) xrealloc(atoms[j].indices,
                                                atoms[j].asize * sizeof(int));
            i += atoms[j].asize * sizeof(int);
        }
    }
    if (dbr_verbosity > 1) /* very verbose */
        dbr_mesg("%i kbytes allocated for sorted atoms, %i kb for not sorted\n",
                i/1024, n*sizeof(dbr_atom)/1024);
    return type_count;
}

void free_dbr_atoms(dbr_atoms* xa)
{
    free(xa->xyz);
}

static
int get_number_of_bins(dbr_real rcut, dbr_real rquanta)
{
    return (int) ceil(rcut / rquanta);
}

/* open id file and write header */
static
FILE* write_header_of_id_file(irdfs rdfs, const char *filename)
{
    FILE *f;
    int i;
    if (dbr_nid != 0)
        return 0;
    assert(filename);
    f = fopen(filename, "wb");
    if (!f) {
        dbr_mesg("Error: can not open file for writing ID: %s\n", filename);
        return NULL;
    }
    fprintf(f, "debyer-id 1\n");
    fprintf(f, "# by debyer ver. " VERSION "\n");
    fprintf(f, "system");
    for (i = 0; i < rdfs.symbol_count; ++i)
        fprintf(f, " %s %i", rdfs.atom_symbols[i], rdfs.atom_counts[i]);
    fprintf(f, "\n");
    fprintf(f, "step %g\n", rdfs.step);
    fprintf(f, "bins %i\n", rdfs.rdf_bins);
    fflush(f);
    return f;
}

/* return 0 on success */
static
int append_pair_data_to_id_file(const irdf *p, int rdf_bins, FILE* f)
{
    int j, r;
    r = fprintf(f, "atoms %s %s", p->at1, p->at2);
    if (r < 0) {
        return r;
    }
    if (p->sample)
        fprintf(f, " sample %i", p->sample);
    fprintf(f, "\n");
    for (j = 0; j < rdf_bins; ++j)
        if (p->nn[j]) /* only non-zero bins are printed */
            fprintf(f, "%i %i\n", j, p->nn[j]);
    fflush(f);
    return 0;
}

/* O(n^2) TODO optimize - perhaps finding diagonal of containing */
dbr_real find_largest_distance(int n, const dbr_atoms* xa)
{
    int i, j, k, m;
    dbr_real max_sq_dist=0;
    for (i = 0; i < n; ++i) {
        for (k = 0; k < i; ++k) /*different atom types*/
            for (j = 0; j < xa[i].count; ++j) {
                for (m = 0; m < xa[k].count; ++m) {
                    if (get_sq_dist(xa[i].xyz[j], xa[k].xyz[m]) > max_sq_dist)
                        max_sq_dist = get_sq_dist(xa[i].xyz[j], xa[k].xyz[m]);
                }
            }
        for (j = 0; j < xa[i].count; ++j)  /*the same atom types*/
            for (m = 0; m < j; ++m) {
                if (get_sq_dist(xa[i].xyz[j], xa[i].xyz[m]) > max_sq_dist)
                    max_sq_dist = get_sq_dist(xa[i].xyz[j], xa[i].xyz[m]);
            }
    }
    return sqrt(max_sq_dist);
}

int dbr_is_atom_in_sector(const dbr_real *xyz, const dbr_picker* picker)
{
    return picker->x_min < xyz[0] && xyz[0] < picker->x_max
        && picker->y_min < xyz[1] && xyz[1] < picker->y_max
        && picker->z_min < xyz[2] && xyz[2] < picker->z_max;
}

static
char* pick_items(const dbr_picker* picker, const dbr_cell *c1)
{
    int i;
    int n = (int) ceil(c1->count * picker->probab / dbr_noprocs);
    int limit = c1->count / dbr_noprocs;
    char *t = (char*) xmalloc(c1->count);
    for (i = 0; i < c1->count; ++i)
        t[i] = (picker->probab == 0. ? 1. : 0.);
    if (picker->probab == 0. || limit == 0)
        ;
    else if (picker->probab < 0.5) { /*small sample*/
        for (i = 0; i < n; ) {
            int k = dbr_noprocs * (rand() % limit) + dbr_nid;
            if (t[k] == 0) {
                t[k] = 1;
                ++i;
            }
        }
    }
    else { /*large sample*/
        for (i = dbr_nid; i < c1->count; i += dbr_noprocs) {
            t[i] = (rand() < picker->probab * RAND_MAX);
        }
    }
    if (picker->cut) {
        for (i = dbr_nid; i < c1->count; i += dbr_noprocs)
            if (!dbr_is_atom_in_sector(c1->atoms[i], picker))
                t[i] = 0;
    }
    return t;
}

/* calculate irdf.nn using cell method
 * O(n), O(rcut^3) */
static
int calculate_irdf_cm(const dbr_picker* picker,
                      const dbr_cells cc1, const dbr_cells cc2,
                      dbr_real rcut, dbr_real rquanta, int nbins, int *t)
{
    int i, j, k;
    dbr_real rcut2 = rcut * rcut;
    int same = (cc1.data == cc2.data);
    int sampled = 0;
    assert(rcut >= 0.);

    for (i = 0; i < cc1.count; ++i) {
        char *picked = 0;
        const dbr_cell *c1 = &cc1.data[i];
        if (!c1->real)
            continue;
        if (!picker->all) {
            picked = pick_items(picker, c1);
            for (k = dbr_nid; k < c1->count; k += dbr_noprocs)
                if (picked[k])
                    sampled++;
        }
        for (j = 0; j < (same && picker->all ? 14 : 27); ++j) {
            const dbr_cell *c2;
            int nn = c1->neighbours[j];
            if (nn == -1)
                continue;
            c2 = &cc2.data[nn];

            if (picker->all) {
                /* Get complete data -- calculate irdf.nn using simple
                 * itaration (two loops) over all atoms.
                 */
                for (k = dbr_nid; k < c1->count; k += dbr_noprocs) {
                    if (rcut2 != 0.)
                        /* We check (r2 < rcut2) before calculating sqrt.
                         */
                        calculate_irdf_innerloop2(c1 == c2 ? k : c2->count,
                                                  nbins, rquanta,
                                                  c1->atoms[k], c2, t, rcut2);
                    else
                        /* We don't check (r2 < rcut2) before calculating sqrt.
                         * It's optimized for the case when most of
                         * the distances are <rcut
                         */
                        calculate_irdf_innerloop(c1 == c2 ? k : c2->count,
                                                 nbins, rquanta,
                                                 c1->atoms[k], c2, t);
                }
            }

            else {
                /* sampling */
                for (k = dbr_nid; k < c1->count; k += dbr_noprocs) {
                    if (!picked[k])
                        continue;
                    if (rcut2 != 0.)
                        calculate_irdf_innerloop2(c2->count, nbins, rquanta,
                                                  c1->atoms[k], c2, t, rcut2);
                    else
                        calculate_irdf_innerloop(c2->count, nbins, rquanta,
                                                 c1->atoms[k], c2, t);
                }
            }

        }
        free(picked);
    }
    if (same && !picker->all) // removing distance 0
        t[0] -= sampled;
    return sampled;
}


void dbr_inverse_3x3_matrix(const dbr_pbc a, double b[3][3])
{
    int i, j;
    b[0][0] =  (a.v11 * a.v22 - a.v12 * a.v21);
    b[0][1] = -(a.v01 * a.v22 - a.v02 * a.v21);
    b[0][2] =  (a.v01 * a.v12 - a.v02 * a.v11);
    b[1][0] = -(a.v10 * a.v22 - a.v12 * a.v20);
    b[1][1] =  (a.v00 * a.v22 - a.v02 * a.v20);
    b[1][2] = -(a.v00 * a.v12 - a.v02 * a.v10);
    b[2][0] =  (a.v10 * a.v21 - a.v11 * a.v20);
    b[2][1] = -(a.v00 * a.v21 - a.v01 * a.v20);
    b[2][2] =  (a.v00 * a.v11 - a.v01 * a.v10);

    dbr_real s = 1. / (a.v00 * b[0][0] + a.v01 * b[1][0] + a.v02 * b[2][0]);

    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            b[i][j] *= s;
}

/* works properly only with PBC in 3 dimensions */
dbr_pbc_prop get_pbc_properties(dbr_pbc pbc)
{
    dbr_pbc_prop p;
    double axb[3], bxc[3], cxa[3];

    p.vectors = pbc;
    p.lengths[0] = dbr_len3(pbc.v00, pbc.v01, pbc.v02);
    p.lengths[1] = dbr_len3(pbc.v10, pbc.v11, pbc.v12);
    p.lengths[2] = dbr_len3(pbc.v20, pbc.v21, pbc.v22);

    p.cosines[0] = (pbc.v00*pbc.v10 + pbc.v01*pbc.v11 + pbc.v02*pbc.v12)
                    / (p.lengths[0]*p.lengths[1]);
    p.cosines[1] = (pbc.v00*pbc.v20 + pbc.v01*pbc.v21 + pbc.v02*pbc.v22)
                    / (p.lengths[0]*p.lengths[2]);
    p.cosines[2] = (pbc.v10*pbc.v20 + pbc.v11*pbc.v21 + pbc.v12*pbc.v22)
                    / (p.lengths[1]*p.lengths[2]);

    /* vector products of cell vectors */
    axb[0]= pbc.v01 * pbc.v12 - pbc.v02 * pbc.v11;
    axb[1]= pbc.v02 * pbc.v10 - pbc.v00 * pbc.v12;
    axb[2]= pbc.v00 * pbc.v11 - pbc.v01 * pbc.v10;
    bxc[0]= pbc.v11 * pbc.v22 - pbc.v12 * pbc.v21;
    bxc[1]= pbc.v12 * pbc.v20 - pbc.v10 * pbc.v22;
    bxc[2]= pbc.v10 * pbc.v21 - pbc.v11 * pbc.v20;
    cxa[0]= pbc.v21 * pbc.v02 - pbc.v01 * pbc.v22;
    cxa[1]= pbc.v00 * pbc.v22 - pbc.v02 * pbc.v20;
    cxa[2]= pbc.v01 * pbc.v20 - pbc.v00 * pbc.v21;

    p.volume = fabs(pbc.v00*bxc[0] + pbc.v01*bxc[1] + pbc.v02*bxc[2]);

    p.widths[0] = p.volume / dbr_len3v(bxc);
    p.widths[1] = p.volume / dbr_len3v(cxa);
    p.widths[2] = p.volume / dbr_len3v(axb);

    return p;
}


/* it free()'s dbr_atoms xa[] (in prepare_cells()) */
irdfs calculate_irdfs(int n, dbr_atoms* xa, dbr_real rcut, dbr_real rquanta,
                      dbr_pbc pbc, const dbr_picker* picker,
                      const char* id_filename)
{
    int i, j, k, counter=0, all_atom_count=0;
    dbr_real max_r;
    irdfs rdfs;
    dbr_cells *cells = NULL;
    FILE *id = 0;
#ifdef USE_MPI
    int *tmp_nn = NULL;
#endif /*USE_MPI*/
    dbr_pbc_prop pbc_prop = get_pbc_properties(pbc);
    if (pbc.v00 != 0. || pbc.v11 != 0. || pbc.v22 != 0.) {
        if (rcut <= 0) {
            dbr_mesg("Error: cut-off must be specified for PBC system\n");
            dbr_abort(1);
        }
        if (rcut >= pbc_prop.widths[0]/2 || rcut >= pbc_prop.widths[1]/2
                                         || rcut >= pbc_prop.widths[2]/2) {
            dbr_mesg("Error: cut-off must be smaller than half of PBC box\n");
            dbr_abort(1);
        }
    }
    max_r = rcut > 0 ? rcut : find_largest_distance(n, xa);
    rdfs.step = rquanta;
    rdfs.rdf_bins = get_number_of_bins(max_r, rquanta);
    rdfs.symbol_count = n;
    rdfs.atom_symbols = (dbr_symbol*) xmalloc(n * sizeof(dbr_symbol));
    rdfs.atom_counts = (int*) xmalloc(n * sizeof(int));
    for (i = 0; i < n; ++i) {
        strncpy(rdfs.atom_symbols[i], xa[i].name, 8);
        rdfs.atom_counts[i] = xa[i].count;
        all_atom_count += xa[i].count;
    }
    rdfs.pair_count = n * (n+1) / 2;
    rdfs.data = (irdf*) xmalloc(rdfs.pair_count * sizeof(irdf));
    rdfs.density = -1.;
    if (pbc.v00 != 0. && pbc.v11 != 0. && pbc.v22 != 0.)
        rdfs.density = all_atom_count / pbc_prop.volume;
    if (id_filename) {
        if (dbr_verbosity > 0)
            dbr_mesg("ID will be written to file: %s.\n", id_filename);
        id = write_header_of_id_file(rdfs, id_filename);
        //if (id && dbr_nid == 0)
        //    write_comments(id, );
        if (id && dbr_nid == 0 && rdfs.density > 0)
            fprintf(id, "# numeric-density %g\n", rdfs.density);
    }

    cells = prepare_cells_all(pbc, max_r, xa, n);

#ifdef USE_MPI
    tmp_nn = (int*) xmalloc(rdfs.rdf_bins * sizeof(int));
#endif /*USE_MPI*/
    if (picker->probab > 0) {
        srand(time(NULL)+133*dbr_nid);
        if (dbr_verbosity > 0)
            dbr_mesg("Sampling about %g%% of atoms.\n", picker->probab*100);
    }
    for (i = 0; i < n; ++i) { /* first atom type */
        for (j = 0; j <= i; ++j) { /* second atom type */
            irdf *p = &rdfs.data[counter++];
            strncpy(p->at1, cells[i].name, 8);
            p->c1 = cells[i].atom_count;
            strncpy(p->at2, cells[j].name, 8);
            p->c2 = cells[j].atom_count;
            p->nn = (int*) xmalloc(rdfs.rdf_bins * sizeof(int));
            for (k = 0; k < rdfs.rdf_bins; ++k)
                p->nn[k] = 0;
            p->sample = calculate_irdf_cm(picker, cells[i], cells[j],
                                          rcut, rquanta, rdfs.rdf_bins, p->nn);

            /* if we are only sampling, we also need to swap i and j */
            if (i != j && !picker->all) {
                p->sample += calculate_irdf_cm(picker, cells[j], cells[i],
                                          rcut, rquanta, rdfs.rdf_bins, p->nn);
            }

#ifdef USE_MPI
            MPI_Reduce(p->nn, tmp_nn, rdfs.rdf_bins, MPI_INT, MPI_SUM,
                       0, MPI_COMM_WORLD);
            MPI_Reduce(&p->sample, &k, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            if (dbr_nid == 0)
                memcpy(p->nn, tmp_nn, rdfs.rdf_bins * sizeof(int));
            p->sample = k;
#endif /*USE_MPI*/
            if (!picker->all) {
                if (i != j)
                    p->sample /= 2;
                for (k = 0; k < rdfs.rdf_bins; ++k) {
                    // rounding errors may change the bin
                    // if p->nn[k] is odd, round either up or down
                    if (p->nn[k] % 2 != 0 && k % 2)
                        p->nn[k]++;
                    p->nn[k] /= 2;
                }
            }
            if (id) { /*executed only for dbr_nid == 0 */
                int r = append_pair_data_to_id_file(p, rdfs.rdf_bins, id);
                if (r != 0)
                    dbr_mesg("Error: can not write to file: %s\n", id_filename);
            }
            if (dbr_verbosity > 0) {
                dbr_mesg("Elapsed %i s. %s-%s ID calculated",
                         dbr_get_elapsed(), p->at1, p->at2);
                if (picker->probab)
                    dbr_mesg(" (sampled %i %s atoms)", p->sample, p->at1);
                dbr_mesg(".\n");
            }
        }
    }
#ifdef USE_MPI
    free(tmp_nn);
#endif /*USE_MPI*/

    free_cells_all(cells, n);

    if (id)
        fclose(id);
    return rdfs;
}

dbr_cells* prepare_cells_all(dbr_pbc pbc, dbr_real rcut, dbr_atoms* xa, int n)
{
    int i;
    dbr_cells *cells = NULL;
    if (dbr_verbosity > 0)
        dbr_mesg("Preparing cells for cell method...\n");
    cells = (dbr_cells*) xmalloc(n * sizeof(dbr_cells));
    for (i = 0; i < n; ++i) {
        cells[i] = prepare_cells(pbc, rcut, xa + i);
    }
    free(xa);
    if (dbr_verbosity > 0) {
        if (cells[0].count > 1)
            dbr_mesg("... %i x %i x %i cells.\n",
                     cells[0].n[0], cells[0].n[1], cells[0].n[2]);
        else
            dbr_mesg("... space _not_ divided into cells\n");
    }
    if (dbr_verbosity > 0)
        dbr_mesg("Elapsed %i s. Cells are ready.\n", dbr_get_elapsed());
    return cells;
}



void free_cells_all(dbr_cells *cells, int n)
{
    int i;
    if (cells) {
        for (i = 0; i < n; ++i)
            free_cells(cells[i]);
        free(cells);
    }
}

void write_irdfs_to_file(irdfs rdfs, const char *filename)
{
    int i, counter=0;
    FILE *id = write_header_of_id_file(rdfs, filename);
    for (i = 0; i < rdfs.pair_count; ++i) { /* first atom type */
        irdf *p = &rdfs.data[counter];
        append_pair_data_to_id_file(p, rdfs.rdf_bins, id);
        counter++;
    }
    fclose(id);
}

static
void skip_to_next_line(FILE *f)
{
    int r;
    while (1) {
        while ((r = fgetc(f)) != '\n' && r != EOF)
            ;
        if (r == EOF)
            break;
        r = fgetc(f);
        if (r == EOF)
            break;
        if (r != '#') {
            ungetc(r, f);
            break;
        }
    }
}

/* return next character */
static
int skip_blanks(FILE *f)
{
    while(1) {
        int r = fgetc(f);
        if (r == EOF)
            return r;
        else if (!isspace(r) || r == '\n') {
            ungetc(r, f);
            return r;
        }
    }
}

irdfs read_irdfs_from_file(const char *filename)
{
    FILE *f;
    irdfs rdfs;
    irdf *p = 0;
    int i, j, k, r, sn, version, count, next;
    dbr_symbol name;
    rdfs.data = 0;
    rdfs.rdf_bins = 0;
    rdfs.pair_count = 0;
    f = fopen(filename, "r");
    if (!f) {
        dbr_mesg("Error: can not open file: %s\n", filename);
        return rdfs;
    }
    /* read and verify first line */
    r = fscanf(f, "debyer-id %i", &version);
    if (r != 1) {
        dbr_mesg("Error: debyer id header not found in file: %s\n", filename);
        return rdfs;
    }
    else if (version != 1) {
        dbr_mesg("Error: debyer id version %i is not supported by debyer "
                        VERSION "\n", version);
        return rdfs;
    }
    skip_to_next_line(f);
    /* read system */
    sn = fscanf(f, " system %7s %i", name, &count);
    if (sn != 2) {
        dbr_mesg("Error: `system' line not found in id file: %s\n",
                        filename);
        return rdfs;
    }
    rdfs.atom_symbols = (dbr_symbol*) xmalloc(1 * sizeof(dbr_symbol));
    rdfs.atom_counts = (int*) xmalloc(1 * sizeof(int));
    strncpy(rdfs.atom_symbols[0], name, 8);
    rdfs.atom_counts[0] = count;
    i = 1;
    next = skip_blanks(f);
    while (next != '\n' && fscanf(f, "%7s %i", name, &count) == 2) {
        ++i;
        rdfs.atom_symbols = (dbr_symbol*) xrealloc(rdfs.atom_symbols,
                                                   i * sizeof(dbr_symbol));
        rdfs.atom_counts = (int*) xrealloc(rdfs.atom_counts, i * sizeof(int));
        strncpy(rdfs.atom_symbols[i-1], name, 8);
        rdfs.atom_counts[i-1] = count;
        next = skip_blanks(f);
    }
    rdfs.symbol_count = i;
    /* read step */
    sn = fscanf(f, " step " DBR_F, &rdfs.step);
    if (sn != 1) {
        dbr_mesg("Error: `step' line not found in id file: %s\n",
                        filename);
        return rdfs;
    }
    skip_to_next_line(f);
    /* read bins */
    sn = fscanf(f, " bins %i", &rdfs.rdf_bins);
    if (sn != 1) {
        dbr_mesg("Error: `bins' line not found in id file: %s\n",
                        filename);
        return rdfs;
    }
    /* read optional numeric density */
    sn = fscanf(f, " # numeric-density "DBR_F, &rdfs.density);
    if (sn == 1 && dbr_verbosity > 1) /* very verbose */
        dbr_mesg("Numeric density read from id file: %g\n", rdfs.density);

    fscanf(f, " ");
    /* read the rest of the file */
    while (1) {
        r = fgetc(f);
        if (r == 'a') {
            ++rdfs.pair_count;
            rdfs.data =
                (irdf*) xrealloc(rdfs.data, rdfs.pair_count*sizeof(irdf));
            p = &rdfs.data[rdfs.pair_count-1];
            p->nn = (int*) xmalloc(rdfs.rdf_bins * sizeof(int));
            for (i = 0; i < rdfs.rdf_bins; ++i)
                p->nn[i] = 0;
            sn = fscanf(f, "toms %7s %7s", p->at1, p->at2);
            if (sn != 2) {
                dbr_mesg("Error: unexpected format of file "
                         "(when trying to read atoms line): %s\n", filename);
                free_irdfs(&rdfs);
                return rdfs;
            }
            p->sample = 0;
            for (i = 0; i < rdfs.symbol_count; ++i) {
                if (!strncmp(p->at1, rdfs.atom_symbols[i], 8))
                    p->c1 = rdfs.atom_counts[i];
                if (!strncmp(p->at2, rdfs.atom_symbols[i], 8))
                    p->c2 = rdfs.atom_counts[i];
            }
            sn = fscanf(f, " sample %i", &p->sample); /* optional "sample" */
        }
        else if (r == '#') {
            skip_to_next_line(f);
        }
        else if (r == EOF) {
            break;
        }
        else {
            if (p == 0) {
                dbr_mesg("Error: unexpected format of id file "
                         "(\"atoms\" line expected): %s\n", filename);
                return rdfs;
            }
            ungetc(r, f);
            sn = fscanf(f, "%i %i ", &j, &k);
            if (sn != 2) {
                free_irdfs(&rdfs);
                dbr_mesg("Error: unexpected format of id file: %s\n", filename);
                return rdfs;
            }
            if (j >= rdfs.rdf_bins) {
                free_irdfs(&rdfs);
                dbr_mesg("Error: unexpected bin %i in id file: %s\n",
                         j, filename);
                return rdfs;
            }
            p->nn[j] = k;
        }
    }

    fclose(f);
    return rdfs;
}


void free_irdfs(irdfs *rdfs)
{
    int i;
    free(rdfs->atom_symbols);
    free(rdfs->atom_counts);
    if (rdfs->data == NULL) {
        rdfs->pair_count = 0;
        return;
    }
    for (i = 0; i < rdfs->pair_count; ++i)
        free(rdfs->data[i].nn);
    free(rdfs->data);
    rdfs->data = NULL;
    rdfs->pair_count = 0;
}

static
dbr_real calculate_avg_b(char weight, const irdfs* rdfs, dbr_real q)
{
    int i, c = 0;
    dbr_real sum = 0.;
    if (weight != 'x' && weight != 'n')
        return 1.;
    for (i = 0; i < rdfs->symbol_count; ++i) {
        int count = rdfs->atom_counts[i];
        const char* symbol = rdfs->atom_symbols[i];
        c += count;
        if (weight == 'x')
            sum += count * get_xray_scattering_factor(symbol, q);
        else if (weight == 'n')
            sum += count * get_neutron_scattering_factor(symbol);
    }
    return sum / c;
}

static
dbr_real count_all_atoms(const irdfs* rdfs)
{
    int i, count = 0;
    for (i = 0; i < rdfs->symbol_count; ++i)
        count += rdfs->atom_counts[i];
    return count;
}

dbr_real* dbr_get_RDF(const irdfs* rdfs, int rdf_index, struct dbr_pdf_args* p)
{
    dbr_real *pattern;
    int i, j, m, n;
    int first_ibin, end_ibin;
    dbr_real avg = calculate_avg_b(p->weight, rdfs, 0);
    int all_count = count_all_atoms(rdfs);

    n = (p->pattern_to - p->pattern_from) / p->pattern_step;
    pattern = (dbr_real*) xmalloc(n*sizeof(dbr_real));
    for (i = 0; i < n; ++i)
        pattern[i] = 0;
    first_ibin = p->pattern_from / rdfs->step;
    end_ibin = p->pattern_to / rdfs->step;
    if (end_ibin > rdfs->rdf_bins)
        end_ibin = rdfs->rdf_bins;

    for (i = 0; i < rdfs->pair_count; ++i) {
        if (rdf_index != -1 && rdf_index != i)
            continue;
        /* prepare weighting factor... */
        dbr_real w = 1.;
        irdf *rdf = &rdfs->data[i];
        dbr_real sampled_fraction = (rdf->sample > 0 ?
                                         (dbr_real) rdf->sample / rdf->c1 : 1.);
        if (p->weight == 'x')
            w = get_xray_scattering_factor(rdf->at1, 0)
                  * get_xray_scattering_factor(rdf->at2, 0) / (avg * avg);
        else if (p->weight == 'n')
            w = get_neutron_scattering_factor(rdf->at1)
                * get_neutron_scattering_factor(rdf->at2) / (avg * avg);
        /* ...and calculate RDF */
        for (j = first_ibin; j < end_ibin; ++j) {
            m = (int) (((j + 0.5) * rdfs->step - p->pattern_from)
                       / p->pattern_step);
            assert(m >= 0 && m < n);
            /* take twice every pair distance (i,j and j,i) */
            pattern[m] += 2 * w * rdf->nn[j] / sampled_fraction;
        }
    }
    for (i = 0; i < n; ++i)
        pattern[i] /= (p->pattern_step * all_count);
    return pattern;
}

dbr_real get_density(dbr_real given_density, dbr_real auto_density)
{
    if (dbr_verbosity > 0 && given_density >= 0. &&
            auto_density > 0 && given_density != auto_density)
        dbr_mesg("Ignoring number density from ID (%g), using %g\n",
                 auto_density, given_density);

    return given_density >= 0. ? given_density : auto_density;
}

static
void add_cutoff_correction(const irdfs* rdfs,
                           const struct dbr_diffract_args* dargs,
                           int n_pattern, dbr_real *pattern)
{
    int j;
    char weight = '1'; // '1' is for output_sf

    if (dargs->c == output_xray)
        weight = 'x';
    else if (dargs->c == output_neutron)
        weight = 'n';
    for (j = 0; j < n_pattern; ++j) {
        dbr_real x = dargs->pattern_from + (j+0.5) * dargs->pattern_step;
        dbr_real q = dargs->lambda <= 0. ? x
                            : 4*M_PI * sin(M_PI/180.*x/2) / dargs->lambda;
        dbr_real avg = calculate_avg_b(weight, rdfs, q);
        dbr_real qc = q * dargs->cutoff;
        pattern[j] += avg * avg * 4 * M_PI * dargs->ro / (q * q) * (
                                   dargs->cutoff * cos(qc) - sin(qc) / q);
    }
}

dbr_real* get_pattern(const irdfs* rdfs, struct dbr_diffract_args* dargs)
{
    dbr_real *pattern;
    int i, j, k, n;
    irdf *p = NULL;
    dbr_real ff = 1.;
    int all_count = count_all_atoms(rdfs);

    assert(dargs->c == output_xray || dargs->c == output_neutron ||
           dargs->c == output_sf);
    for (i = 0; i < rdfs->symbol_count; ++i) {
        const char* at = rdfs->atom_symbols[i];
        if (dargs->c == output_xray && find_in_it92(at) == NULL) {
            dbr_mesg("Error: No x-ray scattering factor for atom: %s\n", at);
            return NULL;
        }
        else if (dargs->c == output_neutron && find_in_nn92(at) == NULL) {
            dbr_mesg("Error: No neutron scattering factor for atom: %s\n", at);
            return NULL;
        }
    }
    if (dargs->cutoff > 0)
        /* use calculated density unless different value is explicitly given */
        dargs->ro = get_density(dargs->ro, rdfs->density);

    n = (dargs->pattern_to - dargs->pattern_from) / dargs->pattern_step;
    pattern = (dbr_real*) xmalloc(n*sizeof(dbr_real));
    for (i = 0; i < n; ++i)
        pattern[i] = 0.;
    for (i = 0; i < rdfs->pair_count; ++i) {
        p = &rdfs->data[i];
        if (dargs->c == output_neutron) {
            ff = get_neutron_scattering_factor(p->at1)
                                * get_neutron_scattering_factor(p->at2);
        }
        for (j = 0; j < n; ++j) {
            dbr_real x = dargs->pattern_from + (j+0.5) * dargs->pattern_step;
            dbr_real q = dargs->lambda <= 0. ? x
                                 : 4*M_PI * sin(M_PI/180.*x/2) / dargs->lambda;
            if (dargs->c == output_xray)
                ff = get_xray_scattering_factor(p->at1, q)
                                    * get_xray_scattering_factor(p->at2, q);
            if (!dargs->sinc_damp) {
                for (k = 0; k < rdfs->rdf_bins; ++k) {
                    if (p->nn[k]) {
                        dbr_real r = (k+0.5) * rdfs->step;
                        /* doubled, because in Debye's formula every pair
                         * is taken twice (n,m and m,n) */
                        pattern[j] += ff * (2 * p->nn[k]) * sin(q*r) / (q*r);
                    }
                }
            } else {
                // shell_volume = 4 pi r^2 * rdfs->step
                // n_cont = shell_volume * ro * c1 * c2 / all_count
                dbr_real t = 4*M_PI * rdfs->step * // r^2 is below
                             dargs->ro * p->c1 * p->c2 / all_count;
                for (k = 0; k < rdfs->rdf_bins; ++k) {
                    dbr_real r = (k+0.5) * rdfs->step;
                    dbr_real r_c = dargs->cutoff;
                    dbr_real damp_factor = sin(M_PI*r/r_c) / (M_PI*r/r_c);
                    dbr_real nn = (2 * p->nn[k] - t*r*r) * damp_factor;
                    pattern[j] += ff * nn * sin(q*r) / (q*r);
                }
            }
            /* in Debye's formula, in \sum_{n,m}, also n==m counts */
            if (!strcmp(p->at1, p->at2)) {
                pattern[j] += ff * p->c1;
            }
        }
    }
    for (i = 0; i < n; ++i)
        pattern[i] /= all_count;

    /* this correction has been already applied if sinc_damp is given */
    if (dargs->cutoff > 0 && !dargs->sinc_damp) {
        dargs->ro = get_density(dargs->ro, rdfs->density);
        if (dbr_verbosity > 0 && dargs->ro >= 0)
            dbr_mesg("Numeric density: %g\n", dargs->ro);
        if (dargs->ro > 0.)
            add_cutoff_correction(rdfs, dargs, n, pattern);
    }

    return pattern;
}

static
FILE* start_pattern_output(const char *ofname)
{
    FILE *f;
    if (!ofname || !*ofname || !strcmp(ofname, "-"))
        f = stdout;
    else {
        f = fopen(ofname, "w");
        if (!f) {
            dbr_mesg("Error: can not open file: %s\n", ofname);
            dbr_abort(1);
        }
    }
    fprintf(f, "#debyer-pattern ");
    return f;
}

static
int process_output_step(dbr_real pattern_from, dbr_real pattern_to,
                        dbr_real pattern_step, FILE* f)
{
    int n = (pattern_to - pattern_from) / pattern_step;
    fprintf(f, "\n#from %g to %g step %g\n",
            pattern_from, pattern_to, pattern_step);
    if (n < 1) {
        dbr_mesg("Error: Incorrect values for pattern from/to/step.\n");
        dbr_abort(1);
    }
    return n;
}

/* returns non-zero on failure */
int write_diffraction_to_file(struct dbr_diffract_args* dargs, irdfs rdfs,
                              const char *ofname)
{
    FILE *f;
    int i, n;
    dbr_real *result=NULL;
    assert(dbr_is_inverse(dargs->c));
    if (dbr_nid != 0)
        return 0;
    f = start_pattern_output(ofname);
    if (dargs->c == output_xray)
        fprintf(f, "x-ray");
    else if (dargs->c == output_neutron)
        fprintf(f, "neutron");
    else if (dargs->c == output_sf)
        fprintf(f, "scattering-function");
    if (dargs->lambda > 0.)
        fprintf(f, " lambda=%g", dargs->lambda);
    else
        fprintf(f, " Q");
    if (dargs->cutoff > 0.) {
        int nb = get_number_of_bins(dargs->cutoff, rdfs.step);
        if (nb < rdfs.rdf_bins) {
            rdfs.rdf_bins = nb;
        }
        else if (nb > rdfs.rdf_bins) {
            dargs->cutoff = rdfs.rdf_bins * rdfs.step;
            dbr_mesg("WARNING: can't set cut-off larger than %g\n",
                     dargs->cutoff);
        }
    }
    fprintf(f, " cut-off=%g", dargs->cutoff > 0 ? dargs->cutoff
                                                : rdfs.rdf_bins * rdfs.step);
    if (dargs->pattern_from <= 0.)
        dargs->pattern_from = dargs->lambda <= 0. ? 0.5 : 5;
    if (dargs->pattern_to <= 0.)
        dargs->pattern_to = dargs->lambda <= 0. ? 10 : 170;
    if (dargs->pattern_step <= 0.) {/*TODO Nyquist */
        dargs->pattern_step = dargs->lambda <= 0. ? 0.01 : 0.1;
    }

    n = process_output_step(dargs->pattern_from, dargs->pattern_to,
                            dargs->pattern_step, f);

    result = get_pattern(&rdfs, dargs);

    for (i = 0; i < n; ++i)
        fprintf(f, "%g %g\n", dargs->pattern_from + (i+0.5)*dargs->pattern_step,
                              result[i]);
    if (f != stdout)
        fclose(f);
    free(result);
    return 0;
}

dbr_real scale_rdf(dbr_real y, OutputKind c, dbr_real r, dbr_real ro)
{
    if (c == output_pdf)
        return y / (4 * M_PI * r * r * ro);
    else if (c == output_rpdf)
        return y / r - (4 * M_PI * r * ro);
    else // output_rdf
        return y;
}

/* returns non-zero on failure */
int write_pdfkind_to_file(struct dbr_pdf_args* pargs, irdfs rdfs,
                          const char *ofname)
{
    FILE *f;
    int i, j, n, pc;
    dbr_real *result = NULL;
    dbr_real **partial_results = NULL;
    dbr_real irdf_max = rdfs.rdf_bins * rdfs.step;
    assert(dbr_is_direct(pargs->c));
    if (dbr_nid != 0)
        return 0;
    f = start_pattern_output(ofname);

    if (pargs->c == output_rdf)
        fprintf(f, "RDF");
    else if (pargs->c == output_pdf)
        fprintf(f, "PDF");
    else if (pargs->c == output_rpdf)
        fprintf(f, "reduced-PDF");

    if (pargs->c != output_rdf) {
        pargs->ro = get_density(pargs->ro, rdfs.density);
        if (pargs->ro <= 0.) {
            dbr_mesg("Error: Unknown value of density (ro).\n");
            dbr_abort(1);
        }
        fprintf(f, " ro=%g", pargs->ro);
    }

    fprintf(f, " weight=%c", pargs->weight);

    if (pargs->pattern_from <= 0.)
        pargs->pattern_from = 0.;
    if (pargs->pattern_to <= 0.)
        pargs->pattern_to = irdf_max;
    else if (pargs->pattern_to > irdf_max) {
        dbr_mesg("WARNING: ID is calculated only to: %g\n", irdf_max);
        pargs->pattern_to = irdf_max;
    }
    if (pargs->pattern_step <= 0.)
        pargs->pattern_step = rdfs.step;

    n = process_output_step(pargs->pattern_from, pargs->pattern_to,
                            pargs->pattern_step, f);

    result = dbr_get_RDF(&rdfs, -1, pargs);
    if (pargs->include_partials) {
        pc = rdfs.pair_count;
        partial_results = (dbr_real**) xmalloc(pc * sizeof(dbr_real*));
        fprintf(f, "# sum");
        for (i = 0; i != pc; ++i) {
            partial_results[i] = dbr_get_RDF(&rdfs, i, pargs);
            fprintf(f, " %s-%s", rdfs.data[i].at1, rdfs.data[i].at2);
        }
        fprintf(f, "\n");
    }

    for (i = 0; i < n; ++i) {
        dbr_real r = pargs->pattern_from + (i+0.5) * pargs->pattern_step;
        dbr_real y = scale_rdf(result[i], pargs->c, r, pargs->ro);
        fprintf(f, "%g %g", r, y);
        if (pargs->include_partials) {
            for (j = 0; j != pc; ++j) {
                y = scale_rdf(partial_results[j][i], pargs->c, r, pargs->ro);
                fprintf(f, " %g", y);
            }
        }
        fprintf(f, "\n");
    }

    if (f != stdout)
        fclose(f);
    free(result);
    if (pargs->include_partials) {
        for (i = 0; i != pc; ++i)
            free(partial_results[i]);
        free(partial_results);
    }
    return 0;
}

/* virtual cells are (if any) with a = 0, a = cells->n[0], etc */
static
int get_cell_nr(int a, int b, int c, dbr_cells *cells)
{
    a += cells->v[0];
    b += cells->v[1];
    c += cells->v[2];
    return (a * cells->a[1] + b) * cells->a[2] + c;
}

static
dbr_cell* get_cell(int a, int b, int c, dbr_cells *cells)
{
    return &cells->data[get_cell_nr(a, b, c, cells)];
}

/* a,b,c point to absolute cell position, ie. a=0 can give virtual cell */
static
dbr_cell* get_abs_cell(int a, int b, int c, dbr_cells *cells)
{
    int pos = (a * cells->a[1] + b) * cells->a[2] + c;
    return &cells->data[pos];
}

static
void make_virtual_cell(dbr_cell *dest, dbr_cell *src,
                       dbr_real dx, dbr_real dy, dbr_real dz)
{
    int i;
    assert (src->real);
    dest->asize = src->asize;
    dest->count = src->count;
    dest->real = 0;
    dest->original = src->neighbours[13];
    dest->atoms = (dbr_xyz*) xmalloc(dest->asize * sizeof(dbr_xyz));
    for (i = 0; i < src->count; ++i) {
        dest->atoms[i][0] = src->atoms[i][0] + dx;
        dest->atoms[i][1] = src->atoms[i][1] + dy;
        dest->atoms[i][2] = src->atoms[i][2] + dz;
    }
    for (i = 0; i < 27; ++i)
        dest->neighbours[i] = -1;
}

static
void set_neighbours(int a, int b, int c, dbr_cells *cells)
{
    dbr_cell *cell = get_cell(a, b, c, cells);
    int i, j, k, counter = 0;
    for (i = a-1; i <= a+1; ++i)
        for (j = b-1; j <= b+1; ++j)
            for (k = c-1; k <= c+1; ++k) {
                if ((!cells->v[0] && (i == -1 || i == cells->n[0]))
                        || (!cells->v[1] && (j == -1 || j == cells->n[1]))
                        || (!cells->v[2] && (k == -1 || k == cells->n[2])))
                    cell->neighbours[counter] = -1;
                else
                    cell->neighbours[counter] = get_cell_nr(i, j, k, cells);
                counter++;
            }
}

void print_cells_memory(dbr_cells cells)
{
    int i, j, k;
    int rc_at=0, vc_at=0, count_empty=0;
    for (i = 0; i < cells.n[0]; ++i) {
        for (j = 0; j < cells.n[1]; ++j) {
            for (k = 0; k < cells.n[2]; ++k) {
                const dbr_cell* cell = get_cell(i, j, k, &cells);
                if (cell->count == 0) {
                    if (count_empty == 0)
                        dbr_mesg("Empty cells for %s:", cells.name);
                    dbr_mesg(" (%i,%i,%i)", i, j, k);
                    count_empty++;
                }
                rc_at += cell->asize;
                vc_at += ((int)(i == 0) + (int)(i == cells.n[0] - 1)
                        + (int)(j == 0) + (int)(j == cells.n[1] - 1)
                        + (int)(k == 0) + (int)(k == cells.n[2] - 1))
                                                            * cell->asize;
            }
        }
    }
    if (count_empty > 0)
        dbr_mesg("\n");
    dbr_mesg("Memory for atoms in %s cells: %i + %i = %i kb\n",
             cells.name,
             rc_at * sizeof(dbr_xyz) / 1024,
             vc_at * sizeof(dbr_xyz) / 1024,
             (rc_at + vc_at) * sizeof(dbr_xyz) / 1024);
}

/* it free()'s dbr_atoms xa[] */
dbr_cells prepare_cells(dbr_pbc pbc, dbr_real rcut, dbr_atoms* xa)
{
    int ini_size, i, j, k, m, n[3];
    dbr_cells cells;
    double inv_v[3][3];
    dbr_pbc_prop pbc_prop = get_pbc_properties(pbc);
    rcut += 1e-3;

    /* virtual cells are used in PBC systems, as images of the first and last
     * cells, only in periodic directions */
    cells.v[0] = (pbc.v00 == 0. ? 0 : 1);
    cells.v[1] = (pbc.v11 == 0. ? 0 : 1);
    cells.v[2] = (pbc.v22 == 0. ? 0 : 1);

    /* calculate number of cells */
    for (i = 0; i < 3; ++i) {
        if (cells.v[i] == 0)  /* no pbc, put all in one cell */
            cells.n[i] = 1;
        else { /* cell dimension must be smaller than rcut */
            int t = (int) floor(pbc_prop.widths[i] / rcut);
            cells.n[i] = t > 0 ? t : 1;
        }
    }
    cells.vectors.v00 = pbc.v00 / cells.n[0];
    cells.vectors.v01 = pbc.v01 / cells.n[0];
    cells.vectors.v02 = pbc.v02 / cells.n[0];
    cells.vectors.v10 = pbc.v10 / cells.n[1];
    cells.vectors.v11 = pbc.v11 / cells.n[1];
    cells.vectors.v12 = pbc.v12 / cells.n[1];
    cells.vectors.v20 = pbc.v20 / cells.n[2];
    cells.vectors.v21 = pbc.v21 / cells.n[2];
    cells.vectors.v22 = pbc.v22 / cells.n[2];

    /* calculate the number of all cells (real + virtual) */
    for (i = 0; i < 3; ++i)
        cells.a[i] = cells.n[i] + 2 * cells.v[i];
    cells.count = cells.a[0] * cells.a[1] * cells.a[2];

    /* initialize cells */
    cells.data = (dbr_cell*) xmalloc(cells.count * sizeof(dbr_cell));
    for (i = 0; i < cells.count; ++i) {
        strncpy(cells.name, xa->name, 8);
        cells.atom_count = xa->count;
        for (j = 0; j < 27; ++j)
            cells.data[i].neighbours[j] = -1;
        cells.data[i].real = 0;
    }

    /* special case, single cell */
    if (cells.count == 1) {
        cells.data[0].neighbours[13] = 0;
        cells.data[0].asize = xa->asize;
        cells.data[0].count = xa->count;
        cells.data[0].atoms = xa->xyz;
        cells.data[0].indices = xa->indices;
        cells.data[0].real = 1;
        cells.data[0].original = -1;
        xa->xyz = NULL;
        xa->indices = NULL;
        return cells;
    }

    /* real cells initializations */
    ini_size = xa->count / (cells.n[0] * cells.n[1] * cells.n[2]);
    if (ini_size < 4)
        ini_size = 4;
    for (i = 0; i < cells.n[0]; ++i) {
        for (j = 0; j < cells.n[1]; ++j) {
            for (k = 0; k < cells.n[2]; ++k) {
                dbr_cell *cell = get_cell(i, j, k, &cells);
                cell->asize = ini_size;
                cell->count = 0;
                cell->atoms = (dbr_xyz*) xmalloc(ini_size * sizeof(dbr_xyz));
                if (xa->indices)
                    cell->indices = (int*) xmalloc(ini_size * sizeof(int));
                cell->real = 1;
                cell->original = -1;
            }
        }
    }
    /* set neighbour numbers */
    for (i = 0; i < cells.n[0]; ++i) {
        for (j = 0; j < cells.n[1]; ++j) {
            for (k = 0; k < cells.n[2]; ++k) {
                set_neighbours(i, j, k, &cells);
            }
        }
    }

    /* put atoms into real cells */
    for (i = 0; i < xa->count; ++i) {
        dbr_cell *cell;
        double vec[3][3] = {
            { cells.vectors.v00, cells.vectors.v01, cells.vectors.v02 },
            { cells.vectors.v10, cells.vectors.v11, cells.vectors.v12 },
            { cells.vectors.v20, cells.vectors.v21, cells.vectors.v22 } };
        dbr_inverse_3x3_matrix(cells.vectors, inv_v);
        for (j = 0; j < 3; ++j) {
            n[j] = (int) floor(  xa->xyz[i][0] * inv_v[0][j]
                               + xa->xyz[i][1] * inv_v[1][j]
                               + xa->xyz[i][2] * inv_v[2][j] );
            if (n[j] < 0 || n[j] >= cells.n[j]) {
                int t = mod(n[j], cells.n[j]);
                int diff = t - n[j];
                xa->xyz[i][0] += diff * vec[0][j];
                xa->xyz[i][1] += diff * vec[1][j];
                xa->xyz[i][2] += diff * vec[2][j];
                n[j] = t;
            }
        }
        cell = get_cell(n[0], n[1], n[2], &cells);
        if (cell->count == cell->asize) {
            if (cell->asize < 1000000)
                cell->asize *= 2;
            else
                cell->asize += cell->asize / 2;
            cell->atoms = (dbr_xyz*) xrealloc(cell->atoms,
                                              cell->asize * sizeof(dbr_xyz));
            if (xa->indices)
                cell->indices = (int*) xrealloc(cell->indices,
                                              cell->asize * sizeof(int));
        }
        memcpy(cell->atoms[cell->count], xa->xyz[i], sizeof(dbr_xyz));
        if (xa->indices)
            cell->indices[cell->count] = xa->indices[i];
        cell->count++;
    }
    /* shrink cell->atoms of real cells */
    for (i = 0; i < cells.count; ++i) {
        dbr_cell *cell = &cells.data[i];
        if (!cell->real)
            continue;
        cell->asize = cell->count;
        if (cell->asize == 0) {
            free(cell->atoms);
            cell->atoms = 0;
        }
        else {
            cell->atoms = (dbr_xyz*) xrealloc(cell->atoms,
                                              cell->asize * sizeof(dbr_xyz));
        }
    }

    free(xa->xyz);
    free(xa->indices);
    if (dbr_verbosity > 1) /* very verbose */
        print_cells_memory(cells);

    /* put atoms into virtual cells */
    for (i = 0; i < cells.a[0]; ++i) {
        for (j = 0; j < cells.a[1]; ++j) {
            for (k = 0; k < cells.a[2]; ++k) {
                if (((i == 0 || i == cells.n[0] + 1) && cells.v[0])
                    || ((j == 0 || j == cells.n[1] + 1) && cells.v[1])
                    || ((k == 0 || k == cells.n[2] + 1) && cells.v[2])) {
                    /* find original cell indices */
                    dbr_real d[3];
                    n[0] = i;
                    n[1] = j;
                    n[2] = k;
                    for (m = 0; m < 3; ++m) {
                        if (cells.v[m] && n[m] == 0)
                            n[m] = cells.n[m];
                        else if (cells.v[m] && n[m] == cells.n[m] + 1)
                            n[m] = 1;
                    }

                    d[0] = (i - n[0]) * cells.vectors.v00
                         + (j - n[1]) * cells.vectors.v10
                         + (k - n[2]) * cells.vectors.v20;
                    d[1] = (i - n[0]) * cells.vectors.v01
                         + (j - n[1]) * cells.vectors.v11
                         + (k - n[2]) * cells.vectors.v21;
                    d[2] = (i - n[0]) * cells.vectors.v02
                         + (j - n[1]) * cells.vectors.v12
                         + (k - n[2]) * cells.vectors.v22;
                    make_virtual_cell(get_abs_cell(i, j, k, &cells),
                                      get_abs_cell(n[0], n[1], n[2], &cells),
                                      d[0], d[1], d[2]);
                }
            }
        }
    }
    return cells;
}

void free_cells(dbr_cells cells)
{
    int i;
    for (i = 0; i < cells.count; ++i)
        free(cells.data[i].atoms);
    free(cells.data);
}


void dbr_print_version()
{
    if (dbr_nid != 0)
        return;
    printf("debyer %s (%s, %s precision%s%s)\n",
           VERSION,
#ifdef USE_MPI
          "MPI",
#else
          "serial",
#endif
          sizeof(dbr_real) == sizeof(float) ?  "single" : "double",
#ifdef HAVE_ZLIB
          ", handles .gz",
#else
          "",
#endif
#ifdef HAVE_BZLIB
          ", handles .bz2"
#else
          ""
#endif
        );
}

