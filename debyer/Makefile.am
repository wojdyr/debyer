EXTRA_DIST = debyer.ggo extend_cmd.ggo

noinst_LIBRARIES = libdebyer.a

BUILT_SOURCES = extend_cmd.c cmdline.c

bin_PROGRAMS = debyer dbr_conv dbr_extend dbr_bonds sic_analysis dbr_symmetry

# programs that will likely be deleted in future
noinst_PROGRAMS = sic_orient dbr_volume dbr_diffus

libdebyer_a_SOURCES = debyer.c lineio.cc fileio.cc atomtables.c \
 		      debyer.h lineio.h fileio.h  iloops.h atomtables.h

debyer_SOURCES = cmdline.c cmdline.h main.cc
debyer_LDADD = libdebyer.a

sic_analysis_SOURCES = sic_analysis.cc
sic_analysis_LDADD = libdebyer.a

sic_orient_SOURCES = sic_orient.cc utils.h
sic_orient_LDADD = libdebyer.a

dbr_diffus_SOURCES = dbr_diffus.cc utils.h
dbr_diffus_LDADD = libdebyer.a

dbr_volume_SOURCES = dbr_volume.cc
dbr_volume_LDADD = libdebyer.a

dbr_symmetry_SOURCES = dbr_symmetry.cc
dbr_symmetry_LDADD = libdebyer.a

dbr_extend_SOURCES = dbr_extend.cc extend_cmd.c extend_cmd.h cells.cc cells.h
dbr_extend_LDADD = libdebyer.a

dbr_bonds_SOURCES = dbr_bonds.cc cells.cc cells.h
dbr_bonds_LDADD = libdebyer.a

dbr_conv_SOURCES = conv.cc
dbr_conv_LDADD = libdebyer.a

cmdline.c: debyer.ggo
	gengetopt -i debyer.ggo

extend_cmd.c: extend_cmd.ggo
	gengetopt -i extend_cmd.ggo

