EXTRA_DIST = debyer.ggo extend_cmd.ggo

noinst_LIBRARIES = libdebyer.a
LDADD = libdebyer.a $(OPENMP_CFLAGS)

BUILT_SOURCES = extend_cmd.c cmdline.c

bin_PROGRAMS = debyer dbr_conv dbr_extend dbr_bonds dbr_sic dbr_symmetry

# programs that will likely be deleted in future
#noinst_PROGRAMS = dbr_volume dbr_diffus

libdebyer_a_SOURCES = debyer.c lineio.cc fileio.cc atomtables.c \
		      debyer.h lineio.h fileio.h atomtables.h
libdebyer_a_CFLAGS = $(OPENMP_CFLAGS)

debyer_SOURCES = cmdline.c cmdline.h main.cc

dbr_sic_SOURCES = sic.cc utils.h

#dbr_diffus_SOURCES = dbr_diffus.cc utils.h

#dbr_volume_SOURCES = dbr_volume.cc utils.h

dbr_symmetry_SOURCES = dbr_symmetry.cc utils.h

dbr_extend_SOURCES = dbr_extend.cc extend_cmd.c extend_cmd.h cells.cc cells.h \
		     utils.h

dbr_bonds_SOURCES = dbr_bonds.cc cells.cc cells.h

dbr_conv_SOURCES = conv.cc

cmdline.c: debyer.ggo
	gengetopt -i $<

extend_cmd.c: extend_cmd.ggo
	gengetopt -i $<

