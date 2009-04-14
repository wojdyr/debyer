EXTRA_DIST = debyer.ggo

noinst_LIBRARIES = libdebyer.a 

bin_PROGRAMS = debyer sic_analysis sic_orient dbr_diffus dbr_volume \
	       dbr_symmetry

libdebyer_a_SOURCES = debyer.c lineio.cc fileio.cc crystables.c \
 		      debyer.h lineio.h fileio.h  iloops.h crystables.h 

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


cmdline.c: debyer.ggo
	gengetopt -uFILE --no-handle-version < debyer.ggo

