AUTOMAKE_OPTIONS = serial-tests
AM_CFLAGS = -I $(top_srcdir)/include

TESTS_ENVIRONMENT = $(srcdir)/runtest.sh
EXTRA_DIST = runtest.sh
TESTS = \
	hsl_mi20ss \
	hsl_mi20ds \
	hsl_mi20ss1 \
	hsl_mi20ds1

check_PROGRAMS = \
	hsl_mi20ss \
	hsl_mi20ds \
	hsl_mi20ss1 \
	hsl_mi20ds1

hsl_mi20ss_SOURCES = hsl_mi20ss.c ../../include/hsl_mi20s.h
hsl_mi20ss_LDADD = $(top_builddir)/src/libhsl_mi20.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_mi20ss.output
hsl_mi20ds_SOURCES = hsl_mi20ds.c ../../include/hsl_mi20d.h
hsl_mi20ds_LDADD = $(top_builddir)/src/libhsl_mi20.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_mi20ds.output
hsl_mi20ss1_SOURCES = hsl_mi20ss1.c ../../include/hsl_mi20s.h
hsl_mi20ss1_LDADD = $(top_builddir)/src/libhsl_mi20.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_mi20ss1.output
hsl_mi20ds1_SOURCES = hsl_mi20ds1.c ../../include/hsl_mi20d.h
hsl_mi20ds1_LDADD = $(top_builddir)/src/libhsl_mi20.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_mi20ds1.output
