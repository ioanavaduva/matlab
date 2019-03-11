lib_LIBRARIES = libhsl_mi20.a

modulesdir = modules
nodist_modules_HEADERS = \
	hsl_mi20_single.mod \
	hsl_mi20_double.mod

libhsl_mi20_a_SOURCES = \
	hsl_mi20s.f90 \
	hsl_mi20d.f90 \
	common90.f90 \
	sdeps.f \
	sdeps90.f90 \
	ddeps.f \
	ddeps90.f90 \
	hsl_mi20s_ciface.f90 \
	hsl_mi20d_ciface.f90

hsl_mi20s.$(OBJEXT): sdeps90.$(OBJEXT)
hsl_mi20s.$(OBJEXT): common90.$(OBJEXT)
hsl_mi20d.$(OBJEXT): ddeps90.$(OBJEXT)
hsl_mi20d.$(OBJEXT): common90.$(OBJEXT)
sdeps90.$(OBJEXT): common90.$(OBJEXT)
ddeps90.$(OBJEXT): common90.$(OBJEXT)
hsl_mi20s_ciface.$(OBJEXT): hsl_mi20s.$(OBJEXT)
hsl_mi20d_ciface.$(OBJEXT): hsl_mi20d.$(OBJEXT)
