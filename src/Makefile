#
# example makefile to compile
#


UNAME := $(shell uname)
HOSTNAME := $(shell hostname)

#
# set compiler variables for clusters
#
ifeq ($(UNAME), Linux)
FC = ifort
FCFLAGS = -O3 -lmkl:sequential
LDFLAGS = -lmkl_lapack95_lp64
# since hipas is missing this include folder for mkl
ifeq ($(HOSTNAME), hipas.lpl.arizona.edu) 
FCFLAGS += -I$(MKLROOT)/include/em64t/lp64
endif
endif


#
# set compiler variables for Mac
#
ifeq ($(UNAME), Darwin)
FC = ifort
#FDEBUG = -warn all -debug -g -fp-stack-check -traceback -gen-interfaces
# FDEBUG=-warn all
# FDEBUG+=-debug -g -fp-stack-check -traceback -gen-interfaces
 #FDEBUG+= -fstack-protector 
 #FDEBUG+= -assume protect_parens 
# FDEBUG+= -implicitnone
 #FDEBUG+= -check
# FDEBUG+= -check uninit
# FDEBUG=-warn all -debug -g -fp-stack-check -check all -traceback -gen-interfaces
#FCFLAGS= -O0 -g -check all -fpe0 -warn -traceback -nogen-interface -mkl=sequential
#FCFLAGS = -O3 -heap-arrays -mkl=sequential
FCFLAGS = -O -heap-arrays -mkl=sequential
#FCFLAGS = -C -heap-arrays -mkl=sequential
#FCFLAGS = -heap-arrays -mkl=sequential
LDFLAGS = -lmkl_lapack95_lp64
FCFLAGS+=$(FDEBUG)
endif

#
# OK, now input all the files to be
#
FILES = precision.f90 constants.f90 global_variables.f90 subs.f90 find_name.f90 vapor.f90 \
	locate.f90 tridag.f90 lubksb.f90 ludcmp.f90 intrp.f90 indexx.f90 wait.f90 read_molecules.f90 \
	read_reactions.f90 read_atmos.f90 read_photo.f90 read_photoB.f90 read_jvals.f90 \
	read_sol.f90 read_sol_hres.f90 norm_sol.f90 \
	read_aerosol.f90 read_rayleigh.f90 ext_prod.f90 subs_photo.f90 photo.f90 paths1D.f90 \
	rateco.f90 difco.f90 soldep1.f90 subs_compo.f90 compo.f90 chemeq.f90 compout.f90 \
	hydrost.f90 solveit.f90 find_bin.f90 difcs_sec.f90 difcs_rees.f90 eldep1.f90 elctrn.f90 pchem1D.f90
OBJS=$(FILES:%.f90=%.o)



all: pchem.x

pchem.x: $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) -o pchem.x $(LDFLAGS); mv pchem.x ../wrk


# make objects

#define dependencies
compo.o: compo.f90 vapor.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

%.o : %.f*
	$(FC) $(FCFLAGS) -c -o $@ $<


# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ *.x
