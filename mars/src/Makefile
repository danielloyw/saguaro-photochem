UNAME := $(shell uname)
HOSTNAME := $(shell hostname)

# --- Linux Compiler Variables ------------------------------------------------
ifeq ($(UNAME), Linux)
FC = ifx
WARN = -warn all -check all
#FC = gfortran
#WARN =  -g -Wall -Wextra -fcheck=all

# MKL variables
FCFLAGS = -I${MKLROOT}/include/mkl/intel64/lp64
LDFLAGS = -qmkl -lmkl_lapack95_lp64
endif


# --- MacOS Compiler Variables ------------------------------------------------
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

# --- Project Files -----------------------------------------------------------
#
FILES = types.f90 constants.f90 global_variables.f90 utils.f90 \
	read_inputs.f90 diffusion.f90 photo.f90 electrons.f90 \
	output.f90 solve.f90 saguaro.f90
OBJS=$(FILES:%.f90=%.o)

all: saguaro-mars

saguaro-mars: $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) -o saguaro-mars $(LDFLAGS)
	mv saguaro-mars ../wrk/

test: $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) -o test $(LDFLAGS)

%.o : %.f*
	$(FC) $(WARN) $(FCFLAGS) -c -o $@ $< $(LDFLAGS)


# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ *.x
