OPTIM    = no
PARALLEL = yes
nb_procs = 4
###################################################################
###################################################################
# python
PC = python3
###################################################################
###################################################################
# Fortran compiler and flags
ifeq ($(PARALLEL),yes)
	FC = mpif90
	FFLAGS_MPI = -specs override-specs.h -Ddo_mpi
else
	FC = gfortran
	FFLAGS_MPI = -specs override-specs.h
endif

ifeq ($(OPTIM),yes)
	FFLAGS = -O3 $(FFLAGS_MPI)
else
	FFLAGS = -O0 -g -fbounds-check -Wfatal-errors -Wall -Wconversion -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-local-zero $(FFLAGS_MPI) #Wall -Wextra
endif
###################################################################
###################################################################
# for f2py
#F2PYFLAGS = --help-fcompiler
#F2PYFLAGS = --no-lower --quiet --f90flags="$(FFLAGS)" --fcompiler=gnu95
F2PYFLAGS = --no-lower --quiet --f90flags="$(FFLAGS)" --f77exec=$(FC) --f90exec=$(FC)
###################################################################
###################################################################
# source files and targets
MOD = propag
FSRC = maths.f90 outdata.f90 intderiv.f90 dataspace.f90 findiff.f90 modelspace.f90 propag.f90 linop.f90
#FSRC = foo.f90
OBJ  = $(FSRC:.f90=.o)
###################################################################
###################################################################
# To compile Fortran subroutine in pure Fortran
all: main

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

main: $(OBJ)  main.o
	$(FC) $(FFLAGS) -o $@ $^

runf:
ifeq ($(PARALLEL),yes)
	time mpirun -n $(nb_procs) -host localhost ./main
else
	time ./main	
endif
###################################################################
###################################################################
# To compile with f2py
fpy:
	f2py -c $(FSRC) -m $(MOD) $(F2PYFLAGS)

run:
ifeq ($(PARALLEL),yes)
	mpirun -np $(nb_procs) -host localhost $(PC) main.py
else
	$(PC) main.py
endif
###################################################################
###################################################################
.PHONY: clean cleangraph cleanall

cleanall: clean cleangraph

clean:
	rm -f main *.mod *.o *~ *.txt~ *.pyf *.so *.npy

cleangraph:
	rm -f out/*.pckl out/*.npy
