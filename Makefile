# ------------------- Variables, environment  ----------------
SHELL	= /bin/sh
MODULEDIR1 = include
MODULEDIR2 = ../include

# ------------------- Variables, compiler     ----------------

FC	= ifort
OMP     = -openmp
#FFLAGS	= -O -w -openmp -static-libcxa -module $(MODULEDIR1) -fpp -qp #-r8 -g -qp
PROFILING= -pg
FFLAGS	= -O -w $(OMP)  -module $(MODULEDIR1) -fpp #-r8 -g -qp
#FFLAGS	= -O -w $(OMP)  -module $(MODULEDIR1) -fpp $(PROFILING)  #-r8 -g -qp
LDFLAGS	= 

# ------------------- Variables, libraries    ----------------

NRECIPES = ../nr/librecipes_f90.a

# MKL_ROOT is set by intel shell script, source it in .bashrc
# source /opt/intel/Compiler/$(INTEL_VERSION)/mkl
# 
# or for separate MKL installation one has to do it by hand
# MKL_VERSION=9.1
# MKLROOT=/opt/intel/mkl/${MKL_VERSION}
# 
MKLARCH = em64t
MKLPATH=$(MKLROOT)/lib/$(MKLARCH)
	                  # Other choices:  32 64

FITSOUT	= 

# ------------------------------------------------------------

# For old standalone MKL
# LAPACK  =  -I$(MKLROOT)/include -L$(MKLROOT)/lib/$(MKLARCH) -lmkl_lapack -lmkl_$(MKLARCH)

LAPACK  = -I$(MKLROOT)/include $(MKLPATH)/libmkl_solver_lp64.a -Wl, --start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl, --end-group  -lmkl_lapack


LIB 	= -Vaxlib $(LAPACK) -lguide -lpthread  -L$(HEALPIX)/lib -lhealpix -lcfitsio  -lfftw3 

INCLUDE = -I$(HEALPIX)/include -I$(MODULEDIR2)

OBJ    	= nml_mod.o Topology_types.o Topology_map_mod.o Topology_map_mod_nel.o nr_minimization.o lm_rotate.o ctpp_eigen_mod.o Topology_Lmarg_mod.o Topology_Lmarg.o

# ------------------------- Rules ----------------------------

default: topmarg

%.o:%.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@

topmarg: $(OBJ) $(NRECIPES)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIB) $(NRECIPES)

clean:
	-rm -f $(OBJ)
tidy: clean
	-rm -f topmarg
