# ------------------- Variables, environment  ----------------
SHELL	= /bin/sh
MODULEDIR1 = include
MODULEDIR2 = ../include

# ------------------- Variables, compiler     ----------------

FC	= ifort
OMP     = -openmp
#CHECKS = -check all
CHECKS =
PROFILING= # -pg
FFLAGS	= -O -w $(OMP)  -module $(MODULEDIR1) -fpp $(CHECKS) $(PROFILING)
LDFLAGS	= 

# ------------------- Variables, libraries    ----------------

NRECIPES = ../nr/librecipes_f90.a

CTPPPROC = ../lib/libctppproc.a

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

#LAPACK  = -I$(MKLROOT)/include $(MKLPATH)/libmkl_solver_lp64.a -Wl, --start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl, --end-group  -lmkl_lapack
LAPACK = -I$(MKLROOT)/include -mkl


LIB 	= -Vaxlib $(LAPACK) -lpthread  -L$(HEALPIX)/lib -lhealpix -lcfitsio  -lfftw3 

INCLUDE = -I$(HEALPIX)/include -I$(MODULEDIR2)

OBJ    	= Topology_types.o Topology_map_mod.o nr_minimization.o lm_rotate.o ctpp_eigen_mod.o basis_modes.o Topology_Lmarg_mod.o Topology_Lmarg.o

OBJMAP  = Topology_types.o Topology_map_mod.o Topology_make_map.o lm_rotate.o ctpp_eigen_mod.o basis_modes.o
# ------------------------- Rules ----------------------------

default: all

all: topmarg topmap

topmarg: $(OBJ) $(NRECIPES) $(CTPPPROC)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIB) $(NRECIPES) $(CTPPPROC)

topmap: $(OBJMAP) $(NRECIPES) $(CTPPPROC)
	$(FC) $(FFLAGS) -o $@ $(OBJMAP) $(LIB) $(NRECIPES) $(CTPPPROC)

Topology_Lmarg.o      : Topology_types.o Topology_Lmarg_mod.o \
			Topology_map_mod.o basis_modes.o
Topology_make_map.o   : Topology_types.o Topology_map_mod.o
Topology_Lmarg_mod.o  : Topology_types.o ctpp_eigen_mod.o nr_minimization.o
Topology_map_mod.o    : Topology_types.o 
ctpp_eigen_mod.o      : Topology_types.o basis_modes.o
basis_modes.o         : Topology_types.o

%.o    : %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@

# ----------------------- Clean up ----------------------------
 
.PHONY : tidy clean cleanall
tidy:
	-rm -f $(OBJ)
clean: tidy
	-rm -f topmarg
cleanall: clean
	/bin/rm -f $(MODULEDIR1)/[!n]*.mod

