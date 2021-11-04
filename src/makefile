#-------------------------------------------------------------------------------
#
# Default host compiler options
#
#-------------------------------------------------------------------------------

FC            = gfortran
FFLAGS        = -O3 -DGNU # -fmax-errors=5
#FFLAGS        = -O0 -DGNU -Wall -fbacktrace -fcheck=all -g #-ffpe-trap=invalid,zero,overflow 
LD            = $(FC)
LDFLAGS       = 
LIBS          =


#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FFLAGS) $*.F90

#-------------------------------------------------------------------------------

name = NStar

default: $(name).x

sources = constants.F90 parameters.F90 EoS_NS.F90 minpackMod1.F90 mesh.F90        \
          variables.F90 boundaries.F90 init_NS.F90 physics.F90 fluxsplit.F90   \
          reconstr.F90 num_flux.F90 sources.F90 evolution.F90                  \
          evol_sm_euler.F90 evol_sm_rgk2.F90 evol_sm_rgk3.F90 evol_sm_rgk4.F90 \
          writedata.F90 main.F90

objects = constants.o parameters.o EoS_NS.o minpackMod1.o mesh.o variables.o      \
          boundaries.o init_NS.o physics.o fluxsplit.o reconstr.o num_flux.o   \
          sources.o evolution.o evol_sm_euler.o evol_sm_rgk2.o                 \
          evol_sm_rgk3.o evol_sm_rgk4.o writedata.o main.o 

files   = $(sources) makefile params.in 



$(name).x: $(objects)
	$(LD) $(LDFLAGS) $(objects) -o $(name).x $(LIBS)

dist: $(files)
	tar czvf $(name).tar.gz $(files)

clean-bak:
	rm -rf *.bak *~

clean-data:
	rm -rf *.dat info.txt ./data/*.dat

clean-exec:
	rm -rf *.x

clean-modules:
	rm -rf *.mod *.smod

clean-objects:
	rm -rf *.o

clean: clean-bak clean-modules clean-objects

clean-all: clean-bak clean-data clean-exec clean-modules clean-objects
