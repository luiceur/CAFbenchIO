FC=		ftn
FFLAGS=		-c 

EXECUTABLE=benchCAFio
# level 1
L1_SRC=		cgca_m1co.f90 benchclock.f90
L1_OBJ=		$(L1_SRC:.f90=.o)

# level 2 modules and submodules

L2_SRC_MOD=     cgca_m2netcdf.f90 cgca_m2mpiio.f90 cgca_m2alloc.f90 \
	 	cgca_m2phys.f90 cgca_m2hdf5.f90
L2_SRC=		$(L2_SRC_MOD) 
L2_OBJ=		$(L2_SRC:.f90=.o)

LTOP_SRC_MOD=	drive.F90
LTOP_OBJ=	$(LTOP_SRC_MOD:.F90=.o)

ALL_OBJ=	$(L1_OBJ) $(L2_OBJ) $(LTOP_OBJ)

.SUFFIXES: .F90 .f90 .o

.f90.o:
	$(FC) $(FFLAGS) $<

.F90.o:
	$(FC) $(FFLAGS) $<

all: $(ALL_OBJ) $(EXECUTABLE)

$(EXECUTABLE): $(ALL_OBJ) 
	$(FC)  $(ALL_OBJ) -o $@

clean:
	rm -f *~ *.o *.mod *.lst ${EXECUTABLE}