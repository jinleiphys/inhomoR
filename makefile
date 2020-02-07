OBJ= precision.o  constants.o gauss.o interpolation.o mesh.o pot.o channels.o input.o  \
coul90.o lagrange_mesh_source.o green.o rmat.o numerov.o inhomoR.o

BASE := $(shell expr $(CURDIR) : "\(.*\)/.*")
VERDATE := $(shell git log -1 --format=%cd  )
# VERDATE := $(shell expr "$(VERDATE)" : '\(.*\)(.*).*')
VERREV := $(shell git log -1 --pretty=format:"%h")
COMPDATE :=$(shell date)




# Laptop
 LIB =-L  /opt/local/lib -llapack
#LIB = -L ../lapack-3.5.0 -lrefblas -llapack
FC = ifort
F90 = ifort
FFLAGS = -O3 -132

.SUFFIXES: .F90 .f90 .f95

all: inhomo


inhomo:  $(OBJ)
	$(FC) -o inhomoR $(FFLAGS) $(OBJ) $(LIB)




.F90.o          :
	$(F90) $(FFLAGS) -c $<

.f95.o          :
	$(F90) $(FFLAGS) -c $<

.F.o          :
	$(F90) $(FFLAGS) -c $<

.f.o          :
	$(F90) $(FFLAGS) -c $<


clean:
	rm -f *.mod core *.o inhomoR
