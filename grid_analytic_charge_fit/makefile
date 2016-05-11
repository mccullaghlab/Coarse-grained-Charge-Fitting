
prefix = /usr/local/bin
mpif90 = /opt/local/lib/openmpi/bin/mpif90
f90 = gfortran
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none
lapacklib = /usr/local/lib

density : grid_esp.f90 dcd.f90 stringlib.f90
	$(f90) -c grid_esp.f90 dcd.f90 stringlib.f90 $(flags) -L$(lapacklib) -llapack -lblas
	$(f90)  grid_esp.o dcd.o stringlib.o -o grid_esp.x  $(flags) -L$(lapacklib) -llapack -lblas
#	cp fit_charges.x /Volumes/LaCie/voth/hydrophobic_fm/villin/charge_fitting

clean:
	rm -f *.o *.mod *.x

