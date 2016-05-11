
prefix = /usr/local/bin
mpif90 = /opt/local/lib/openmpi/bin/mpif90
f90 = gfortran
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none
lapacklib = /usr/local/lib

density : analytic_charge_min.f90 dcd.f90 stringlib.f90
	$(f90) -c analytic_charge_min.f90 dcd.f90 stringlib.f90 $(flags) -L$(lapacklib) -llapack -lblas
	$(f90)  analytic_charge_min.o dcd.o stringlib.o -o analytic_charge_min.x  $(flags) -L$(lapacklib) -llapack -lblas
#	cp fit_charges.x /Volumes/LaCie/voth/hydrophobic_fm/villin/charge_fitting

clean:
	rm -f *.o *.mod *.x

