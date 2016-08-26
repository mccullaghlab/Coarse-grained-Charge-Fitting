
f90 = gfortran
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none
lapacklib = /usr/local/lib

density : analytic_charge_min.f90 dcd.f90 stringlib.f90
	$(f90) -c analytic_charge_min.f90 dcd.f90 stringlib.f90 $(flags) -L$(lapacklib) -llapack -lblas
	$(f90)  analytic_charge_min.o dcd.o stringlib.o -o analytic_charge_min.x  $(flags) -L$(lapacklib) -llapack -lblas

clean:
	rm -f *.o *.mod *.x

