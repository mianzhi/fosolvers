FC=gfortran
FFLAGS=-Wall -fopenmp -O2 -mtune=native
objects=heat3dFVM.o fileio.o tadv.o misc.o

%.o:%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

heat3dFVM:$(objects)
	$(FC) $(FFLAGS) -o heat3dFVM $(objects)

clean:
	rm *.o *.mod

install:
	./installer.sh
