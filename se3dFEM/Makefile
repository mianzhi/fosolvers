FC=gfortran
FFLAGS=-Wall -fopenmp -O2 -mtune=native
objects=se3dFEM.o fileio.o stima.o asm.o misc.o

%.o:%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

se3dFEM:$(objects)
	$(FC) $(FFLAGS) -o se3dFEM $(objects)

clean:
	rm *.o *.mod

install:
	./installer.sh
