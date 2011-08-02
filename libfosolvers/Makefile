objects=grid.o fileio.o misc.o
target=libfosolvers.a

FC=gfortran
FFLAGS=-Wall -fopenmp -O2 -mtune=native

%.o:%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

$(target):$(objects)
	ar cr $(target) $(objects)

clean:
	rm *.o *.mod
