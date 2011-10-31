objects=grid.o fileio.o mp.o misc.o
target=libfosolvers.a

FC=mpif90
FFLAGS=-Wall -fopenmp -O2 -mtune=native

%.o:%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

$(target):$(objects)
	ar cr $(target) $(objects)

clean:
	rm *.o *.mod
