objects=misc.o utility.o grid.o fileio.o
target=libfosolvers.a

FC=mpif90
FFLAGS=-Wall -fopenmp -O3 -mtune=native

%.o:%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

$(target):$(objects)
	ar cr $(target) $(objects)

clean:
	rm *.o *.mod
