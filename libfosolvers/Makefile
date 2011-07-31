FC=gfortran
FFLAGS=-Wall -fopenmp -O2 -mtune=native
objects=grid.o
target=libfosolvers.a

%.o:%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

$(target):$(objects)
	ar cr $(target) $(objects)

clean:
	rm *.o *.mod
