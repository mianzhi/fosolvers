objects=demo.o
libpath=../
libs=-lfosolvers -llapack
target=demo

FC=mpif90
FFLAGS=-Wall -fopenmp -O3 -mtune=native -I$(libpath) -L$(libpath)

%.o:%.f90
	$(FC) -c $(FFLAGS) -o $@ $^

$(target):$(objects)
	$(FC) $(FFLAGS) -o $(target) $(objects) $(libs)

clean:
	rm *.o *.mod
