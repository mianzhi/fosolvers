objects=demo.o
libpath=../
libs=-lfosolvers -llapack
target=demo

FC=gfortran
FFLAGS=-Wall -fopenmp -O2 -mtune=native -I$(libpath) -L$(libpath)

%.o:%.f90
	$(FC) -c $(FFLAGS) -o $@ $^

$(target):$(objects)
	$(FC) $(FFLAGS) -o $(target) $(objects) $(libs)

clean:
	rm *.o *.mod
