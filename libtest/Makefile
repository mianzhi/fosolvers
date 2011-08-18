objects=libtest.o
libpath=../libfosolvers/
libs=-lfosolvers -llapack
target=libtest

FC=gfortran
FFLAGS=-Wall -fopenmp -O2 -mtune=native -I$(libpath) -L$(libpath)

%.o:%.f90
	$(FC) -c $(FFLAGS) -o $@ $^

$(target):$(objects)
	$(FC) $(FFLAGS) -o $(target) $(objects) $(libs)

clean:
	rm *.o *.mod

install:
	./installer.sh
