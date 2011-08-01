objects=libtest.o
libpath=../libfosolvers/
libs=libfosolvers.a
target=libtest

FC=gfortran
FFLAGS=-Wall -fopenmp -O2 -mtune=native -I$(libpath)

%.o:%.f90 $(libpath)$(libs)
	$(FC) -c $(FFLAGS) -o $@ $^

$(target):$(objects) $(libpath)$(libs)
	$(FC) $(FFLAGS) -o $(target) $(objects) $(libpath)$(libs)

clean:
	rm *.o *.mod

install:
	./installer.sh
