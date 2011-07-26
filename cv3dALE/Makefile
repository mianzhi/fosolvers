main:main.o cv3dALE.o utility.o misc.o fileio.o
	gfortran -fopenmp -O2 main.o cv3dALE.o utility.o misc.o fileio.o -o main

main.o:main.f90
	gfortran -c -fopenmp -O2 main.f90
cv3dALE.o:cv3dALE.f90
	gfortran -c -fopenmp -O2 cv3dALE.f90
utility.o:utility.f90
	gfortran -c -fopenmp -O2 utility.f90
misc.o:misc.f90
	gfortran -c -fopenmp -O2 misc.f90
fileio.o:fileio.f90
	gfortran -c -fopenmp -O2 fileio.f90
	
clean:
	rm main.o cv3dALE.o utility.o misc.o fileio.o globalvar.mod
