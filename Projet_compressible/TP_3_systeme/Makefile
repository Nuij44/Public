a.out : main.o fonction.o
	gfortran -o exe main.o fonction.o
	
main.o : main.f90 fonction.o
	gfortran -c main.f90

fonction.o : fonction.f90
	gfortran -c fonction.f90

