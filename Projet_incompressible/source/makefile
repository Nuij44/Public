exe : main.o recup_affichage.o decla_type.o mef.o
	gfortran -O2 -fbacktrace --warn-all main.o recup_affichage.o decla_type.o mef.o

main.o : main.f90 recup_affichage.o decla_type.o mef.o
	gfortran -c main.f90

recup_affichage.o : recup_affichage.f90	
	gfortran -c recup_affichage.f90

decla_type.o : decla_type.f90
	gfortran -c decla_type.f90

mef.o : mef.f90
	gfortran -c mef.f90