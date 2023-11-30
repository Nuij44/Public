exe : main.o outil.o solex.o coefficients.o frontiere.o fredholm.o
	gcc -O2 main.o outil.o solex.o coefficients.o frontiere.o fredholm.o

outil.o : outil.c
	gcc -c outil.c

solex.o : solex.c
	gcc -c solex.c

coefficients.o : coefficients.c outil.h solex.h
	gcc -c coefficients.c

frontiere.o : frontiere.c outil.h coefficients.h solex.h
	gcc -c frontiere.c

fredholm.o : fredholm.c outil.h
	gcc -c fredholm.c

main.o : main.c outil.h coefficients.h solex.h frontiere.h fredholm.h
	gcc -c main.c