exe: test.o fredholm.o outil.o
	gcc -O2 test.o fredholm.o outil.o

test.o : test.c fredholm.h outil.h
	gcc -c test.c

fredholm.o : fredholm.c outil.h
	gcc -c fredholm.c

outil.o : outil.c
	gcc -c outil.c