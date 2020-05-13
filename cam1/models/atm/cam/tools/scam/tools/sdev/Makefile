all: sdev.o
	g++  sdev.o  -L/contrib/lib -lnetcdf -o sdev

sdev.o: sdev.C
	g++ -g -c -Wall -I/contrib/include sdev.C

