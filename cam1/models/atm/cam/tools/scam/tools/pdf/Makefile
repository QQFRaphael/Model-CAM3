all: pdf

pdf: pdf.o
	g++  pdf.o  -L/contrib/lib -lnetcdf -o pdf

pdf.o: pdf.C
	g++ -O -c -Wall -I/contrib/include pdf.C


