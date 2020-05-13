all: diurnal_ave.o
	g++  diurnal_ave.o  -L/contrib/lib -lnetcdf -o diurnal_ave

diurnal_ave.o: diurnal_ave.C
	g++ -O2 -c -Wall -I/contrib/include diurnal_ave.C

