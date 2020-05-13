INCLUDES = /contrib/include
CPPFLAGS = -I$(INCLUDES)
CXX = g++
LIBS = -L/contrib/lib -lnetcdf
CXXFLAGS = -g -Wall

.SUFFIXES: .C
.C.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $<

all: post

post: ncfile.o post.o field.o 
	g++ ncfile.o  post.o field.o $(LIBS) -o post 

ncfile.o: field.H ncfile.H

post.o: field.H ncfile.H 

field.o: field.H

clean:
	rm *.o post