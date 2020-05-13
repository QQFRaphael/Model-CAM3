#ifndef FIELD_H
#define FIELD_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

class Field;

#include "ncfile.H"

#define  MAX_FIELDS 1000

class TimeSlice {
    friend class Field;
    float* data;
    int size;
    int ntime;
    enum { bumpsize = 10 };    // amount to increase size in Resize()
public:
    TimeSlice() { 
        ntime = 1;
        size = 1;
        data = new float[ntime];
    }

    TimeSlice(int ntime) { 
        data = new float[ntime];
        this->size = ntime;
        this->ntime = ntime;
    }
    TimeSlice operator=(TimeSlice& t) {
//        cout << "TimeSlice() operator= "<< endl ;
        size = t.size;
        ntime = t.ntime;
        if ( data != NULL )     
            delete[] data;
        data = new float[size];
        for ( int i=0; i<size; i++ ) 
            data[i] = t.data[i];
        return *this;
    }
    
    TimeSlice(TimeSlice& t) {
//        cout << "TimeSlice() copy constructor = "<< endl ;
        size = t.size;
        ntime = t.ntime;
        data = new float[size];
        for ( int i=0; i<size; i++ ) 
            data[i] = t.data[i];
    }
    
    ~TimeSlice() { 
//        cout<<"delete ~TimeSlice "<<endl; 
        delete[] data;
    }
    void Resize() {
//        cout << "TimeSlice::resize()" << endl;
        float *tmp = data;
        data = new float[size+bumpsize];
        for ( int i=0;i<size;i++ )
            data[i] = tmp[i];
        delete[] tmp;
        size += bumpsize;
    }

    float& operator[](int i) { 
        if ( i<0 ) {
            cerr<<"error: level - illegal index " << i << endl;
            exit ( 1 );
        }
        if ( i >= ntime ) {
            ntime = i+1;
            while ( i >= size )
                Resize();
        }
        return data[i];
    }
};

class Field
{
    friend ostream& operator<<(ostream& o, Field& f );
public:
    Field( int nlev, int ntime );
    Field( char* descrip, char* format, int nlev, int ntime, int id, float unitConversion = 1.0, int group = 0 );
    const char* GetDescrip() { return descrip; }
    void   WriteData( ofstream& o, int time );
    void   WriteData( ofstream& o );
    int    GetGroup() { return group; }
    int    GetId() { return id; }
    int    GetNlev() { return nlev; }
    int    GetNtime() { return ntime; }
    void   Resize();
    void   SetDataIsMissing( bool isMissing );
    bool   DataIsMissing() { return dataIsMissing; }
    void   SetMissingValueString();
    TimeSlice&  operator[](int i); 
    Field  operator=( Field& f );
private:
    enum { bumpsize = 10 };
    TimeSlice* timedim;
    char   descrip[256];
    int    id;
    int    group;
    float  unitConversion;
    char   format[10];
    char   missingValue[256];
    int    size;
    int    nlev;
    int    ntime;
    bool   dataIsMissing;
};

#endif // FIELD_H
