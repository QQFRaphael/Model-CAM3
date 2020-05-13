#include <string.h>
#include "field.H"

#define PR(x) cout << #x " = " << x << endl;

Field::Field( char* descrip, char* format, int nlev, int ntime, int id, float unitConversion, int group )
{
    strcpy( this->descrip, descrip );
    strcpy( this->format, format );
    this->nlev = nlev;
    this->ntime = ntime;
    this->group = group;
    this->id = id;
    size = nlev;
    this->unitConversion = unitConversion;
    timedim = new TimeSlice[nlev];
    dataIsMissing = false;
//    cout << "new field: " <<  *this << endl;
}

Field::Field( int nlev, int ntime )
{
    strcpy( descrip, "n/a" );
    strcpy( format, "" );
    this->nlev = nlev;
    this->ntime = ntime;
    id = -1;
    size = nlev;
    unitConversion = 1;
    timedim = new TimeSlice[nlev];
    dataIsMissing = false;
}

ostream& 
operator<<(ostream& o, Field& f ) 
{
    int nlev = f.GetNlev();
    int ntime = f.GetNtime();
    o << "Desc = " << f.GetDescrip() << endl;
    o << "Id = " << f.GetId() << endl;
    o << "Nlev = " << nlev << endl;
    o << "Ntime = " << ntime << endl;
//      o << "Data: " << endl;
//      for ( int t=0; t<ntime; t++ ) {
//          for ( int l=0; l<nlev; l++ ) {
//              o << f[l][t] << ", ";
//          }
//          o << "\n\n";
//      }

    o << "\n -------------------------\n" << endl;
    return o;
}


TimeSlice& 
Field::operator[](int i)
{
    if ( i<0 ) {
        cerr<<"error: level - illegal index " << i << endl;
        exit( 1 );
    }
    if ( i>=nlev) { 
        while ( i >= size )
            Resize();
        nlev = i+1;
     }
    return timedim[i];
}

Field
Field::operator=(Field& f)
{
    delete[] timedim;
    strcpy( descrip, f.descrip );
    strcpy( this->format, format );
    nlev = f.nlev;
    ntime = f.ntime;
    group = f.group;
    id = f.id;
    size = f.nlev;
    unitConversion = f.unitConversion;
    timedim = new TimeSlice[nlev];
    for ( int i=0; i<nlev; i++ ) 
        timedim[i] = f.timedim[i];
    return *this;
}

void
Field::Resize() {
    TimeSlice* tmp = new TimeSlice[size+bumpsize];
    for ( int i=0;i<size;i++ )
        tmp[i] = timedim[i];
    delete[] timedim;
    timedim = tmp;
    size += bumpsize;
}

void
Field::SetDataIsMissing( bool isMissing )
{
    dataIsMissing = isMissing;
    if ( isMissing == true )
        SetMissingValueString();
}

void
Field::SetMissingValueString()
{
    
    missingValue[0]='-';
    for (int i=1;i<(int)sizeof(missingValue)-1;i++)
        missingValue[i]='9';
    char w[32];
    char p[32];
    int width;
    int precision;
    
    int wLen = strcspn( format, "." ) - 1;
    int pLen = strcspn( format, "dfe" ) - wLen - 2;
    
    strncpy( w, &format[1], wLen );
    w[wLen] = 0;
    strncpy( p, &format[wLen+2], pLen );
    p[pLen] = 0;
    char style = format[wLen+pLen+2];
    
    width = atoi(w);
    precision = atoi(p);

    switch( style ) {
    case 'd':
        break;
    case 'f':
        missingValue[width-precision-1]='.'; 
        break;
    case 'e':
//         missingValue[2]='.';
//         missingValue[width-4]='E';
//         missingValue[width-3]='-';
        missingValue[width-precision-1]='.'; 
        break;
    default:
        cerr << "ERROR - GetMissingValue(): can't handle format \""<<
            "\"" << endl;
    }
    missingValue[width]=0;
}

void
Field::WriteData( ofstream& o )
{
    char fData[256];            // formatted numerical string
    
    for ( int l=0; l<GetNlev(); l++ ) {
        for ( int t=0; t<GetNtime(); t++ ) {
            if ( dataIsMissing )
                o << missingValue;
            else {
                sprintf( fData, format, (*this)[l][t] * unitConversion );
                o << " " << fData;
            }
            if ( ! ((t+1)%10) && t!=GetNtime()-1 )
                o << endl;
        }
        o << endl << endl;
    }
}

void
Field::WriteData( ofstream& o, int timeIdx )
{
    char fData[256];            // formatted numerical string
    for ( int l=0; l<GetNlev(); l++ ) {
        if ( dataIsMissing )
            o << missingValue;
        else {
            sprintf( fData, format, (*this)[l][timeIdx] * unitConversion );
            o << fData;
        }
        if ( !((l+1)%10) && l!=GetNlev()-1 )
            o << endl;
    }
}

