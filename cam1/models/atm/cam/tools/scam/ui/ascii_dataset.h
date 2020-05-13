/*------------------------------------------------------------------------*
 * File: ascii_dataset.h 
 * $Author: jet $
 * $Id: ascii_dataset.h,v 1.1.6.1 2004/05/18 20:51:50 jet Exp $ *
 *------------------------------------------------------------------------*/
#ifndef ASCII_DATASET_H
#define ASCII_DATASET_H

#include "realtype.h"
#include "ncfile.h"

class Ascii_Dataset 
{

public:
    Ascii_Dataset( string filename, int type )
        : _name(filename),_type(type) {}
    Ascii_Dataset::Ascii_Dataset( int type )
        : _type(type) {}
    const string& name() const {return _name;}
    int type(){return _type;}
    
protected:
    string   _name;
    int   _type;
};

#endif // ASCII_DATASET_H



