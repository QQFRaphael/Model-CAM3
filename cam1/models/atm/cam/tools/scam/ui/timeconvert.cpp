#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "manager.h"
#include "timeconvert.h"

#define YEAR(yymmdd)	 ((yymmdd) / 10000)
#define MONTH(yymmdd)	(((yymmdd) % 10000) / 100)
#define DAY(yymmdd)	((yymmdd) % 100)
#define SEC_PER_DAY     (86400)
#define SEC_PER_YEAR    (SEC_PER_DAY * 365)


TimeConverter::TimeConverter()
{
    baseDate = MANAGER.BaseDate(); 
}

TimeConverter::TimeConverter( int steps )
{
    baseDate = MANAGER.BaseDate(); 
    this->steps = steps;
}

TimeConverter::TimeConverter( int baseDate, int steps )
{

    this->baseDate = baseDate;
    this->steps = steps;
}

TimeConverter::TimeConverter( const string& timeString )
{
    baseDate = MANAGER.BaseDate(); 
    SetTime( timeString );
}

void 
TimeConverter::SetStep( int steps )
{
    this->steps = steps;
}

void 
TimeConverter::SetTime( const string& timestring, Time format )
{
    real_t hours, days;
    int date;

    switch( format ) {
    case STEPS:
        steps = atoi( timestring.c_str() );
        break;
    case HOURS:
        hours = atof( timestring.c_str() );
        steps = (int) (( hours * 3600.0) / MANAGER.StepLen());
        break;
    case DAYS:
        days =  atof( timestring.c_str() );
        steps = (int)( (days * 86400.0) / MANAGER.StepLen());
        break;
   case DATE:
       date =  atoi( timestring.c_str() );
       steps = (int)((( DaysBetweenDates( baseDate, date ) *
                        86400.0 )) / MANAGER.StepLen() );
       break;
    default:
        cerr << "ERROR: " __FILE__ ":" << __LINE__
             << " TimeConverter::SetTime():  bad time format:" << (int)format << endl;
        exit( -1 );
    }
}

void 
TimeConverter::SetTime( const string& timestring )
{
    real_t hours, days;
    int date;

    switch( MANAGER.TimeFormat() ) {
    case STEPS:
        steps = atoi( timestring.c_str() );
        break;
    case HOURS:
        hours = atof( timestring.c_str() );
        steps = (int) (( hours * 3600.0) / MANAGER.StepLen());
        break;
    case DAYS:
        days =  atof( timestring.c_str() );
        steps = (int)( (days * 86400.0) / MANAGER.StepLen());
        break;
   case DATE:
       date =  atoi( timestring.c_str() );
       steps = (int)((( DaysBetweenDates( baseDate, date ) *
                        86400.0 )) / MANAGER.StepLen() );
       break;
    default:
        cerr << "ERROR: "__FILE__":"<< __LINE__ 
             << " TimeConverter::SetTime():  bad time format: " << int(MANAGER.TimeFormat()) << endl;
        exit( -1 );
    }
}

int TimeConverter::Steps()
{
    return steps;
}

real_t
TimeConverter::Days()
{
    return steps * MANAGER.StepLen() / 3600.0;
}

int
TimeConverter::Date()
{
    int secs, date;
    SecondsToDate( baseDate, steps*MANAGER.StepLen(), date, secs );
    return date;
}

real_t 
TimeConverter::Hours()
{
    return steps * MANAGER.StepLen() / 3600.0;
}


string
TimeConverter::TimeToString( Time format )
{
    char timeString[256];
    real_t days, hours;
    int date;

    switch( format ) {
    case STEPS:
        sprintf( timeString,"%d", steps );
        break;
    case HOURS:
        hours  = steps * MANAGER.StepLen() / 3600.0;
        sprintf( timeString,"%.1f", hours );
        break;
    case DAYS:
        days = steps * MANAGER.StepLen() /  86400.0;
        sprintf( timeString,"%.1f", days );
        break;
   case DATE:
       date = DaysToDate( baseDate,(int)rint( steps* MANAGER.StepLen()/86400.0) ); 
       sprintf( timeString,"%d", date );
       break;        
    default:
        cerr <<  "ERROR: "__FILE__", line " << __LINE__ 
             << " TimeConverter::TimeToString(): - bad time format: "<< int(format) << endl;
        exit( -1 );
    }
    return timeString;
}

string
TimeConverter::TimeToString()
{
    char timeString[256];
    real_t days, hours;
    int date;

    switch( MANAGER.TimeFormat() ) {
    case STEPS:
        sprintf( timeString,"%d", steps );
        break;
    case HOURS:
        hours  = steps * MANAGER.StepLen() / 3600.0;
        sprintf( timeString,"%.1f", hours );
        break;
    case DAYS:
        days = steps * MANAGER.StepLen() /  86400.0;
        sprintf( timeString,"%.1f", days );
        break;
   case DATE:
       date = DaysToDate( baseDate,(int)rint(steps* MANAGER.StepLen()/86400.0) ); 
       sprintf( timeString,"%d", date );
       break;        
    default:
        cerr <<  "ERROR: "__FILE__", line " << __LINE__ 
             << ": TimeConverter::TimeToString(): - bad time format: "<< int(MANAGER.TimeFormat()) << endl;
    }
    return timeString;
}

int 
TimeConverter::DaysToDate( int bdate, int days )
{
    int date;
    int secs;

    SecondsToDate( days * 86400, bdate, date, secs );
    return date;
}

int 
TimeConverter::DateToJulianDay(int bdate)
{
    int julianDay[] = {0, 31,59,90,120,151,181,212,243,273,304,334 };

    int bday, jday;
    int bmonth;

    if ( bdate == 0 )
        return 0;
    
    bday = bdate%100;
    bmonth = (bdate%10000)/100;
        
    jday = julianDay[bmonth-1] + bday;
    return jday;
}

int
TimeConverter::DaysBetweenDates(int date1, int date2) 
{
    int jday1 = DateToJulianDay( date1 );
    int jday2 = DateToJulianDay( date2 );
    int year1 = (date1%1000000)/10000;
    int year2 = (date2%1000000)/10000;
    int days = (year2 - year1) * 365 + jday2 - jday1;
    return days < 0 ? -days : days;
}

void
TimeConverter::SecondsToDate( int seconds, int basedate ,
                              int& outDate, int& outSecs )
{
    int jdcon[] = {0,0,31,59,90,120,151,181,212,243,273,304,334, 365};
    int YY,MM,DD,SS;
    int jday, jsec;            // julian day, julian sec
    
    // Convert the base date to a julian day 
    jday = jdcon[MONTH( basedate )] + DAY( basedate );
    // find the Julian seconds and add the offset of seconds.
    jsec = ( jday-1 ) * SEC_PER_DAY + seconds;
    YY = YEAR( basedate ) + ( jsec / SEC_PER_YEAR );
    jsec = jsec % SEC_PER_YEAR;
    for ( MM = 0; ( jsec - ( jdcon[MM] ) * SEC_PER_DAY ) >= 0 ; MM++ )
        ; //  Do Nothing 
    --MM;
    jsec = jsec - jdcon[MM] * SEC_PER_DAY;
    DD = jsec / SEC_PER_DAY +1;
    jsec = jsec % SEC_PER_DAY;
    SS = jsec;
    outDate =  ( YY * 10000 )+ ( MM * 100 ) + DD;
    outSecs = SS;
}







