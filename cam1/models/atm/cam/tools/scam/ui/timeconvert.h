/*------------------------------------------------------------------------*
 * File: timeconvert.h 
 * $Author: jet $
 * $Id: timeconvert.h,v 1.1.6.1 2004/05/18 20:52:03 jet Exp $ *
 *------------------------------------------------------------------------*/
#ifndef TIMECONVERT_H
#define TIMECONVERT_H

#include <string>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

class TimeConverter
{
public:
    enum Time { STEPS, HOURS, DAYS, DATE };
    TimeConverter();
    TimeConverter( const string& timeString );
    TimeConverter( int steps );
    TimeConverter( int basedate, int steps );
    void SetTime( const string& timestring, Time format );
    void SetTime( const string& timestring );
    void SetStep( int steps );
    int Steps();
    real_t Days();
    int Date();
    real_t Hours();
    string TimeToString(); // return time in current time format
    string TimeToString( Time format ); 
    
    static void SecondsToDate( int secs, int baseDate, int& outDate, int& outSecs );
    static int DaysToDate( int bdate, int days );
    static int DateToJulianDay(int bdate);
    static int DaysBetweenDates(int date1, int date2);
    
private:
    int steps;
    int baseDate;
    static char daysStr[256];
    static char hoursStr[256];
    static char stepsStr[256];
    static char dateStr[256];

};

// Time enumeration was made part of the TimeConverter class
//  to avoid potential name clashes with the enumerated types
typedef TimeConverter::Time Time;



#endif // TIMECONVERT_H



