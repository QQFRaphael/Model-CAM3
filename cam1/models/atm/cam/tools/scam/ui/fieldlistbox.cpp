#include <iostream>
#include <string>
#include <ctype.h>

#include <qapp.h>
#include <qlistbox.h>
#include <qwidget.h>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

#include "fieldlistbox.h"
#include "max.h"

FieldListBox::FieldListBox( QWidget* parent, const char* name ) : QListBox( parent, name )
{
      // respond to keyboard events
    setFocusPolicy( QWidget::StrongFocus );
}

string
FieldListBox::FieldName( int index )
{
      // truncate the name after the short name (look for space)
    string textField = text(index).ascii();   // text is a member of QListBox 
    return textField.substr( 0, textField.find_first_of( " " ) );
}

void
FieldListBox::keyPressEvent( QKeyEvent* keyEvent )
{
    switch ( keyEvent->key() ) {
          // allow base class listbox to receive signals that
          // it regularly handles.
          // (e.g., so that it will emit "selected" signal when 
          //  receiving returns)
    case Key_Return:
        if ( isMultiSelection() ) // toggle the selection
            setSelected( currentItem(), !isSelected( currentItem() ) );
        Reset();                // presumably the desired item was found
    case Key_Down:
    case Key_Up:
    case Key_PageUp:
    case Key_PageDown:
        QListBox::keyPressEvent(keyEvent);
        return;
    case Key_Backspace: 
    case Key_Delete:
        if ( searchStr.size() > 0 ) 
            searchStr.erase( searchStr.end() - 1 );
        else
            QApplication::beep();
        break;
    case Key_Escape:
        searchStr = "";
        break;
    case Key_A:
        if ( keyEvent->state() == AltButton ) {
            SelectAll();
            break;
        }
    case Key_C:
        if ( keyEvent->state() == AltButton ) {
            ClearAll();
            break;
        }
    default:
        int asc = keyEvent->ascii();
        if ( !isalnum( asc )) {
            keyEvent->ignore();
            return;
        }
          // convert search string to uppercase
        asc = toupper( asc );
        searchStr += asc;
    }
    
    size_t matchIdx = 0;

    while( matchIdx < count() ) {
          // look at next field in list
        int comp = searchStr.compare( FieldName( matchIdx));
	//        int comp = searchStr.compare( FieldName( matchIdx), 0, searchStr.size() );
	//        int comp = searchStr.compare( 0, searchStr.size(), FieldName( matchIdx) );
        if( comp > 0 ) 
            matchIdx++;
        else if ( comp == 0 )
            break;              // found it
        else {
            QApplication::beep();
            searchStr.erase( searchStr.end()-1 );
            break;
        }
    }
      // check if we've gone past the end of the list
    if ( matchIdx == count() ) { 
        QApplication::beep();
        searchStr.erase( searchStr.end()-1 );
        matchIdx--; // match the end of the list
    }
    setCurrentItem( matchIdx );
    if ( ! itemVisible( matchIdx ) )
         centerCurrentItem();
}

void
FieldListBox::Reset()
{
    searchStr = "";
}

void 
FieldListBox::SelectAll()
{
    for ( int i=0; i< (int)count(); i++ ) 
        setSelected( i, TRUE );
    Reset();
}

void
FieldListBox::ClearAll()
{
    for ( int i=0; i< (int)count(); i++ ) 
        setSelected( i, FALSE );
    Reset();
}
