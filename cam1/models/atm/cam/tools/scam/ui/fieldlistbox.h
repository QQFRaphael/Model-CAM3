#ifndef FieldListBox_H
#define FieldListBox_H

#include <string>
#include <qlistbox.h>
#include <qkeycode.h>
#include <qpainter.h>
#include <qpen.h>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

class FieldListBoxText: public QListBoxText
{
public:
    FieldListBoxText( string text, bool modifiable ) : QListBoxText( text.c_str() ), _modifiable( modifiable ) {}
    void FieldListBoxText::paint( QPainter* p ) {
        if ( _modifiable == true ) {
            QPen pen( Qt::red, 2 );
            p->setPen( pen );
        }
        QListBoxText::paint( p );
    }
private:
    bool _modifiable;
};


class FieldListBox : public QListBox
{

public:
    FieldListBox( QWidget* parent, const char* name );
    string FieldName( int index );
protected:
    void ClearAll();            // deselects all entries
    void SelectAll();           // selects all entries
    void Reset();               // resets the search
      // Overridden QListBox methods
    void keyPressEvent( QKeyEvent* keyEvent );
    void focusOutEvent( QFocusEvent* ) { } // override base class implementation

private:
    string             searchStr;
};

#endif  // FieldListBox_H
