#ifndef LIST_H
#define LIST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// this whole class could be replaced with an off-the-shelf
// list template class, but there are a couple of problems 
// with the available alternatives: if standard template library
// classes are used, then the gcc standard template library must
// be installed wherever the GUI is run (because of dynamic linking).
// Qt has a list class, which works in a very similar method, but
// if the single column model is made to work without Qt, which is
// not really hard to do (designed specifically for this) then you
// wouldn't have access to the QList type anymore.

template<class type> class ListItr;

template <class type>
class ListNode {
public:
    ListNode() { next = 0; value = 0; }
    ListNode( type* value ) { next = 0; this->value = value; }
    ListNode* next;             // pointer to next node
    type* value;                // data "contained" in this node
};


template <class type>
class List
{
      // elements are stored in a linked list
    friend class ListItr<type>;
public:
    List() {
        head = new ListNode<type>;
        head->next = 0;
        count = 0;
    }
    ~List() { 
        Clear(); 
    }

    type*   operator[]( int index ) {
        if ( index < 0 || index >= count )
            return (type*) 0;
        ListItr<type> it( this );
        for( int i=0; i<index; ++it )
            ; // do nothing
        return it.Current();
    }

    void   Insert( type* value ){
          // insert at front of linked list
        ListNode<type>* tmp = new ListNode<type>( value ); 
        tmp->next = head->next;
        head->next = tmp;
        count++;
    }

    void     Clear(){
        ListNode<type>* cursor = head->next;
        while ( cursor != 0 ) {
            ListNode<type>* tmp = cursor->next;
            if ( autoDelete == true ) 
                delete cursor->value;
            delete cursor;
            cursor = tmp;
        }
        head->next = 0;
        count = 0;
    }

    int      Count() { return count; }
    void     Remove( type* t ) {     
        ListNode<type>* cursor = head;    
        ListNode<type>* tmp = head;    
        while ( (cursor = cursor->next) ) {
            if ( cursor->value == t  ) {
                if ( autoDelete == true )
                    delete cursor->value;
                tmp->next = cursor->next;
                delete cursor;
                count--;
                return;
            }
            tmp = cursor;
        }
    }

    void     SetAutoDelete( bool enable ) { autoDelete = enable; }

protected:
    ListNode<type>* head;       // head of linked list
    int         count;          // number of elements contained in list
    bool        autoDelete;     // delete elements when removed or list is destroyed
};

template <class type>
class ListItr
{
public:
    ListItr( const List<type>* l ) { list = l; cursor = l->head->next; }
    type* Current() const { if( cursor != 0 ) return cursor->value; else return 0; }
    void  ToFirst() { cursor = list->head->next; }
    void  operator++() { cursor = cursor->next; }
private:
    const List<type>* list;
    ListNode<type>* cursor;
};

#endif // LIST_H



