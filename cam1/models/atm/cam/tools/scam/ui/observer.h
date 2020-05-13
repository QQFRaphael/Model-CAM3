#ifndef OBSERVER_H
#define OBSERVER_H

#include <list>

#ifndef _RWSTD_NO_NAMESPACE
  using namespace std;
#endif

class Observer {
public:
    virtual void Update() {};
};

class Observable {
    typedef list<Observer*>		    ObsList;
    typedef list<Observer*>::iterator       ObsListItr;
    typedef list<Observer*>::const_iterator ObsListConstItr;
    
public:
             Observable() { isDirty = false; }
    virtual ~Observable() { observers.clear(); }
    void     AddObserver( Observer* o ) { observers.push_back( o ); }
    void     ClearChanged() { isDirty = false; }
    void     DeleteObserver( Observer* o ) { observers.remove( o ); }
    void     NotifyObservers() ;
    void     DeleteObservers() { observers.clear(); }
    bool     HasChanged() const { return isDirty; }
    int      NumObservers() const { return observers.size(); }
    void     SetChanged() { isDirty = true; }

protected:
    bool     isDirty;
    ObsList  observers;
};


inline void Observable::NotifyObservers() {
    ObsListItr it = observers.begin();
    for ( ;it != observers.end(); ++it )
        (*it)->Update();
}

#endif

