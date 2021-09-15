#ifndef TIMERS_H
#define TIMERS_H


#include <string>
#include <memory>
#include <map>
#include "omp.h"

class timer
{
    public:

    virtual void start()=0;
    virtual void stop()=0;

    virtual double timeElapsed()=0;

};


struct openMPClock
{
    static double time() {return omp_get_wtime() ;}
};



template<class clock_t=openMPClock>
class clockTimer : public timer
{
    public:
    clockTimer() :
    tStart(0),
    tStop(0),
    tElapsed(0) {}

    virtual void start()
    {
        tStart=clock_t::time();
    }

    virtual void stop()
    {
        tStop=clock_t::time();
        tElapsed+=tStop-tStart;
    }

    virtual double timeElapsed() { return tElapsed;}

    private:

    double tElapsed;
    double tStart;
    double tStop;
};




class timers
{
    public:

       static timers& getInstance( );



    void add( std::string label , std::shared_ptr<timer>  timerObj ) {   
        _timers[label]=timerObj;
      }


    void add(std::string label)
    {
        add(label,std::make_shared<clockTimer<openMPClock> >() );
    }

    timer & operator[] (const std::string & label) {return *(_timers.at(label));}

    std::string report() const;

    void start(std::string label)
    {
        if (_timers.count( label ) == 0 )
        {
            add(label);
        }
        (*this)[label].start();
    }
    
    private:

    timers(){}

    static timers* _singleton;

    std::map<std::string,std::shared_ptr<timer>  > _timers;

};


#define START_TIMER(name) timers::getInstance().start(name)
#define STOP_TIMER(name) timers::getInstance()[name].stop()
#define ADD_TIMER(name) timers::getInstance().add(name)




#endif