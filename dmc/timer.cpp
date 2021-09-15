#include "timer.h"

void timer::start(){_start = std::chrono::high_resolution_clock::now(); }
void timer::stop(){_stop =  std::chrono::high_resolution_clock::now();addElapsedTime(); }

void timer::addElapsedTime(){duration+=_stop- _start;}

timer::timer() : duration(0){}


void timers::add(std::string name) {_timers[name]=timer();}
timer & timers::get(const std::string & name) {return _timers.at(name);}

timers*  timers::get()    {
  if (!singleton)
    singleton = new timers();
  return singleton;
}

timers* timers::singleton = nullptr;

timers&  getTimers() {return *(timers::get()) ;}
