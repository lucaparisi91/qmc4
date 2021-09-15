#include "timers.h"
#include <iostream>
#include <sstream> 

timers & timers::getInstance()
{
    if (_singleton == nullptr )
    {
        _singleton= new timers();
    }

    return (*_singleton);
}

timers*  timers::_singleton = nullptr;




std::string timers::report() const
{
    std::stringstream ss;

    for (auto it : _timers)

    {
        ss << "Timer " <<it.first <<": " << (it.second)->timeElapsed() << std::endl; 
    }
    
    return ss.str();
}