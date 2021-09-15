#include "factory.h"

factory&  getFactory() {return *(factory::get()) ;}


factory*  factory::get()    {
  if (!singleton)
    singleton = new factory;
  return singleton;
  
}

factory* factory::singleton = nullptr;
