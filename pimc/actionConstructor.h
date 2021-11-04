#include "action.h"
#ifndef ACTION_CREATOR_H
#define ACTION_CREATOR_H

namespace pimc{

    struct actionConstructor
    {
        virtual std::shared_ptr<action> create(const json_t & j) { 
            throw std::runtime_error("Creation of action not implemented");
             return std::make_shared<nullPotentialAction>();};

    };


};

#endif