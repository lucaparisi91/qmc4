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






class nullPotentialActionConstructor : public actionConstructor
{
    public:
    nullPotentialActionConstructor() {}

    void setGeometry( const geometry_t & geo ) { _geo=geo;}

    void setTimeStep( Real timeStep) { _timeStep = timeStep; }



    std::shared_ptr<action> create( const json_t & j )
    {
        auto S = std::make_shared<nullPotentialAction>(_timeStep,_geo);
        return S;
    }

    private:

    geometry_t _geo;
    Real _timeStep;
};



};

#endif