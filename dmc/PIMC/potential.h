#include "../tools.h"
#include "../traits.h"

namespace pimc
{

    template<class F_t>
    class potential
    {

        Real operator()( const states_t & states )=0;
        Real operator()(const states_t & states)=0;
        Real operator()(const distances_t & distances )=0;
        Real operator()(const states_t & states,const distances_t & distances, int iModified, int iStart, int iEnd )=0; // interaction potential change when one particle is modifed due to interactions with (iStart,iEnd) others

    private:
        geometry_t & geo;
        F_t V;  
    };

}