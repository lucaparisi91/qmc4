#ifndef NCONNECTEDCHAINS_H
#define NCONNECTEDCHAINS_H

#include "pimcConfigurations.h"

namespace pimc
{
    using configurations_t = pimcConfigurations ;

    class nConnectedChains
    {
        public:
        nConnectedChains() {}
        size_t count(const configurations_t & configurations, int set);

        private:
        std::vector<bool> _visited;
    };


};

#endif