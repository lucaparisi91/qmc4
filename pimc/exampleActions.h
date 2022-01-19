#ifndef EXAMPLEACTION_H
#define EXAMPLEACTION_H

#include "action.h"
#include "pimcConfigurations.h"
#include "propagators.h"
#include "pairProductKernel.h"
#include "actionTwoBody.h"

namespace pimc
{
    pimc::firstOrderAction createFreeCaoBerneAction(Real radius, Real timeStep , const geometry_t & geo,range_t groups, int nChains,int M);
};

#endif