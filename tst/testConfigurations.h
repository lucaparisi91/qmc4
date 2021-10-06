#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "../pimc/action.h"
#include "../pimc/pimcConfigurations.h"
#include "../pimc/moves.h"
#include "../pimc/action.h"

#define TOL 1e-7

class configurationsTest : public ::testing::Test {
protected:
    configurationsTest() {
        open=NULL;
        close=NULL;
    }

    void SetRandom(const std::array<Real,DIMENSIONS> & lBox={ TRUNCATE_D(1,1,1) });

    void SetUp( int N_ , int M_ , Real Beta_,std::array<Real,getDimensions()> lBox={ TRUNCATE_D(3000,3000,3000)}) ;

    void SetUp( std::vector<int> Ns , int M_ , Real Beta_,std::array<Real,getDimensions()> lBox={ TRUNCATE_D(3000,3000,3000)});
    

    void SetSeed(int seed)
    {
        randG=randomGenerator_t(seed);
    }

    void SetUpFreeParticleAction();

    void SetUpTwoBodyInteractionHarmonic();


    void SetUpTwoBodyInteractionGaussian( Real V0 , Real alpha );
    
    void SetUpTwoBodyInteractionHarmonicInTrap();


    void SetUpNonInteractingHarmonicAction(int order=2);
    

    void SetGrandCanonicalEnsamble(Real chemicalPotential)
    {
        configurations.setEnsamble(pimc::ensamble_t::grandCanonical);
        configurations.setChemicalPotential({chemicalPotential});
    }

    void SetSemiGrandCanonicalEnsamble(Real deltamu)
    {
        std::vector<Real> deltamus{deltamu};

        configurations.setDifferenceChemicalPotential(deltamus);
    }


   

    void addOpenCloseFixedSegment(Real C , int l , int t0)
    {
        open=new pimc::openMove(C, 0, l );
        close=new pimc::closeMove(C, 0, l );

        open->setStartingBead(t0);
        close->setStartingBead(t0);

        open->setLengthCut(l);
        close->setLengthCut(l);

        tab.push_back(open,1,pimc::sector_t::diagonal);
        tab.push_back(close, 1, pimc::sector_t::offDiagonal);
    }



    
    template<class T>
    void accumulate(int nBurns,int nTrials,const T & f, int correlationSteps=1, pimc::sector_t sector= pimc::sector_t::diagonal)
    {
        for(int k=0;k<nBurns;k++)
        {
            bool accept= tab.attemptMove(configurations,S,randG);    
        }

        for(int k=0;k<nTrials;k++)
        {
            for (int kk=0;kk<correlationSteps;kk++)
            {
                 bool accept= tab.attemptMove(configurations,S,randG);

            }


            nMoves+=1;
            if ( configurations.isOpen() )
            {
                    

                    nOpen+=1;

                    if (sector == pimc::sector_t::offDiagonal)
                    {
                        f(configurations,S);
                    }
                    
                
            }
            else
            {
                nClosed+=1;

                  if (sector == pimc::sector_t::diagonal)
                    {
                        f(configurations,S);
                    }


            }

           

        }


    }



    virtual void TearDown() {
        if (open != NULL )
        {
            delete open;
        }

        if (close != NULL )
        {
            delete close;
        }

    }

    void resetCounters()
    {
        nMoves=0;
        nOpen=0;
        nClosed=0;
    }


    int N;
    int M;
    Real Beta;
    Real timeStep;
    pimc::pimcConfigurations configurations;
    pimc::firstOrderAction S;
    randomGenerator_t randG;
    pimc::geometryPBC_PIMC geo;
    pimc::tableMoves tab;

    pimc::openMove *open;
    pimc::closeMove *close;

    Real nOpen=0;
    Real nClosed=0;
    Real nMoves=0;

};


Real accumulateLengthSquare(const pimc::configurations_t & configurations, const std::array<int,2> & particleRange, std::array<int,2> timeRange,const geometryPBC & geo);



