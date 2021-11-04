#include "moves.h"
#include "action.h"
#include "qmcExceptions.h"
#include "toolsPimcTest.h"

namespace pimc
{      
    std::array<int,2> timeSliceGenerator::operator()(randomGenerator_t & randG, int nBeads , int maxBeadLength)
    {
        int t0=std::floor( uniformRealNumber(randG)*nBeads );
        int length = uniformRealNumber(randG)*maxBeadLength;

        int t1 = t0 + length;

        return {t0,t1};
    }
    
levyMove::levyMove(int maxBeadLength_,int set) : _levy(maxBeadLength_) , uniformRealNumber(0,1),maxBeadLength(maxBeadLength_) , buffer(( maxBeadLength_+1)*2,getDimensions() ) ,singleSetMove::singleSetMove(set){}

bool levyMove::attemptMove( configurations_t & confs, firstOrderAction & ST,randomGenerator_t & randG)
{

    int nChains = confs.nChains();
    int nBeads = confs.nBeads();
    auto & S = ST.getPotentialAction();

    auto & geo = S.getGeometry();

    int iChain=confsSampler.sampleChain(confs,getSet(),randG);

    if (iChain < 0)
    {
        return false;
    }


    auto timeRange = tGen(randG,nBeads,maxBeadLength);

    int M= confs.nBeads();


    std::array<int,2> timeRange1 { timeRange[0] , std::min(timeRange[1],M) - 1 } ;
    std::array<int,2> timeRange2 { 0  , timeRange[1] - M -1  };


    const auto & currentChain = confs.getChain(iChain);

    if (timeRange[0] <= currentChain.tail)
    {
        return false;
    }
    int iChainNext = currentChain.next;

    if (
        (timeRange[1] > currentChain.head and  (currentChain.next == -1 ) )
        or (   ( currentChain.next != -1 ) and
            ( confs.getChain(iChainNext).head < timeRange2[1] + 1   )
    )
        )
    {
        return false; // do not accept if the time range is not included inside a chain (including time PBC)
    }
    
    auto & data = confs.dataTensor();

    
    auto sOld = S.evaluate(confs, timeRange1, iChain);
    sOld += S.evaluate(confs, timeRange2, iChainNext);


    // copy to internal buffer beads to move
     confs.copyDataToBuffer(buffer,{timeRange1[0] , timeRange1[1]+ 1},iChain);

     if (iChainNext!=-1)
     {
        confs.copyDataToBuffer(buffer,{0,timeRange2[1]+1},iChainNext, timeRange1[1]  - timeRange1[0] + 2);
     }


    if (timeRange[1] > confs.nBeads() )//copy the end bead in second chain to the first chain
    {
        confs.copyData( { timeRange2[1] + 1 , timeRange2[1] + 1 } , iChainNext, timeRange[1],iChain );

        for(int d=0;d<getDimensions();d++)
        {
            data(iChain,d,timeRange[1])-= data(iChainNext,d,0) - data(iChain,d,confs.nBeads() );
        }// ensures to reconstruct along a continuos path
    }
    
    _levy.apply(confs,timeRange,iChain,S,randG);
    
    if (iChainNext!=-1)
     {
         confs.copyData( { timeRange1[1]+1 , timeRange[1]-1 } , iChain, 0 ,iChainNext ); // time periodic boundary conditions

        for(int t=0;t<=timeRange2[1];t++)
            for(int d=0;d<getDimensions();d++)
            {
                data(iChainNext,d,t)+=
                data(iChainNext,d,timeRange2[1] + 1 ) - data(iChain,d,timeRange[1] );

                //std::cout << data(iChainNext,d,timeRanges[1][1]) - data(iChain,d,timeRange[1] )<<std::endl;
            }// ensures to reconstruct along a continuos path
     }

    bool accepted = S.checkConstraints(confs,timeRange1,{iChain,iChain});
    accepted = accepted and S.checkConstraints(confs,timeRange2,{iChainNext,iChainNext});
    
    if (accepted)
    {
        auto  sNew= S.evaluate(confs,timeRange1, iChain) ;
        sNew+=S.evaluate(confs,timeRange2, iChainNext) ;
        const auto actionDifference = sNew - sOld;

        accepted = sampler.acceptLog(-actionDifference,randG);

    }

    if (! accepted)
    {
        // copy back old beads

        confs.copyDataFromBuffer(buffer,{timeRange1[0] , timeRange1[1]+ 1},iChain);

        if (iChainNext!=-1)
         {
            confs.copyDataFromBuffer(buffer,timeRange2,iChainNext, timeRange1[1]  - timeRange1[0] + 2);
         }

    }
    return accepted;
}

swapMove::swapMove(int maxStepLength_,int maxN,int set) :
 maxStepLength(maxStepLength_),buffer(maxStepLength_*2,getDimensions() )
 , uniformRealNumber(0,1) ,_levy(maxStepLength_+2),particleSampler(maxN)
, singleSetMove(set){
    setRandomLength();
}


bool swapMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    if ( confs.getEnsamble() == ensamble_t::canonical)
    {
        return attemptCanonicalMove(confs,S,randG);
    }
    else if ( confs.getEnsamble() == ensamble_t::grandCanonical)
    {
        return attemptGrandCanonicalMove(confs,S,randG);
    }

    return false;
}


bool swapMove::attemptCanonicalMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG)
{
    if ( ! confs.isOpen(getSet()) )
    {
        throw invalidState("Swap move can only be done in the G sector.");
    }
     auto & data = confs.dataTensor();
    // selects an head at random

    const auto & group = confs.getGroups()[getSet()];

    const auto & heads =  confs.getGroups()[getSet()].heads;
    int iiChainHead =  std::floor(  uniformRealNumber( randG )*heads.size() ) ;
    int iChainHead = heads[iiChainHead];

    auto & Spot = S.getPotentialAction();

    // selects a time slice length at random. Refuse if reconstructed time slice crosses over the end bead

    int l = std::floor(uniformRealNumber(randG) * maxStepLength) + 1;
    //int l = maxStepLength;

    
    int iHead=confs.nBeads();

    
    // tower sampling a particle i with gaussian weights on the relative distances
   
    
    auto & geo = S.getGeometry();

    std::array<Real, getDimensions()> distance;
    particleSampler.reset();

    Real weightForwardMove = 0;
    
    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        
        for(int d=0;d<getDimensions();d++)
        {
            distance[d]=geo.difference(  data(iChainHead,d,iHead) - data(i,d,0),d) - data(i,d,l) + data(i,d,0)    ;
        }


        particleSelectionWeight=exp(freeParticleLogProbability(distance,S.getTimeStep()*l,group.mass));
        weightForwardMove+=particleSelectionWeight;
        particleSampler.accumulateWeight(particleSelectionWeight);
    }

    int iPartner=particleSampler.sample(randG) + group.iStart;

    if ( confs.getChain(iPartner).hasTail() )
    {
        return false;
    }

    // metropolic test based on ratio of forward and backward move
    Real weightBackwardMove= 0;

    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        Real norm=0;
        for(int d=0;d<getDimensions();d++)
        {
            distance[d]=geo.difference(  data(iPartner,d,0) - data(i,d,0),d) + data(i,d,0) - data(i,d,l);
        }


        weightBackwardMove+=exp(freeParticleLogProbability(distance,S.getTimeStep()*l,group.mass));
    }


    

    Real deltaS=-(log(weightForwardMove) - log(weightBackwardMove)) ;

    
    const auto  partnerChain = confs.getChain(iPartner);

    
    deltaS-=Spot.evaluate(confs,{0,l-1},iPartner);

    confs.copyDataToBuffer(buffer,{0,l-1},iPartner);


    // performs levy reconstruction between the head and the bead
    for(int d=0;d<getDimensions();d++)
        {
            data(iPartner,d,0)=//data(iChainHead,d,iHead);
            geo.difference (- data (iPartner,d,0) + data(iChainHead,d,iHead) , d) -  data(iPartner,d,l) + data(iPartner,d,0);
            data(iPartner,d,0)+=data (iPartner,d,l);
        }

    _levy.apply(confs,  {0,l },  iPartner, S ,randG );
    deltaS+=Spot.evaluate(confs,{0,l-1},iPartner);
    
    bool accept = metropolisSampler.acceptLog(-deltaS,randG);
    if (accept) 
    {
        confs.join(iChainHead,iPartner );
    }
    else
    {
        confs.copyDataFromBuffer(buffer,{0,l-1},iPartner);
    }    

    return accept;
}







bool swapMove::attemptTestMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG )
{
       if ( ! confs.isOpen(getSet()) )
    {
        throw invalidState("Swap move can only be done in the G sector.");
    }

     auto & data = confs.dataTensor();
    // selects an head at random

    const auto & group = confs.getGroups()[getSet()];

    const auto & heads =  confs.getGroups()[getSet()].heads;

    
    int iiChainHead =  std::floor(  uniformRealNumber( randG )*heads.size() ) ;
    int iChainHead = heads[iiChainHead];

    auto & Spot = S.getPotentialAction();

    // selects a time slice length at random. Refuse if reconstructed time slice crosses over the end bead


    int l = sampleLength();

    //int l = maxStepLength;
    int M= confs.nBeads();
    
    int tHead=confs.getChain(iChainHead).head;
    int tTail=confs.getChain(iChainHead).tail;

    int tJoin = (tHead + l) %M;
    
    // tower sampling a particle i with gaussian weights on the relative distances
    
    auto & geo = S.getGeometry();

    const auto & tags = confs.getTags();

    int iPartner=2;
    std::array<Real,getDimensions()> distance;

    for(int d=0;d<getDimensions();d++)
        {
            distance[d]=  data(iChainHead,d,tHead) - data(iPartner,d,tJoin)   ;
        }
    
    int iNewChainHead = iChainHead == 1 ? 2 : 1;


    Real weightForwardMove=exp(freeParticleLogProbability(distance,S.getTimeStep()*l,1)) ;



    for(int d=0;d<getDimensions();d++)
        {
            distance[d]=  data(iNewChainHead,d,tHead) - data(iPartner,d,tJoin) ;
        }

    Real weightBackwardMove=exp(freeParticleLogProbability(distance,S.getTimeStep()*l,1)) ;




    Real deltaS=-(log(weightForwardMove) - log(weightBackwardMove)) ;

    deltaS-=Spot.evaluate(confs,{tHead,std::min(tHead + l ,M) -1},iNewChainHead);
    deltaS-=Spot.evaluate(confs,{0,tHead + l - M  - 1},iPartner);

    confs.copyDataToBuffer(buffer,{0,tHead + l - M  },iPartner,0);

     for(int d=0;d<getDimensions();d++)
        {
            data(iChainHead,d,tHead + l)=data(iPartner,d,tJoin);
        }
    _levy.apply(confs,  {tHead,tHead + l },  iChainHead, S ,randG );
    confs.copyData( { M, tHead + l} , iChainHead, 0,iPartner );
    
    confs.setHead( iChainHead , std::min(tHead + l,M) );
    confs.setHead(iNewChainHead,tHead);

    deltaS+=Spot.evaluate(confs,{tHead,std::min(tHead + l ,M) -1},iChainHead);
    deltaS+=Spot.evaluate(confs,{0,tHead + l - M  - 1},iPartner);

     bool accept = metropolisSampler.acceptLog(-deltaS,randG);

     if (accept)
     {
        confs.join(iChainHead,iPartner );

     }
     else
     {
        confs.copyDataFromBuffer(buffer,{0,tHead + l - M  },iPartner,0);
        confs.setHead( iChainHead , tHead );
        confs.setHead( iNewChainHead , M );
        confs.join(iNewChainHead,iPartner);

     }

    return accept; 

}


bool swapMove::attemptGrandCanonicalMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG)
{
    if ( ! confs.isOpen(getSet()) )
    {
        throw invalidState("Swap move can only be done in the G sector.");
    }

     auto & data = confs.dataTensor();
    // selects an head at random

    const auto & group = confs.getGroups()[getSet()];

    const auto & heads =  confs.getGroups()[getSet()].heads;

    
    int iiChainHead =  std::floor(  uniformRealNumber( randG )*heads.size() ) ;
    int iChainHead = heads[iiChainHead];



    auto & Spot = S.getPotentialAction();

    // selects a time slice length at random. Refuse if reconstructed time slice crosses over the end bead


    int l = sampleLength();

    //int l = maxStepLength;
    int M= confs.nBeads();

    
    int tHead=confs.getChain(iChainHead).head;
    int tTail=confs.getChain(iChainHead).tail;

    
    int tJoin = (tHead + l) %M;
    
    // tower sampling a particle i with gaussian weights on the relative distances
    
    auto & geo = S.getGeometry();

    const auto & tags = confs.getTags();


    std::array<Real, getDimensions()> distance;
    particleSampler.reset();

    Real weightForwardMove = 0;
    
    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        
        for(int d=0;d<getDimensions();d++)
        {
            distance[d]=geo.difference(  data(iChainHead,d,tHead) - data(i,d,tJoin),d)    ;
        }

        particleSelectionWeight=exp(freeParticleLogProbability(distance,S.getTimeStep()*l,group.mass))*tags(i,tJoin) ;
        weightForwardMove+=particleSelectionWeight;

        
        particleSampler.accumulateWeight(particleSelectionWeight);

    }

    if (std::abs(weightForwardMove) < 1e-20)
    {
        return false;
    }

    int iPartner=particleSampler.sample(randG) + group.iStart;
    #ifndef NDEBUG
    assert(tags(iPartner,tJoin) != 0  );
    assert( ( iPartner >= group.iStart) and (iPartner <= group.iEnd)    );
    int tTailPartner=confs.getChain(iPartner).tail;
    int tHeadPartner=confs.getChain(iPartner).head;
    #endif
    

    int iPartnerNext = confs.getChain(iPartner).next;
    

    int iNewChainHead=-1;
    if (tHead + l >= M)
    {
        iNewChainHead=confs.getChain(iPartner).prev;

        if (iNewChainHead == -1)
        {
            return false;
        }

       
       
    }
    else
    {
        iNewChainHead=iPartner;
        
    }


     if ( confs.getChain(iNewChainHead).tail >= tHead - 1  )
        {
            return false;
        }


    assert(iNewChainHead != -1);


    // metropolic test based on ratio of forward and backward move
    Real weightBackwardMove= 0;

    for(int i=group.iStart;i<=group.iEnd;i++)
    {
        Real norm=0;
        for(int d=0;d<getDimensions();d++)
        {
            distance[d]=geo.difference(  data(iNewChainHead,d,tHead) - data(i,d,tJoin),d) ;
        }


        weightBackwardMove+=exp(freeParticleLogProbability(distance,S.getTimeStep()*l,group.mass))*tags(i,tJoin);
    }
    

    Real deltaS=-(log(weightForwardMove) - log(weightBackwardMove)) ;
    
    bool acceptChain = metropolisSampler.acceptLog(-deltaS,randG);

    if (not acceptChain)
    {
        return false;
    }

    deltaS=0;

    
    const auto  partnerChain = confs.getChain(iPartner);


    deltaS-=Spot.evaluate(confs,{tHead,std::min(tHead + l ,M) -1},iNewChainHead);
    deltaS-=Spot.evaluate(confs,{0,tHead + l - M  - 1},iPartner);


    confs.copyDataToBuffer(buffer,{0,tHead + l - M  },iPartner);

    std::array<Real,3> difference;
    std::array<Real,getDimensions()> delta;

    // performs levy reconstruction between the head and the bead
    
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=
               geo.difference( 
                 - data(iChainHead,d,tHead) + data(iPartner,d,tJoin),d
            );
        data(iChainHead,d,tHead + l)=data(iChainHead,d,tHead) + difference[d];
        delta[d]=data(iChainHead,d,tHead + l) - data(iPartner,d,tJoin);
    }


    
    _levy.apply(confs,  {tHead,tHead + l },  iChainHead, S ,randG );
    confs.copyData( { M, tHead + l} , iChainHead, 0,iPartner );
    
    confs.setHead( iChainHead , std::min(tHead + l,M) );
    confs.setHead(iNewChainHead,tHead);

    bool accept= Spot.checkConstraints(confs,{tHead,std::min(tHead + l ,M) -1},iChainHead);
    accept= accept and Spot.checkConstraints(confs,{0,tHead + l - M  - 1},iPartner);

    if (accept)
    {
        deltaS+=Spot.evaluate(confs,{tHead,std::min(tHead + l ,M) -1},iChainHead);
        deltaS+=Spot.evaluate(confs,{0,tHead + l - M  - 1},iPartner);

        accept = metropolisSampler.acceptLog(-deltaS,randG);
    }

    ;
    if (accept) 
    {
        if (tHead + l >=M )
        {
            confs.join(iChainHead,iPartner );
            confs.translateData({tJoin + 1 , M}, {iPartner,iPartner},delta);
        }
        else
        {
            confs.copyData({tJoin+1,M},iPartner,iChainHead);
            confs.translateData({tJoin + 1 , M}, {iChainHead,iChainHead},delta);
            confs.setHead(iChainHead,M);
            confs.join(iChainHead,iPartnerNext);
        }
        assert( confs.getGroups()[getSet()].heads[0] == iNewChainHead);
        assert( confs.getGroups()[getSet()].heads.size() == 1);
        

    }
    else
    {
        confs.copyDataFromBuffer(buffer,{0,tHead + l - M  },iPartner);
        confs.setHead( iChainHead , tHead );
        confs.setHead( iNewChainHead , M );

        if (tHead + l >= M)
        {
            confs.join(iNewChainHead,iPartner);
        }
        else
        {
            confs.join(iNewChainHead,iPartnerNext);
        }

        #ifndef NDEBUG
        TESTCHAIN(confs,iChainHead,tHead,tTail);
        TESTCHAIN(confs,iPartner,tHeadPartner,tTailPartner);
        #endif

        assert( confs.getGroups()[getSet()].heads[0] == iChainHead);
        assert( confs.getGroups()[getSet()].heads.size() == 1);
        
    }    

    assert(group.isOpen());


    return accept;
}


std::ostream & sectorTableMoves::operator>> (std::ostream & os)
{
    for (int i=0;i<_moves.size();i++)
    {
        if (_nTrials[i] > 0 )
        {
            os << _names[i] << ":\t" << acceptanceRatio(i) << std::endl;
        }
        else
        {
            os << _names[i] << ":\t" << "-" << std::endl;
        }        
    }

    return os;
}


bool sectorTableMoves::attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG)
{

    if (_moves.size() == 0)
    {
        return false;
    }

    int iMove=sample(randG);


    auto  move = (_moves[ iMove ]);
    bool success=move->attemptMove(confs,S,randG);


    _nTrials[iMove]++;

    if (success)
    {
    _nSuccess[iMove]++;
    }

    return success;
};


void sectorTableMoves::resetCounters()
{
   
    for (int iMove=0;iMove< _nTrials.size() ; iMove++ )
    {
    _nTrials[iMove]=0;
    _nSuccess[iMove]=0;
    }
    
};

void tableMoves::resetCounters()
{
   
    for ( auto & tab: openTabs)
    {
        tab.resetCounters();
    }

    for ( auto & tab: closedTabs)
    {
        tab.resetCounters();
    }
};

void sectorTableMoves::push_back(move * move_,Real weight,const std::string & name)
{
    sampler.accumulateWeight(weight);
    _moves.push_back(move_);
    _names.push_back(name);
    _nTrials.push_back(0);
    _nSuccess.push_back(0);
    _weights.push_back(weight);

    totalWeight+=weight;


};


std::string reportMoves(const sectorTableMoves & tab)
{
    const auto & names = tab.getNames();
    const auto & weights = tab.getWeights();

    const auto totWeight = tab.getTotalWeight();

    std::stringstream ss;

    for (int i=0;i<names.size();i++)
    {
        ss << names[i] << " , probability: " << weights[i]/totWeight << std::endl;
    }

    return ss.str();
    
}



int sectorTableMoves::sample(randomGenerator_t & randG)
{
    int iMove = sampler.sample(randG);
    return iMove;
};

openMove::openMove(Real C_ , int set,int maxReconstructedLength_) : C(C_), _levy(maxReconstructedLength_+2) ,  _maxReconstructedLength(maxReconstructedLength_+2) ,buffer(2*(maxReconstructedLength_+2),getDimensions()),
gauss(0,1),uniformRealNumber(0,1),
setStartingBeadRandom(true),
setLengthCutRandom(true),
startingBead(-1),
lengthCut(-1),
startingChain(-1),
setStartingChainRandom(true),
singleSetMove(set)
{}


openMoveTest::openMoveTest(Real C_ , int set,int maxLength, int startingBead_) : C(C_), _levy(maxLength+2) ,
gauss(0,1),uniformRealNumber(0,1),
startingBead(startingBead_),
order(2),
singleSetMove(set)
{}


closeMoveTest::closeMoveTest(Real C_ , int set,int maxLength, int startingBead_) : C(C_), _levy(maxLength+2) ,
gauss(0,1),uniformRealNumber(0,1),
startingBead(startingBead_),
order(2),
singleSetMove(set)
{}



translateMove::translateMove(Real max_delta, int maxBeads, int set) : _max_delta(max_delta),buffer(maxBeads+1,getDimensions())  , distr(-1.,1.),singleSetMove(set)
{

}

bool translateMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{

    auto & Spot = S.getPotentialAction();
    std::array<int,2> timeRange={0,confs.nBeads()-1};

    // sample a chain
    int iChain = confSampler.sampleChain(confs,randG);
    auto & data = confs.dataTensor();

    if (iChain == -1)
    {
        return false;
    }
    auto chainsInThePolimer = confs.buildPolimerList(iChain); 


    Real deltaS = 0;

    // evaluate the action in the old configuration

     for (auto iCurrentChain : chainsInThePolimer)
     {
        deltaS-= Spot.evaluate(confs,timeRange,iCurrentChain);   
     }



    // translate the whole 
    for (int d=0;d<getDimensions();d++)
    {
        delta[d]=(distr(randG))*_max_delta;
    }

    int iSeq=0; // ith chain in the list
    for (auto iCurrentChain : chainsInThePolimer)
    {

        // save old data to buffer
        confs.copyDataToBuffer(buffer,{0,confs.nBeads()},iCurrentChain,iSeq*(confs.nBeads()+1));



        for(int t=0;t<=confs.nBeads();t++)
            {
                for (int d=0;d<getDimensions();d++)
                {
                    data(iCurrentChain,d,t)+=delta[d];
                }    

            }

        iSeq++; 
    }
    
    bool accept=true;
    for (auto iCurrentChain : chainsInThePolimer)
     {
        accept=  S.checkConstraints(confs,timeRange,{iCurrentChain,iCurrentChain});
        if (not accept)
        {
            break;
        }
     }


    if (accept)
    {
        //evaluate the action in the new configuration
        for (auto iCurrentChain : chainsInThePolimer)
        {
            deltaS+= Spot.evaluate(confs,timeRange,iCurrentChain);

        }



        accept = sampler.acceptLog(-deltaS,randG);
    }

    if ( not accept)
    {
        // copy back configurations in al chains of the permutation cycle
        int iSeq=0;
        for (auto iCurrentChain : chainsInThePolimer)
        {
          confs.copyDataFromBuffer(buffer,{0,confs.nBeads()},iCurrentChain,iSeq*(confs.nBeads()+1));
          iSeq++;
        }
    }

    return accept;
}

bool openMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    _levy.setReconstructorBoundaries(chainBoundary::fixed,chainBoundary::fixed);
    if ( confs.getEnsamble() == ensamble_t::canonical)
    {

        return attemptCanonicalMove(confs,S,randG);
    }
    else if ( confs.getEnsamble() == ensamble_t::grandCanonical)
    {
        
        return attemptGrandCanonicalMove(confs,S,randG);
    }
    else if ( confs.getEnsamble() == ensamble_t::semiGrandCanonical)
    {
        _levy.setReconstructorBoundaries(chainBoundary::free,chainBoundary::free);
        return attemptSemiGrandCanonicalMove(confs,S,randG);
    }

    return false;
}


bool openMove::attemptCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    Real timeStep = S.getTimeStep();

    if ( confs.isOpen(getSet()) )
    {
        throw invalidState("The configuration is already open.");
    }

    int iChain = confsSampler.sampleChain(confs,getSet(),randG);
    int iChainTail=confs.getChain(iChain).next;

    if (iChain < 0)
    {
        return false;
    }

    int iHead=confs.nBeads();
    
    //int l= std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ; // distance from itime where the head is formed

    int l = _maxReconstructedLength - 2;

    int t0= iHead - l;
    Real deltaS=0;
    std::array<int,2> timeRange={t0  , iHead -1 };

    auto & data = confs.dataTensor();

    std::array<Real,getDimensions()> difference;
    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);

    confs.copyDataToBuffer(buffer,{t0,iHead},iChain,0);

    // generates the head
    
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> headPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,t0);
    }

    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=
               geo.difference( 
                data(iChainTail,d,0)-data(iChain,d,t0),d
            );                                

            //data(iChain,d,iHead)-data(iChain,d,t0);

        if (    
            std::abs(data(iChain,d,iHead)-data(iChain,d,t0) ) > geo.getLBox(d)*0.5 )
            {
                std::cout << "Warning: too long segment for open chain" << std::endl;
                return false;
            }
        
        
   
    }

    
    Real mass = confs.getGroupByChain(iChain).mass;


    confsSampler.sampleFreeParticlePosition(headPosition,startPosition,timeStep*l,randG,mass);



    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,iHead)=data(iChain,d,t0)  + geo.difference( headPosition[d] - data(iChain,d,t0) ,d);
    }


    // perform levy reconstruction on l beads
    _levy.apply(confs,{t0,iHead},iChain,S,randG);

    // evaluates the action
    deltaS+=sPot.evaluate(confs,timeRange,iChain) ;

    // compute the acceptance ratio
  
    auto propRatio = -deltaS - freeParticleLogProbability(difference,S.getTimeStep()*l,mass) + log(C);


    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {

        confs.setHead(iChain,iHead);
    }
    else
    {
        confs.copyDataFromBuffer(buffer,{t0,iHead},iChain,0);
    }


    return accept;
};

createWorm::createWorm(Real C_, int set, int maxAdvanceLength_, Real sigma_) :
_maxReconstructedLength(maxAdvanceLength_+2) , _levy((maxAdvanceLength_+2)*2),gauss(0,sigma_),uniformRealNumber(0,1),singleSetMove(set),
sigma2(sigma_*sigma_),C(C_),
randomLength(true),randomInitialBead(true),
initialBead(0),
uniformDistributionFirstParticle(false)
{

}




int createWorm::sampleInitialBead(randomGenerator_t & randG, int M) 
{
    if (randomInitialBead)
    {
        return  std::floor(uniformRealNumber(randG) * M )  ;
    }
    else
    {
        return initialBead;
    }
    
}


int createWorm::sampleLength(randomGenerator_t & randG) 
{
    if (randomLength)
    {
        return std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ;
    }
    else
    {
        return (_maxReconstructedLength -2);
    }
}


bool createWorm::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    Real timeStep = S.getTimeStep();
    int M = confs.nBeads();
    const auto & geo = S.getGeometry();

    int l= sampleLength(randG);

    int lMax= _maxReconstructedLength -2;
    if (confs.isOpen(getSet()) )
    {
        throw std::runtime_error("Worm cannot be created in the open sector.");
    }
    #ifndef NDEBUG
    int oldSize = confs.getGroups()[getSet()].size();
    #endif

    auto & data = confs.dataTensor();
    
    
    int t0 = sampleInitialBead( randG,  M);
    int tTail = t0 - 1;


    //std::cout << t0 << " -> " << t0+l << std::endl;

    int iChainTail=confs.pushChain(getSet());
    int iChainHead=iChainTail;
    int tHead=t0+l;
    confs.setHeadTail(iChainTail,std::min(t0 + l,M),tTail);

    if (t0 + l  >= M )
    {
        iChainHead=confs.pushChain(getSet() );
        tHead=t0 + l - M; 
        confs.setHeadTail(iChainHead,tHead,-1);

    }

    // generates the tail according to a gaussian distribution
    std::array<Real,getDimensions() > x1;
    
    sampleFirstParticle(geo,x1,randG);


    // generates the head
    std::array<Real,getDimensions()> xl;

    Real mass = confs.getGroupByChain(iChainTail).mass;
    Real var=2*D*timeStep/ mass;

    confsSampler.sampleFreeParticlePosition(xl,x1,timeStep*l,randG,mass);

    // copy head and tail to the data
    for (int d=0;d<getDimensions();d++)
    {
        data(iChainTail,d,t0)=x1[d];
        data(iChainTail,d,t0+ l)=xl[d];   
    }

    auto & sPot = S.getPotentialAction();
    Real deltaS=0;


    std::array<int,2> timeRange{t0, std::min(t0 + l,M) - 1  };
    std::array<int,2> timeRange2{0,t0 + l - M -1};


    //deltaS-=sPot.evaluate(confs,timeRange,iChainTail);
    //deltaS-=sPot.evaluate(confs,timeRange2,iChainHead);

    _levy.apply(confs,{t0,t0 + l},iChainTail,S,randG);

    if (t0 + l >= M)
    {
        confs.copyData( { M , t0 + l } , iChainTail, 0,iChainHead ); 
    }

    bool accept=sPot.checkConstraints(confs,timeRange,iChainTail);

    if (accept)
    {
        accept=sPot.checkConstraints(confs,timeRange2,iChainHead);
    }

    if (accept)
    {
        deltaS+=sPot.evaluate(confs,timeRange,iChainTail);
        deltaS+=sPot.evaluate(confs,timeRange2,iChainHead);

        auto propRatio = -deltaS + confs.getChemicalPotential(getSet() )*l*timeStep - probabilityInitialPosition(geo,x1) + log(getCoefficientAccept(C,M,lMax ) );

        accept = sampler.acceptLog(propRatio,randG);
        
    }


    
    
    if ( accept)
    {
        
        if (t0+l >= M)
        {
            confs.join(iChainTail,iChainHead);
            #ifndef NDEBUG
            assert(confs.getGroups()[getSet() ].size() == oldSize + 2);
            #endif

        }
        else
        {
            #ifndef NDEBUG
            assert(confs.getGroups()[getSet() ].size() == oldSize + 1);
            TESTCHAIN(confs,iChainTail,tHead,tTail);
            #endif

        }
#ifndef NDEBUG        
        assert(confs.getChain(iChainTail).hasTail() );
        assert(confs.getChain(iChainHead).hasHead() );
        assert(confs.isOpen(getSet()));            
#endif
    }
    else
    {
        
        if (t0+l >= M)
        {
            confs.setTail(iChainHead,-1);

            confs.removeChains(iChainTail,iChainHead);
            
        }
        else
        {
            confs.removeChain(iChainHead);
            #ifndef NDEBUG
            
            TESTCHAIN(confs,iChainTail,tHead,tTail);
            
            #endif
            
        }

        
        #ifndef NDEBUG
        assert(not confs.isOpen(getSet()));
        assert(confs.getGroups()[getSet()].size() == oldSize );
        #endif

    }

    return accept;
}

deleteWorm::deleteWorm(Real C_, int set, int maxAdvanceLength_, Real sigma_) :
_maxReconstructedLength(maxAdvanceLength_+2) ,gauss(0,sigma_),uniformRealNumber(0,1),singleSetMove(set),
sigma2(sigma_*sigma_),C(C_),
randomLength(true),
randomInitialBead(true),
uniformDistributionFirstParticle(false)
{

}

bool deleteWorm::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    Real timeStep = S.getTimeStep();
    int M = confs.nBeads();
    const auto & geo = S.getGeometry();

    int lMax=_maxReconstructedLength -2;

    if (not confs.isOpen(getSet()) )
    {
        throw std::runtime_error("No worm to delete in closed sector.");
    }


    const auto & heads = confs.getGroups()[getSet()].heads;
    const auto & tails = confs.getGroups()[getSet()].tails;

    int iChainHead=heads[std::floor(uniformRealNumber(randG) * heads.size() )];
    int iChainTail=tails[std::floor(uniformRealNumber(randG) * tails.size() )];

    int tHead=confs.getChain(iChainHead).head;
    int tTail = confs.getChain(iChainTail).tail;

   
    // compute distance between head and tail
    int l = tHead - (tTail+1);

    if (iChainHead != iChainTail)
    {
        if ( confs.getChain(iChainTail).next != iChainHead ) // must be on same slice or successive slices
        {
            return false;
        }

        l += M;
  
    }

    // refuse if the chain is longer than the maximum chain which can be created in a createWorm move
    if (l > lMax )
    {
        return false;
    }

    auto & sPot = S.getPotentialAction();
    Real deltaS=0;

    int t0=tTail + 1;

    std::array<int,2> timeRange{t0, std::min(t0 + l,M) - 1  };
    std::array<int,2> timeRange2{0,t0 + l - M -1};


    deltaS-=sPot.evaluate(confs,timeRange,iChainTail);
    deltaS-=sPot.evaluate(confs,timeRange2,iChainHead);


    // copy t0 position to temporary to compute the probability
    auto & data = confs.dataTensor();
    std::array<Real,getDimensions()> x1;


    for (int d=0;d<getDimensions();d++)
    {
        x1[d]=data(iChainTail,d,t0);   
    }

    auto propRatio = -deltaS - confs.getChemicalPotential(getSet() )*l*timeStep + probabilityInitialPosition(geo,x1) - log(getCoefficientAccept(C,M,lMax) );
    ;

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        if ( iChainTail == iChainHead)
        {
            #ifndef NDEBUG
            TESTCHAIN(confs,iChainHead,tHead,tTail);
            #endif
            confs.removeChain(iChainHead);
        }
        else
        {
            confs.setTail(iChainHead,-1);
            confs.removeChains(iChainTail,iChainHead);
            
        }
    }
    else
    {

    }


    #ifndef NDEBUG
    int nAfterClose = confs.nParticles();

    #endif
    
    return accept;
}




Real openMove::openCloseRatioCoefficient(int N,int M)
    {
        Real coeff=C;
        
        if (setStartingChainRandom)
        {
            coeff*=N;
        }
        

        if (setStartingBeadRandom)
        {
            coeff*=M;
            
        }
        if (setLengthCutRandom)
        {
            coeff*=(_maxReconstructedLength - 2);
        }

        return coeff;
    }





bool openMove::attemptGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    Real timeStep = S.getTimeStep();

    int N = confs.nParticles(getSet());


  /*   std::cout << "Before open" << std::endl;
    confs >> std::cout;
    std::cout << std::endl; */

    const auto & geo = S.getGeometry();

    if ( confs.isOpen(getSet()) )
    {
        throw invalidState("The configuration is already open.");
    }

    int iChain=startingChain;

    if (setStartingBeadRandom)
    {
        iChain = confsSampler.sampleChain(confs,getSet(),randG);
    }

    if (iChain < 0)
    {
        return false;
    }



    int l = lengthCut;


    if (setLengthCutRandom)
    {
        l = std::floor(uniformRealNumber(randG) * (_maxReconstructedLength -2   ) ) + 1;
    }
    
    int M = confs.nBeads();


    int tHead = startingBead;


    if (setStartingBeadRandom)
    {
        

        tHead =  std::floor(  uniformRealNumber( randG )*confs.nBeads() );
    }
    

    int t1 = std::min(tHead+l,M);

    Real deltaS=0;
    std::array<int,2> timeRange={tHead  , t1 - 1 };

    auto & data = confs.dataTensor();

    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);

    // comput deltaS on the next chain ( if tTail overflows to the next chain)

    int t0_2=0;
    int t1_2=(tHead + l - M  );
    std::array<int,2> timeRange2={t0_2,t1_2-1};
    int iChainNext=confs.getChain(iChain).next;
    deltaS-=sPot.evaluate(confs,timeRange2,iChainNext);


    Real mass = confs.getGroupByChain(iChain).mass;

    int tTail = -1000;
    int iChainTail=-100;


    if (tHead+l < M)
    {
        tTail=timeRange[1];
        iChainTail=iChain;
    }
    else 
    {
        tTail=timeRange2[1];
        iChainTail=iChainNext;
    }


    std::array<Real,getDimensions()> difference; 
    for (int d=0;d<getDimensions();d++)
    {
        // difference along the ring (no pbc)
        difference[d]=
                data(iChainTail,d,tTail+1)-data(iChain,d,tHead);
         
        if ( tHead + l >= M)
        {
            difference[d]+= data(iChain,d,M)-data(iChainTail,d,0);
        }

           if (    
            std::abs(difference[d] ) >= geo.getLBox(d)*0.5 )
            {

                return false;
            } 

    }


    auto propRatio = -deltaS - freeParticleLogProbability(difference,S.getTimeStep()*l,mass)  -  confs.getChemicalPotential(getSet())*l*timeStep + log( openCloseRatioCoefficient(N,M) );


    bool accept = sampler.acceptLog(propRatio,randG);


    if ( accept)
    {
        if ( tHead + l >=M)
        {
            confs.setHead(iChain,tHead);
            confs.setTail(iChainTail,tTail);

#ifndef NDEBUG
            if (iChain != iChainTail)
            {
                TESTCHAIN(confs,iChain,tHead,-1);
                TESTCHAIN(confs,iChainTail,M,tTail);
            }
            else
            {
                TESTCHAIN(confs,iChain,tHead,tTail);
            }
#endif 

        }
        else
        {
            confs.setHead(iChain,tHead);
            int iChainTail=confs.pushChain(getSet() );
            confs.setHeadTail(iChainTail,M,tTail);
            confs.join(iChainTail,iChainNext);
            confs.copyData({tTail+1,M}  , iChain, iChainTail  );


#ifndef NDEBUG
                TESTCHAIN(confs,iChainTail,M,tTail);
                TESTCHAIN(confs,iChain,tHead,-1);
                TESTCHAINSJOINED(confs,iChainTail,iChainNext);
#endif 

        }
    }
    else
    {
        
    }

    /* std::cout << "After open" << std::endl;
    confs >> std::cout;
    std::cout << std::endl;
 */
    return accept;
};

bool openMove::attemptSemiGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    /* deleting from the setA and creating a worm in another set */


    Real timeStep = S.getTimeStep();

    // choose a matching set
     int setB= (getSet() + int( (confs.getGroups().size() -1) * uniformRealNumber(randG)  ) + 1 )   % confs.getGroups().size();

    //int setB= samplePartnerSet(randG) ;

    int setA=getSet();


    int NA = confs.nParticles(setA);
    int NB = confs.nParticles(setB);

    const auto & geo = S.getGeometry();

    if ( confs.isOpen(getSet()) )
    {
        throw invalidState("The configuration is already open.");
    }

    int iChainA=startingChain;

    if (setStartingBeadRandom)
    {
        iChainA = confsSampler.sampleChain(confs,getSet(),randG);
    }


    if (iChainA < 0)
    {
        return false;
    }
    

    int l = lengthCut;


    if (setLengthCutRandom)
    {
        l = std::floor(uniformRealNumber(randG) * (_maxReconstructedLength -2   ) ) + 1;
    }
    
    int M = confs.nBeads();


    int tHeadA = startingBead;


    if (setStartingBeadRandom)
    {
        

        tHeadA =  std::floor(  uniformRealNumber( randG )*confs.nBeads() );
    }
    
    // evaluate action on segment to delete in setA
    
    Real deltaS=0;
    std::array<int,2> timeRangeA={tHeadA , std::min(tHeadA+l,M) - 1 };

    auto & data = confs.dataTensor();

    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRangeA,iChainA);

    // comput deltaS on the next chain ( if tTail overflows to the next chain)
    int t0_2=0;
    int t1_2=(tHeadA + l - M  );
    std::array<int,2> timeRangeA2={t0_2,t1_2-1};
    int iChainANext=confs.getChain(iChainA).next;
    deltaS-=sPot.evaluate(confs,timeRangeA2,iChainANext);
    
    

    
    confs.deleteBeads(   timeRangeA, iChainA );
    confs.deleteBeads(   timeRangeA2, iChainANext );


    // creates the worm
    int tTailB=tHeadA - 1;
    int t0B = tTailB + 1;

    int iChainTailB=confs.pushChain(setB);
    int iChainHeadB=iChainTailB;
    int tHeadB=tTailB + l;
    confs.setHeadTail(iChainTailB,std::min(t0B + l,M),tTailB);
    if (t0B + l  > M)
    {
        iChainHeadB=confs.pushChain( setB );
        tHeadB=t0B + l - M; 
        confs.setHeadTail(iChainHeadB,tHeadB,-1);
    }


    _levy.apply(confs,{t0B,t0B + l},iChainTailB,S,randG);
    if (t0B + l > M)
    {
        confs.copyData( { M , t0B + l } , iChainTailB, 0,iChainHeadB ); 
    }

    std::array<int,2> timeRangeB{t0B, std::min(t0B + l,M) - 1  };
    std::array<int,2> timeRangeB2{0,t0B + l - M -1};

    // evaluate the action due to the presence of the worm
    deltaS+=sPot.evaluate(confs,timeRangeB,iChainTailB);
    deltaS+=sPot.evaluate(confs,timeRangeB2,iChainHeadB);

    
    // compute the distance between the head and the tail

    Real mass = confs.getGroupByChain(iChainA).mass;

    int tTailA = -1000;
    int iChainTailA=-100;


    if (tHeadA+l < M)
    {
        tTailA=timeRangeA[1];
        iChainTailA=iChainA;
    }
    else 
    {
        tTailA=timeRangeA2[1];
        iChainTailA=iChainANext;
    }

    std::array<Real,getDimensions()> difference; 
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=
               geo.difference( 
                data(iChainTailA,d,tTailA+1)-data(iChainA,d,tHeadA),d
            );
        
           if (    
            std::abs(data(iChainA,d,tTailA + 1)-data(iChainA,d,tHeadA) ) > geo.getLBox(d)*0.5 )
            {
                //std::cout <<  "Warning: l too long for open/close"<<std::endl;

                return false;
            } 
    }


    std::array<Real, getDimensions()> x0;
    for (int d=0;d<getDimensions();d++)
    {
        x0[d]=data(iChainTailB,d,t0B);
    }
    
    auto propRatio = -deltaS - freeParticleLogProbability(difference,S.getTimeStep()*l,mass)  -  ( confs.getChemicalPotential(setA) - confs.getChemicalPotential(setB)   )*l*timeStep + log( openCloseRatioCoefficient(NA,M) ) - _levy.probabilityInitialPosition(geo,x0);


    bool accept = sampler.acceptLog(propRatio,randG);


    if ( accept)
    {
        if ( tHeadA + l >=M)
        {
            confs.setHead(iChainA,tHeadA);
            confs.setTail(iChainTailA,tTailA);
        }
        else
        {
            confs.setHead(iChainA,tHeadA);
            int iChainTailA=confs.pushChain(getSet() );
            confs.setHeadTail(iChainTailA,M,tTailA);
            confs.join(iChainTailA,iChainANext);
            confs.copyData({tTailA+1,M}  , iChainA, iChainTailA  );
        }

         if (t0B +l > M)
        {
            confs.join(iChainTailB,iChainHeadB);
        }
        else
        {

        }

    }
    else
    {
        confs.createBeads(   timeRangeA, iChainA );
        confs.createBeads(   timeRangeA2, iChainANext );

        if (tTailB + 1 + l > M)
        {
            confs.setTail(iChainHeadB,-1);
            confs.removeChains(iChainTailB,iChainHeadB);
        }
        else
        {
            confs.removeChain(iChainHeadB);       
        }

        
    }

    return accept;
};


bool openMoveTest::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    Real timeStep = S.getTimeStep();


    int N = 1;



    const auto & geo = S.getGeometry();
    auto & data = confs.dataTensor();


    if ( confs.isOpen(getSet()) )
    {
        throw invalidState("The configuration is already open.");
    }

    int iChain = 0;


    int l = 1;


    
    int M = confs.nBeads();


    int tHead = M-1;
    
    
    Real mass = 1;

    Real deltaS=0;

    Real deltax2=0;
    std::array<Real,getDimensions()> difference; 
    for (int d=0;d<getDimensions();d++)
    {
        assert( std::abs(data(0,d,tHead+1) - data(0,d,0))<1e-9);

        difference[d]=data(0,d,tHead+1)-data(0,d,tHead);
        deltax2+=difference[d]*difference[d];

      
        deltaS-= 0.5 * 0.5 *( std::pow(data(0,d,tHead+1),2 ) + std::pow(data(0,d,tHead),2 )  )*timeStep;
      
        
    }
    
    Real var=timeStep;

    Real rhoFree= - 0.5 * deltax2 / var -getDimensions()*0.5 * log(2*M_PI*var) ;



    auto propRatio = -deltaS - rhoFree + log(C);

    bool accept = sampler.acceptLog(propRatio,randG);

    

    if ( accept)
    {  
        confs.setHead(0,tHead);
        //int iChainTail=confs.pushChain(getSet() );
        //assert(iChainTail == 1);
        //confs.setHeadTail(1,M,tHead);
        //confs.join(1,0);
        //confs.copyData({tHead+1,M}  , 0, 1  );

    }

    return accept;
};

bool closeMoveTest::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    Real timeStep = S.getTimeStep();

    int N = 1;

    const auto & geo = S.getGeometry();
    auto & data = confs.dataTensor();

    if ( ! confs.isOpen(getSet()) )
    {
        throw invalidState("The configuration is already closed.");
    }

    int iChainHead = 0;
    int iChainTail = 1;


    int l = 1;


    
    int M = confs.nBeads();

    int tHead = M-1;
    
    
    Real mass = 1;

    Real deltaS=0;

    std::array<Real,getDimensions()> difference;
    Real deltax2=0; 
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=data(0,d,0)-data(0,d,tHead);
        deltax2+=difference[d]*difference[d];
        deltaS+= 0.5 * 0.5 *( std::pow(data(0,d,0),2 ) + std::pow(data(0,d,tHead),2 )  )*timeStep;
        


        
    }
    

    Real var=timeStep;

    Real rhoFree= - 0.5 * deltax2 / var -0.5 * getDimensions() * log(2*M_PI*var) ;


    auto propRatio = -deltaS + rhoFree   - log( C );


    bool accept = sampler.acceptLog(propRatio,randG);


    if ( accept)
    {  

        confs.setHead(0,M);
        //confs.setHeadTail(0,M,-1);
        //confs.copyData({tHead+1,M}  , 1, 0  );
        
        for (int d=0;d<getDimensions();d++)
        {
            data(0,d,M)=data(0,d,0);
        }
        confs.join(0,0);

        //confs.removeChain(1);
    }

    return accept;
};



closeMove::closeMove(Real C_ , int set,int maxReconstructionLength) : C(C_),_levy(2*(maxReconstructionLength+2)),_maxLength(maxReconstructionLength+2),buffer((maxReconstructionLength+2)*2,getDimensions()),gauss(0,1),uniformRealNumber(0,1),
setStartingBeadRandom(true),
setLengthCutRandom(true),
startingBead(-1),
lengthCut(-1),
startingChain(-1),
setStartingChainRandom(true),
singleSetMove(set),
enforceAtomeNumberOnClose(false),
atomNumberOnClose(0)
{}


bool closeMove::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    if ( confs.getEnsamble() == ensamble_t::canonical)
    {
        return attemptCanonicalMove(confs,S,randG);
    }
    else if ( confs.getEnsamble() == ensamble_t::grandCanonical)
    {
        return attemptGrandCanonicalMove(confs,S,randG);
    }
    else if ( confs.getEnsamble() == ensamble_t::semiGrandCanonical)
    {
        return attemptSemiGrandCanonicalMove(confs,S,randG);
    }

    return false;
}


bool closeMove::attemptCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    
     if (! confs.isOpen(getSet()) )
    {
        throw invalidState("The configuration is already closed.");
    }


    const auto & heads = confs.getGroups()[getSet()].heads;
    const auto & tails = confs.getGroups()[getSet()].tails;

    

    int iChainHead=heads[std::floor(uniformRealNumber(randG) * heads.size() )];
    int iChainTail=tails[std::floor(uniformRealNumber(randG) * tails.size() )];
        
    auto timeStep = S.getTimeStep();
    auto & geo = S.getGeometry();

    int iHead=confs.nBeads();
    int iTail = 0;

    
    
    //int l= std::floor( uniformRealNumber(randG) * (_maxLength -2) ) + 1 ; // distance from itime where the head is formed
    int l = _maxLength - 2;
    
    auto & sPot = S.getPotentialAction();


    Real deltaS=0;
    int t0= iHead - l;

    std::array<int,2> timeRange={iHead - l , iHead -1 };

    auto & data = confs.dataTensor();

    //std::cout << data(iChainHead,0,iHead)-data(iChainTail,0,iTail) << std::endl;


    std::array<Real,getDimensions()> difference;

    deltaS-=sPot.evaluate(confs,{t0,iHead-1},iChainHead);

    confs.copyDataToBuffer(buffer,{t0,iHead },iChainHead,0);

    // copy first valid bead of tail of the head
    Real distanceSquared=0;
    std::array<Real,3> headPosition;
    std::array<Real,3> startPosition;
    std::array<Real,3> distance;

    for (int d=0;d<getDimensions();d++)
    {
        //std::cout << "close: " <<geo.difference( - data(iChainHead,d,t0) + data(iChainTail,d,iTail),d) << std::endl;

        
        data(iChainHead,d, iHead )=data(iChainHead,d,t0)  + geo.difference( -  data(iChainHead,d,t0)  + data(iChainTail,d,iTail),d);



        //data(iChainHead,d, iHead )=data(iChainTail,d,iTail);


        //data(iChainHead,d, iHead )+=data(iChainHead,d,t0);


        //data(iChainHead,d, iHead )= - data(iChainHead,d,t0) + data(iChainTail,d,iTail);

        
        
        //std::cout<<  ( - data(iChainHead,d,t0) + data(iChainHead,d,iHead),d) << std::endl;

        //std::cout <<   - data(iChainHead,d,0) + data(iChainHead,d,iHead) << std::endl;



    }


    Real mass = confs.getGroupByChain(iChainHead).mass;
    

    // perform levy reconstruction on l beads
    _levy.apply(confs,{t0,iHead},iChainHead,S,randG);

    // evaluates the action
    deltaS+=sPot.evaluate(confs,{t0,iHead-1},iChainHead) ;



    // compute the acceptance ratio
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=
            //geo.difference( 
            //    data(iChainHead,d,t0)-data(iChainHead,d,iHead),d
            //);
         data(iChainHead,d,t0)-data(iChainHead,d,iHead);
    }

/* 
     int winding =  (data(iChainHead,0,iHead)-data(iChainTail,0,iTail))/2;
     
     {
        std::cout << winding << std::endl;
     } */



    auto propRatio = -deltaS + freeParticleLogProbability(difference,S.getTimeStep()*l,mass) -log(C);


    //std::cout << propRatio << std::endl;


    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        confs.join(iChainHead,iChainTail);
        //confs.fillHead(iChainHead);

        /* for(int d=0;d<getDimensions();d++)
        {
            std::cout<< - data(iChainHead,d,iHead) + data(iChainTail,d,iTail) << " ";
        }
        std::cout << std::endl; */
        
    }
    else
    {
        confs.copyDataFromBuffer(buffer,{t0,iHead },iChainHead,0);
    }
    
    return accept;

};

Real closeMove::openCloseRatioCoefficient(int N,int M)
    {
        Real coeff=C;
        
        if (setStartingChainRandom)
        {
            coeff*=N;
        }
        
        if (setStartingBeadRandom)
        {
            coeff*=M;
            
        }
        if (setLengthCutRandom)
        {
            coeff*=(_maxLength - 2 );
        }



        return coeff;
    }

void closeMove::setAtomNumberOnClose(size_t N,const std::vector<int> & sets)
{
    setsMaxParticleNumber=sets;
    enforceAtomeNumberOnClose=true;
    atomNumberOnClose=N;


}

bool closeMove::checkConstraintsOnClose(const configurations_t & configurations)
{

    if (enforceAtomeNumberOnClose)
    {
        int N=0;

        for (auto set : setsMaxParticleNumber)
        {
            N+=nParticlesOnClose(configurations, set );
        };

        if (N != atomNumberOnClose)
        {
            return false;
        }
        else
        {
            return true;

        }
    }
    else
    {
        return true;
    }    
}

bool closeMove::attemptGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
  /*   std::cout << "Before close" << std::endl;
    confs >> std::cout;
    std::cout << std::endl; */


    if ( not checkConstraintsOnClose(confs))
    {
        return false;
    }

     if (! confs.isOpen(getSet()) )
    {
        throw invalidState("The configuration is already closed.");
    }

    const auto & heads = confs.getGroups()[getSet()].heads;
    const auto & tails = confs.getGroups()[getSet()].tails;


    int iChainHead=heads[std::floor(uniformRealNumber(randG) * heads.size() )];
    int iChainTail=tails[std::floor(uniformRealNumber(randG) * tails.size() )];
    
    auto timeStep = S.getTimeStep();

    auto & geo = S.getGeometry();

    int tHead=confs.getChain(iChainHead).head;

    int tTail = confs.getChain(iChainTail).tail;

    int lMax = _maxLength - 2;

    int M = confs.nBeads();
    
    int N = confs.nParticles(getSet() );


    //int l= std::floor( uniformRealNumber(randG) * (_maxLength -2) ) + 1 ; // distance from itime where the head is formed
    int l = tTail - tHead + 1 ;

    if (l<=0)
    {
        l+= M;
    }

    if (l > lMax)
    {
        return false;
    }

    auto & sPot = S.getPotentialAction();

    Real deltaS=0;
    int t0= tHead;
    int t1=std::min(tHead + l,M); 

    std::array<int,2> timeRange={t0,t1-1 };

    auto & data = confs.dataTensor();

    //std::cout << data(iChainHead,0,iHead)-data(iChainTail,0,iTail) << std::endl;


    //deltaS-=sPot.evaluate(confs,timeRange,iChainHead);

    t0=0;
    t1=tHead + l - M  ; 
    std::array<int,2> timeRange2={t0,t1-1 };


    //deltaS-=sPot.evaluate(confs,timeRange2,iChainTail);
    std::array<Real,3> diff;

    std::array<Real,getDimensions()> difference;
    std::array<Real,getDimensions()> delta;
    

    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=
               geo.difference( 
                data(iChainTail,d,tTail + 1)-data(iChainHead,d,tHead),d
            );
        data(iChainHead,d,tHead + l) =                           data(iChainHead,d,tHead) + difference[d];
        delta[d]=  data(iChainHead,d,tHead + l) -  data(iChainTail,d,tTail + 1);  
    }

    _levy.apply(confs,{tHead,tHead + l},iChainHead,S,randG);

    configurations_t::copyData(confs, {M,tHead + l}, {iChainHead,iChainHead}, confs, 0 , iChainTail   );



    confs.setHead(iChainHead,std::min(tHead+l,M));

    if (tHead + l >= M)
    {
       
        confs.setTail(iChainTail,-1);
              
        confs.join(iChainHead,iChainTail);
    }

    //confs.fillHeads();

    bool accept = sPot.checkConstraints(confs,timeRange , {iChainHead,iChainHead});
    accept = accept and  sPot.checkConstraints(confs,timeRange2 , {iChainTail,iChainTail});
    
    if (accept )
    {
        deltaS+=sPot.evaluate(confs,timeRange,iChainHead);
        deltaS+=sPot.evaluate(confs,timeRange2,iChainTail);

        Real mass = confs.getGroupByChain(iChainHead).mass;

        int nNew = N + 1;

        if ( (tHead + l >=M) and (iChainHead != iChainTail ) )
        {    
            nNew= N + 2;    
        }
        

        auto propRatio = -deltaS + freeParticleLogProbability(difference,S.getTimeStep()*l,mass) -log( openCloseRatioCoefficient(nNew,M)   ) + confs.getChemicalPotential(getSet() )*l*timeStep ;

        accept = sampler.acceptLog(propRatio,randG);

    }


   

    if ( accept)
    {
        
        if (( tHead + l) < M)
        {
            confs.copyData({tTail + 1,M},iChainTail,iChainHead);
            confs.translateData({tTail + 1,M},{iChainHead,iChainHead},delta);
            

            int iChainNext=confs.getChain(iChainTail).next;
            confs.setHead(iChainHead,M);
            confs.join(iChainHead,iChainNext);
            #ifndef NDEBUG
            TESTCHAIN(confs,iChainHead,M,-1);
            TESTCHAINSJOINED(confs,iChainHead,iChainNext);
            #endif
            confs.removeChain(iChainTail);
            
        }
        else
        {
            confs.translateData({tTail + 2, M},{iChainTail,iChainTail},delta);

#ifndef NDEBUG
        TESTCHAIN(confs,iChainHead,M,-1);
        TESTCHAIN(confs,iChainTail,M,-1);
        TESTCHAINSJOINED(confs,iChainHead,iChainTail);
#endif
        }
    }
    else
    {
        confs.setHead(iChainHead,tHead);
        confs.setTail(iChainTail,tTail);
        if (tHead + l >= M)
        {
            std::array<Real,getDimensions()> deltaOpposite;
    
            for(int d=0;d<getDimensions();d++)
                {
                deltaOpposite[d]=-delta[d];
            }
            
            confs.translateData({tTail + 1,tTail + 1},{iChainTail,iChainTail},deltaOpposite);
        }
        
        
    }

   /*  std::cout << "After close" << std::endl;
    confs >> std::cout;
    std::cout << std::endl; */

    return accept;

};








bool closeMove::attemptSemiGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{


     if (! confs.isOpen(getSet()) )
    {
        throw invalidState("The configuration is already closed.");
    }

    const auto & headsA = confs.getGroups()[getSet()].heads;
    const auto & tailsA = confs.getGroups()[getSet()].tails;



    int iChainHeadA=headsA[std::floor(uniformRealNumber(randG) * headsA.size() )];
    int iChainTailA=tailsA[std::floor(uniformRealNumber(randG) * tailsA.size() )];

    int tHeadA=confs.getChain(iChainHeadA).head;

    int tTailA = confs.getChain(iChainTailA).tail;

    int setB= (getSet() + int( (confs.getGroups().size() -1) * uniformRealNumber(randG)  ) + 1 )   % confs.getGroups().size();


    //int setB= samplePartnerSet(randG) ;


    const auto & headsB = confs.getGroups()[setB].heads;
    const auto & tailsB = confs.getGroups()[setB].tails;

    int iChainHeadB=headsB[std::floor(uniformRealNumber(randG) * headsB.size() )];
    int iChainTailB=tailsB[std::floor(uniformRealNumber(randG) * tailsB.size() )];


    int tHeadB=confs.getChain(iChainHeadB).head;
    int tTailB = confs.getChain(iChainTailB).tail;

     if (iChainHeadB != iChainTailB )
    {
        assert( confs.getChain(iChainTailB).next == iChainHeadB );

    }






    
    auto timeStep = S.getTimeStep();

    auto & geo = S.getGeometry();

    
    int lMax = _maxLength - 2;

    int M = confs.nBeads();
    
    int NA = confs.nParticles(getSet() );


    //int l= std::floor( uniformRealNumber(randG) * (_maxLength -2) ) + 1 ; // distance from itime where the head is formed
    int l = tTailA - tHeadA + 1 ;

    if (l<=0)
    {
        l+= M;
    }

    if (l > lMax)
    {
        return false;
    }


    auto & sPot = S.getPotentialAction();

    Real deltaS=0;
    int t0A= tHeadA;
    int t1A=std::min(tHeadA + l,M); 

    std::array<int,2> timeRangeA={t0A,t1A-1 };

    auto & data = confs.dataTensor();

    //std::cout << data(iChainHead,0,iHead)-data(iChainTail,0,iTail) << std::endl;


    deltaS-=sPot.evaluate(confs,timeRangeA,iChainHeadA);

    t0A=0;
    t1A=tHeadA + l - M  ; 
    std::array<int,2> timeRangeA2={t0A,t1A-1 };


    deltaS-=sPot.evaluate(confs,timeRangeA2,iChainTailA);


    int t0B=tTailB + 1;

    std::array<int,2> timeRangeB{t0B, std::min(t0B + l,M) - 1  };
    std::array<int,2> timeRangeB2{0,t0B + l - M -1};


    deltaS-=sPot.evaluate(confs,timeRangeB,iChainTailB);
    deltaS-=sPot.evaluate(confs,timeRangeB2,iChainHeadB);

    








    std::array<Real,3> diff;

    std::array<Real,getDimensions()> difference;
    std::array<Real,getDimensions()> delta;
    

    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=
               geo.difference( 
                data(iChainTailA,d,tTailA + 1)-data(iChainHeadA,d,tHeadA),d
            );
        data(iChainHeadA,d,tHeadA + l) =                           data(iChainHeadA,d,tHeadA) + difference[d];
        delta[d]=  data(iChainHeadA,d,tHeadA + l) -  data(iChainTailA,d,tTailA + 1);  
    }

    _levy.apply( confs,{tHeadA,tHeadA + l}, iChainHeadA,S,randG);

    configurations_t::copyData(confs, {M,tHeadA + l}, {iChainHeadA,iChainHeadA} , confs, 0 , iChainTailA   );

    confs.setHead(iChainHeadA, std::min(tHeadA+l,M));

    if (tHeadA + l >= M)
    {
       
        confs.setTail(iChainTailA,-1);
        confs.join(iChainHeadA,iChainTailA);
    }

    deltaS+=sPot.evaluate(confs,timeRangeA,iChainHeadA);
    deltaS+=sPot.evaluate(confs,timeRangeA2,iChainTailA);


    Real mass = confs.getGroupByChain(iChainHeadA).mass;

    int nNewA = NA + 1;

    if ( (tHeadA + l >=M) and (iChainHeadA != iChainTailA ) )
    {
        
        nNewA= NA + 2;    
    }
    


    std::array<Real,getDimensions()> x0B;
    for (int d=0;d<getDimensions();d++)
    {
        x0B[d]=data(iChainTailB,d,t0B);   
    }


    auto propRatio = -deltaS + freeParticleLogProbability(difference,S.getTimeStep()*l,mass) -log( openCloseRatioCoefficient(nNewA,M)   ) + ( confs.getChemicalPotential( getSet() ) - confs.getChemicalPotential( setB )  )*l*timeStep + _levy.probabilityInitialPosition(geo,x0B) ;



    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        
        if (( tHeadA + l) < M)
        {
            confs.copyData({tTailA + 1,M},iChainTailA,iChainHeadA);
            confs.translateData({tTailA + 1,M},{iChainHeadA,iChainHeadA},delta);
            

            int iChainNextA=confs.getChain(iChainTailA).next;
            confs.setHead(iChainHeadA,M);
            confs.join(iChainHeadA,iChainNextA);
            confs.removeChain(iChainTailA);
            
        }
        else
        {
            confs.translateData({tTailA + 2, M},{iChainTailA,iChainTailA},delta);

        }


         if ( iChainTailB == iChainHeadB)
        {
            confs.removeChain(iChainHeadB);
        }
        else
        {
            confs.setTail(iChainHeadB,-1);
            confs.removeChains(iChainTailB,iChainHeadB);
        }

    }
    else
    {
        confs.setHead(iChainHeadA,tHeadA);
        confs.setTail(iChainTailA,tTailA);
        if (tHeadA + l >= M)
        {
            std::array<Real,getDimensions()> deltaOpposite;
    
            for(int d=0;d<getDimensions();d++)
                {
                deltaOpposite[d]=-delta[d];
            }
            
            confs.translateData({tTailA + 1,tTailA + 1},{iChainTailA,iChainTailA},deltaOpposite);
        }

        
    }

   /*  std::cout << "After close" << std::endl;
    confs >> std::cout;
    std::cout << std::endl; */

    return accept;

};


advanceHead::advanceHead(int maxAdvanceLength_,int set) :
_maxReconstructedLength(maxAdvanceLength_+2) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(set) ,
enforceMaxParticleNumber(false),_nMax(1e+9)
{
    setRandomLength();
}

recedeTail::recedeTail(int maxAdvanceLength_,int set) :
_maxReconstructedLength(maxAdvanceLength_+2) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(set)
{
    setRandomLength();
}


advanceHeadTest::advanceHeadTest(int l_) :
 _levy( l_*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(0),l(l_)
{

}

recedeHeadTest::recedeHeadTest(int l_) :
 _levy( l_*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(0),l(l_)
{

}


int advanceHead::sampleLength(randomGenerator_t & randG) 
{
    if (_setRandomLength)
    {

        return  std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ; // distance from itime where the head is formed
    }
    else
    {
        return _maxReconstructedLength - 2;
    }
}


int recedeTail::sampleLength(randomGenerator_t & randG) 
{
    if (_setRandomLength)
    {

        return  std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ; // distance from itime where the head is formed
    }
    else
    {
        return _maxReconstructedLength - 2;
    }
}


int recedeHead::sampleLength(randomGenerator_t & randG) 
{
    if (_setRandomLength)
    {

    return  std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ; // distance from itime where the head is formed
    }
    else
    {
        return _maxReconstructedLength - 2;
    }
}


int advanceTail::sampleLength(randomGenerator_t & randG) 
{
    if (_setRandomLength)
    {
    return  std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ; // distance from itime where the head is formed
    }
    else
    {
        return _maxReconstructedLength - 2;
    }
}

bool advanceHead::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();
    
    const auto & heads = confs.getGroups()[getSet()].heads;

    int iChain=heads[std::floor(uniformRealNumber(randG) * heads.size() )];
    

    int tHead=confs.getChain(iChain).head;



    int M = confs.nBeads();
    
    int l= sampleLength(randG);
    
    Real deltaS=0;
    
    std::array<int,2> timeRange={tHead  , std::min(tHead + l,M) - 1 };
    std::array<int,2> timeRange2={0  , tHead + l - M  - 1 };


    auto & data = confs.dataTensor();

    std::array<Real,3> difference;
    
    auto & sPot = S.getPotentialAction();

    //deltaS-=sPot.evaluate(confs,timeRange,iChain);


    // generates the head
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> headPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,tHead);
    }

    Real mass = confs.getGroupByChain(iChain).mass;
    Real var=2*D*timeStep/ mass;

    confsSampler.sampleFreeParticlePosition(headPosition,startPosition,timeStep*l,randG,mass);

    assert(l>0);



    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,tHead+ l)=headPosition[d];
    }
    
    confs.setHead(iChain,std::min(tHead+l,M));

    // performs levy reconstruction on l beads
    _levy.apply(confs,{tHead,tHead+l},iChain,S,randG);

    
    
    bool accept=sPot.checkConstraints(confs,timeRange,iChain);
    int iChainNew=-1100;

    if (accept )
    {
        deltaS+=sPot.evaluate(confs,timeRange,iChain);


        if (tHead + l >= M)
        {
            // creates a new chain
            iChainNew=confs.pushChain(getSet());
            confs.setHeadTail(iChainNew,timeRange2[1] + 1,-1);
            confs.copyData( { M , tHead + l } , iChain, 0,iChainNew );
            accept= sPot.checkConstraints(confs,timeRange2,iChainNew);
            if (accept)
            {
                deltaS+=sPot.evaluate(confs,timeRange2,iChainNew);

            }
            confs.join(iChain,iChainNew);
        }

    }

    if (accept)
    {
        auto propRatio = -deltaS + confs.getChemicalPotential(getSet())*l*timeStep;
        accept=sampler.acceptLog(propRatio,randG);
    }
    

    if ( enforceMaxParticleNumber)
    {   
        int nParticlesAfterClose=nParticlesOnClose(confs,getSet());
        if ( nParticlesAfterClose > _nMax )
        {
            accept=false;
        }
    }


    if ( accept)
    {
        
    }
    else
    {
        confs.setHead(iChain,tHead);
        if (tHead + l >= M)
        {
            confs.removeChain(iChainNew);
        }
    }

    return accept;
}

bool recedeTail::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();

    const auto & tails = confs.getGroups()[getSet()].tails;

    int iChain=tails[std::floor(uniformRealNumber(randG) * tails.size() )];

    int tTail=confs.getChain(iChain).tail;

    int M = confs.nBeads();
    
    int l= sampleLength(randG);
    
    Real deltaS=0;

    int iChainNewTail = iChain;
    int tNewTail = tTail - l;

    if (tTail + 1 - l < 0)
    {
        iChainNewTail=confs.pushChain( getSet());
        tNewTail=tTail - l + M;
        confs.setHeadTail(iChainNewTail,M,tNewTail);
    }

    std::array<int,2> timeRange={ std::max(0, tTail +1  - l )  , tTail };
    std::array<int,2> timeRange2={ tTail + 1 - l + M  , M - 1 };

    auto & data = confs.dataTensor();    

    std::array<Real,3> difference;

    auto & sPot = S.getPotentialAction();
    // generates the head
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> tailPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,tTail+1);
    }

    Real mass = confs.getGroupByChain(iChain).mass;
    Real var=2*D*timeStep/ mass;

    confsSampler.sampleFreeParticlePosition(tailPosition,startPosition,timeStep*l,randG,mass);

    assert(l>0);

    for (int d=0;d<getDimensions();d++)
    {
        data(iChainNewTail,d,tNewTail + 1)=tailPosition[d];
        data(iChainNewTail,d,tNewTail + l + 1)=data(iChain,d,tTail + 1);
    }



    confs.setTail(iChain,std::max(tTail- l,-1));

    // performs levy reconstruction on l beads
    _levy.apply(confs,{tNewTail+1,tNewTail+l+1},iChainNewTail,S,randG);

    if (tTail + 1 - l < 0)
    {
        confs.copyData( { M ,  tNewTail + l + 1 } , iChainNewTail, 0, iChain ); 
        deltaS+=sPot.evaluate(confs,timeRange2,iChainNewTail);
    }

    deltaS+=sPot.evaluate(confs,timeRange,iChain);

    auto propRatio = -deltaS + confs.getChemicalPotential(getSet() )*l*timeStep;


    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        if (tTail  +  1 - l < 0 )
        {
            confs.join(iChainNewTail,iChain);
        }
    }
    else
    {
        confs.setTail(iChain,tTail);

        if (tTail + 1 - l < 0 )
        {
            confs.removeChain(iChainNewTail);
        }
    }

    return accept;
}


bool recedeHeadTest::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();

    const auto & heads = confs.getGroups()[0].heads;


    int iChain=heads[std::floor(uniformRealNumber(randG) * heads.size() )];
    

    int tHead=confs.getChain(iChain).head;

    int M = confs.nBeads();

    
    Real deltaS=0;

    
    

    auto & data = confs.dataTensor();

    
    auto & sPot = S.getPotentialAction();

    

    deltaS-=sPot.evaluate(confs,{tHead - l ,tHead},iChain);

    

    auto propRatio = -deltaS - confs.getChemicalPotential(getSet())*l*timeStep;


    bool accept = sampler.acceptLog(propRatio,randG);
    
    if ( accept)
    {
        confs.setHead(iChain,tHead - l);
    }
    else
    {
        
    }

    return accept;
}


bool advanceHeadTest::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();

    const auto & heads = confs.getGroups()[0].heads;


    int iChain=heads[std::floor(uniformRealNumber(randG) * heads.size() )];
    

    int tHead=confs.getChain(iChain).head;

    int M = confs.nBeads();

    
    Real deltaS=0;
    

    auto & data = confs.dataTensor();

    auto & sPot = S.getPotentialAction();

    

    // generates the head
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> headPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,tHead);
    }

    Real mass = confs.getGroupByChain(iChain).mass;
    Real var=2*D*timeStep/ mass;

    confsSampler.sampleFreeParticlePosition(headPosition,startPosition,timeStep*l,randG,mass);
    
    assert(l>0);


    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,tHead+ l)=headPosition[d];
    }
    
    confs.setHead(iChain,tHead + l);

    // performs levy reconstruction on l beads
    _levy.apply(confs,{tHead,tHead+l},iChain,S,randG);


    deltaS+=sPot.evaluate(confs,{tHead,tHead+l},iChain);
    

    auto propRatio = -deltaS + confs.getChemicalPotential(getSet() )*l*timeStep;

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        
    }
    else
    {
        confs.setHead(iChain,tHead);
    }

    return accept;
}

recedeHead::recedeHead(int maxAdvanceLength_,int set) :
_maxReconstructedLength(maxAdvanceLength_+2) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(set),
_nMin(0),enforceMinParticleNumber(false)
{
    setRandomLength();   
}


advanceTail::advanceTail(int maxAdvanceLength_,int set) :
_maxReconstructedLength(maxAdvanceLength_+2) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(set)
{
    setRandomLength();
}

bool recedeHead::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();
    const auto & heads = confs.getGroups()[getSet()].heads;

    int iChain=heads[std::floor(uniformRealNumber(randG) * heads.size() )];    

    assert(iChain !=-1);

    int tHead=confs.getChain(iChain).head;
    
    int M = confs.nBeads();

    int l= sampleLength(randG);

    // do not accept propose which would destry the chain
    int tTail = -100000;
    int iChainTail=-100000;
    if (confs.getChain(iChain).hasTail() )
    {
        tTail=confs.getChain(iChain).tail;
        iChainTail=iChain;

        if (tHead - l <= tTail + 1)
        {
            return false;
        }
    }
    else
    {
        iChainTail=confs.getChain(iChain).prev;
        tTail=confs.getChain(iChainTail).tail;

        if (tHead - l + M <= tTail + 1)
        {
            return false;
        }

    }



    Real deltaS=0;

    std::array<int,2> timeRange={ std::max(tHead - l,0 ) , tHead - 1  };
    std::array<int,2> timeRange2={ tHead - l + M  , M - 1};


    auto & data = confs.dataTensor();

    std::array<Real,3> difference;
    
    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);

    if (tHead - l < 0 )
    {
        deltaS-=sPot.evaluate(confs,timeRange2,iChainTail);
    }

    confs.setHead(iChain,std::max(tHead-l,0));

    // evaluates the action
    //deltaS+=sPot.evaluate(confs,timeRange,iChain);

    if (tHead - l < 0 )
    {
        confs.setHead(iChainTail,tHead - l + M);
        //deltaS+=sPot.evaluate(confs,timeRange2,iChainTail);
    }

    auto propRatio = -deltaS -  confs.getChemicalPotential(getSet() )*l*timeStep;


    bool accept = sampler.acceptLog(propRatio,randG);
/* 
    if ( (iChain<=0) and (tHead<=9) )
    {
        accept=false;
    } */



     if ( accept)
    {
        if (tHead - l < 0 )
        {
            confs.join(iChain,iChain);
        }
    }
    
    if ( enforceMinParticleNumber)
    {   
        int nParticlesAfterClose=nParticlesOnClose(confs,getSet());
        if (tHead - l < 0 )
        {
            nParticlesAfterClose--;
        }

        if ( nParticlesAfterClose < _nMin )
        {
            accept=false;
        }
    }

    if ( accept)
    {
        if (tHead - l < 0 )
        {
            confs.setHeadTail(iChain,M-1,-1);
            confs.removeChain(iChain);
        }
    }
    else
    {
       
        if (tHead - l >= 0)
        {
            confs.setHead(iChain,tHead);
        }
        else
        {
            confs.setHead(iChain,tHead);
            confs.setHead(iChainTail,M);
            confs.join(iChainTail,iChain);
        }
    }

    return accept;

}


bool advanceTail::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();
    const auto & tails = confs.getGroups()[getSet()].tails;

    int iChain=tails[std::floor(uniformRealNumber(randG) * tails.size() )];

    assert(iChain !=-1);

    int tTail=confs.getChain(iChain).tail;

    int M = confs.nBeads();

    int l= sampleLength(randG);
    int tHead=confs.getChain(iChain).head;

    // do not accept propose which would destroy the chain
    int iChainNewTail=iChain;
    int tNewTail=(tTail + 1 +  l)%M -1;

    if (confs.getChain(iChain).hasHead() )
    {
        if ( tTail  + l >= tHead - 1)
        {
            return false;
        }
    }
    else
    {
        int iChainNext=confs.getChain(iChain).next;

        if (tTail + l + 1 >= M )
        {
            iChainNewTail=iChainNext;
            if (tNewTail >= confs.getChain(iChainNext).head - 1  )
            {
                return false;
            }
        }
        
    }



    Real deltaS=0;

    std::array<int,2> timeRange={ tTail + 1 , std::min(tTail + l + 1 ,M) - 1 };
    std::array<int,2> timeRange2={ 0  , tNewTail };

    auto & data = confs.dataTensor();
    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);
    if ( tTail + 1 + l >=M )
    {
        deltaS-=sPot.evaluate(confs,timeRange2,iChainNewTail);
    }

    auto propRatio = -deltaS -  confs.getChemicalPotential(getSet() )*l*timeStep;

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        confs.setTail(iChainNewTail,tNewTail);
        if (tTail + 1 + l >= M )
        {
            confs.setHeadTail(iChain,M,-1);
            confs.removeChain(iChain);
        }
        
    }
    else
    {
        
    }

    return accept;

}

moveHead::moveHead(int maxAdvanceLength_,int set) :
_maxReconstructedLength(maxAdvanceLength_+2),buffer((maxAdvanceLength_+2)*2,getDimensions()) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(set)
{

}


moveTail::moveTail(int maxAdvanceLength_,int set) :
_maxReconstructedLength(maxAdvanceLength_+2),buffer((maxAdvanceLength_+2)*2,getDimensions()) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(set)
{
    
}


bool moveHead::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();

    const auto & heads = confs.getGroups()[getSet()].heads;


    int iChain=heads[std::floor(uniformRealNumber(randG) * heads.size() )];
    

    int iHead=confs.getChain(iChain).head;

    
    int l= std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ; // distance from itime where the head is formed

    int t0= iHead - l;



    if ( t0 <= confs.getChain(iChain).tail )
    {
        return false;
    }
    

    
    Real deltaS=0;
    std::array<int,2> timeRange={t0  , iHead -1 };

    auto & data = confs.dataTensor();

    std::array<Real,3> difference;
    
    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);

    confs.copyDataToBuffer(buffer,{t0,iHead},iChain,0);

    // generates the head
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> headPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,t0);
    }

    Real mass = confs.getGroupByChain(iChain).mass;
    Real var=2*D*timeStep/ mass;

    confsSampler.sampleFreeParticlePosition(headPosition,startPosition,timeStep*l,randG,mass);

    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,iHead)=headPosition[d];
    }
    
    // perform levy reconstruction on l beads
    _levy.apply(confs,{t0,iHead},iChain,S,randG);


    bool accept =sPot.checkConstraints(confs,timeRange,iChain);


    if (accept)
    {
        // evaluates the action
    deltaS+=sPot.evaluate(confs,timeRange,iChain);


    auto propRatio = -deltaS;

    accept = sampler.acceptLog(propRatio,randG);

    }

    
    if ( accept)
    {
        
    }
    else
    {
        confs.copyDataFromBuffer(buffer,{t0,iHead},iChain,0);
    }

    return accept;

}

bool moveTail::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();

    const auto & tails = confs.getGroups()[getSet()].tails;

    if (not confs.isOpen(getSet()) )
    {
        throw std::runtime_error("Cannot move tail in closed sector.");
    }

    
    int iChain=tails[std::floor(uniformRealNumber(randG) * tails.size() )];
    int iChain2=confs.getChain(iChain).next;


    if (not confs.getChain(iChain).hasTail() )
    {
        throw invalidState("No tail to move.");
    }

    int iTail=confs.getChain(iChain).tail + 1;

    
    int l= std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) )  + 1; // distance from itime where the head is formed
    //int l=_maxReconstructedLength - 2;

    int M = confs.nBeads();

    int t1= iTail + l;
    int iChainStart=iChain;
    int tStart=t1;

    if (t1 > M )
    {
        if (iChain2 == -1)
        {
            return false;
        }

        iChainStart=iChain2;
        tStart=t1 - M;
    }


    if ( tStart > confs.getChain(iChainStart).head )
    {
        return false;
    }


    Real deltaS=0;
    std::array<int,2> timeRange={iTail  , std::min(t1,M) - 1 };
    std::array<int,2> timeRange2={0  , t1  - M -1 };

    auto & data = confs.dataTensor();

    std::array<Real,3> difference;
    
    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);
    deltaS-=sPot.evaluate(confs,timeRange2,iChain2);

    int len1= timeRange[1] - timeRange[0] + 2;
    confs.copyDataToBuffer(buffer,{timeRange[0],timeRange[1]+1},iChain,0);
    if (t1 > M)
    {
        confs.copyDataToBuffer(buffer,{timeRange2[0], timeRange2[1] + 1 },iChain2, len1    );

    }
    
    // generates the tail
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> tailPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChainStart,d,tStart);
    }


    auto confs_const = const_cast<pimcConfigurations_t*>(&confs);

    Real mass = (*confs_const).getGroupByChain(iChain).mass;

    Real var=2*D*timeStep/ mass;

    confsSampler.sampleFreeParticlePosition(tailPosition,startPosition,timeStep*l,randG,mass);

    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,iTail)=tailPosition[d];
    }

     for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,t1)=startPosition[d];
    }
 
    
    // perform levy reconstruction on l beads
    _levy.apply(confs,{iTail,t1},iChain,S,randG);
    if (t1>M)
    {
        confs.copyData( { M, t1} , iChain, 0,iChain2 );
    }
    
    bool accept = sPot.checkConstraints(confs,timeRange,iChain);
    accept=accept and sPot.checkConstraints(confs,timeRange2,iChain2);

    if (accept)
    {
        // evaluates the action
        deltaS+=sPot.evaluate(confs,timeRange,iChain);
        deltaS+=sPot.evaluate(confs,timeRange2,iChain2);

    
        
        auto propRatio = -deltaS;

        accept = sampler.acceptLog(propRatio,randG);

    }
    
    if ( accept)
    {

    }
    else
    {
        confs.copyDataFromBuffer(buffer,{timeRange[0],timeRange[1]+1},iChain,0);
        if (t1 > M)
        {
            confs.copyDataFromBuffer(buffer,{timeRange2[0], timeRange2[1] + 1 },iChain2, len1    );
        }

    }

    return accept;

}

void tableMoves::push_back( move * move_,Real weight,sector_t sector,const std::string & name)
{
    const auto & sets = move_->getSets();

    for ( auto set : sets)
    {
        if (set >=nRegisteredSets() )
        {
           resize(set+1);
        } 
        if (sector == sector_t::diagonal)
        {
            closedTabs[set].push_back(move_,weight,name);

        }
        else
        {
            openTabs[set].push_back(move_,weight,name);
        }
    }    

}

bool tableMoves::attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG)
{
    int iGroup= _chainSampler.sampleGroup(confs,randG);

    assert(iGroup < nRegisteredSets() );

    if (  confs.isOpen(iGroup) )
    {
        nOpenSectorMoves[iGroup]++;
        return openTabs[iGroup].attemptMove(confs,S,randG);
    }
    else
    {
        nClosedSectorMoves[iGroup]++;
        return closedTabs[iGroup].attemptMove(confs,S,randG);
    }

}

void tableMoves::resize(int nSets)
{
    closedTabs.resize(nSets);
    openTabs.resize(nSets);

    nOpenSectorMoves.resize(nSets);
    nClosedSectorMoves.resize(nSets);
};


std::ostream & tableMoves::operator>> (std::ostream & os)
{
    for (int iSet=0;iSet<nRegisteredSets() ;iSet++ )
    {
        os << "----Sector: " << iSet << "---------" << std::endl ; 
        os << "----------Open Sector" << "-----------" << std::endl;
        openTabs[iSet] >> os;
        os << "----------Closed Sector" << "-----------" << std::endl;
        closedTabs[iSet] >> os;
        os << "----------------" << std::endl;

        if ( (nOpenSectorMoves[iSet] + nClosedSectorMoves[iSet]) > 0 )
        {
            os << "Open Sector fraction: " << nOpenSectorMoves[iSet]/(nOpenSectorMoves[iSet] + nClosedSectorMoves[iSet]) << std::endl;
        }
    }

    return os;
}

   bool sampleSingleClosedHarmonicChain::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
    {

        auto & data=confs.dataTensor();
        int M=confs.nBeads();
        Real timeStep = S.getTimeStep();

        int tHead=M;


        Real sum=0;
        std::array<Real,3> x0;
        
        for(int t=0;t<=tHead;t++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                buffer(0,d,t)=data(iChain,d,t);
            }
        }

        

        for(int t=0;t<=tHead;t++)
        {
            Real prefactor= (t==0 or t==tHead) ? 0.5 : 1;
            for(int d=0;d<getDimensions();d++)
            {
                sum-=prefactor*data(iChain,d,t)*0.5*data(iChain,d,t);
            }
        }



        for(int d=0;d<getDimensions();d++)
            {
                x0[d]=data(iChain,d,0);

                data(iChain,d,0)=gauss(randG)*std::sqrt(var);
                data(iChain,d,tHead)=data(iChain,d,0);
            }


        _levy.apply(confs,{0,tHead},iChain,S,randG);

        // evaluate potential after
        
        for(int t=0;t<=tHead;t++)
        {
            Real prefactor= (t==0 or t==tHead) ? 0.5 : 1;
            for(int d=0;d<getDimensions();d++)
            {
                sum+=prefactor*data(iChain,d,t)*0.5*data(iChain,d,t);
            }
        }

        sum*=timeStep;

        Real logProbOld=0;
        Real logProbNew=0;

        for(int d=0;d<getDimensions();d++)
            {
                logProbOld+=- 0.5 * x0[d]*x0[d] / (var) ;
                logProbNew+=- 0.5 * data(iChain,d,0)*data(iChain,d,0) / ( var); 
            }


        Real prob = - logProbNew + logProbOld - sum;

        bool accept = sampler.acceptLog(prob,randG);

        if (accept)
        {
            
        }
        else
        {

            for(int t=0;t<=tHead;t++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    data(iChain,d,t)=buffer(0,d,t);
                }
            }

        }


        return accept;

    }


   bool sampleSingleOpenHarmonicChain::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
    {

        auto & data=confs.dataTensor();
        int M=confs.nBeads();
        Real timeStep = S.getTimeStep();

        int tHead=confs.getChain(iChain).head;
        int l = tHead;


        Real sum=0;
        std::array<Real,3> x0;
        
        for(int t=0;t<=tHead;t++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                buffer(0,d,t)=data(iChain,d,t);
            }
        }

        

        for(int t=0;t<=tHead;t++)
        {
            Real prefactor= (t==0 or t==tHead) ? 0.5 : 1;
            for(int d=0;d<getDimensions();d++)
            {
                sum-=prefactor*data(iChain,d,t)*0.5*data(iChain,d,t);
            }
        }



        for(int d=0;d<getDimensions();d++)
            {
                x0[d]=data(iChain,d,0);
                data(iChain,d,0)=gauss(randG)*std::sqrt(var);

                
                data(iChain,d,tHead)=gauss(randG)*std::sqrt(timeStep*l) + data(iChain,d,0);
            }


        _levy.apply(confs,{0,tHead},iChain,S,randG);

        // evaluate potential after
        
        for(int t=0;t<=tHead;t++)
        {
            Real prefactor= (t==0 or t==tHead) ? 0.5 : 1;
            for(int d=0;d<getDimensions();d++)
            {
                sum+=prefactor*data(iChain,d,t)*0.5*data(iChain,d,t);
            }
        }

        sum*=timeStep;

        Real logProbOld=0;
        Real logProbNew=0;

        for(int d=0;d<getDimensions();d++)
            {
                logProbOld+=- 0.5 * x0[d]*x0[d] / (var) ;
                logProbNew+=- 0.5 * data(iChain,d,0)*data(iChain,d,0) / ( var); 
            }


        Real prob = - logProbNew + logProbOld - sum;

        bool accept = sampler.acceptLog(prob,randG);

        if (accept)
        {
            
        }
        else
        {

            for(int t=0;t<=tHead;t++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    data(iChain,d,t)=buffer(0,d,t);
                }
            }

        }


        return accept;

    }


sampleSingleClosedHarmonicChain::sampleSingleClosedHarmonicChain(Real var_,int iChain_, int M) : 
var(var_),
iChain(iChain_),
_levy(2*M+1),
buffer(1,getDimensions(),M + 1) ,
singleSetMove(0),
gauss(0,1)
{

}




sampleSingleOpenHarmonicChain::sampleSingleOpenHarmonicChain(Real var_,int iChain_, int M) : 
var(var_),
iChain(iChain_),
_levy(2*M+1),
buffer(1,getDimensions(),M + 1) ,
singleSetMove(0),
gauss(0,1)
{

}


}

