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

    void levyReconstructor::apply (configurations_t & configurations, std::array<int,2> timeRange,int iChain ,const action & S,randomGenerator_t & randG)
    {
        /* Does not take into account time periodic boundary conditions
        Reconstruct between timeRange[1] -1 and timeRange[0] + 1
        timeRange[0] and timeRange[1] are used as boundary conditions
        */

        Real timeStep=S.getTimeStep();
        const auto & geo = S.getGeometry();

        int l = timeRange[1] - timeRange[0];

        auto & data = configurations.dataTensor() ;


        const auto & currentChain = configurations.getChain(iChain);
        int iChainNext=currentChain.next;


        //int lastBeadChain= timeRange[1] > configurations.nBeads() ? iChainNext : iChain ;
        //int lastBeadTime= timeRange[1] > configurations.nBeads() ? timeRanges[1][1] : timeRange[1] ;


 /*         if (lastBeadChain == -1 )
        {
            throw invalidState("Crossing a head in reconstruction.");
        }

        for (int d=0;d<getDimensions();d++)
            {
                data(iChain,d,timeRange[1])=data(lastBeadChain,d,lastBeadTime);
            }
  */

        std::array<Real,getDimensions()> xBottom;
        std::array<Real,getDimensions()> xTop;
        
        // performs the actual copy
        for (int t=0;t<l-1;t++)
        {

            for (int d=0;d<getDimensions();d++)
                {

                Real eta = gauss(randG);

                //xBottom[d]=geo.difference( data(iChain,d,t+timeRange[0]) - data(iChain,d,t+1+timeRange[0]),d) + data(iChain,d,t+1+timeRange[0]) ; 
                
                //xTop[d]=geo.difference( data(iChain,d,l+timeRange[0]) - data(iChain,d,t+1+timeRange[0]),d) + data(iChain,d,t+1+timeRange[0]) ; 

                xBottom[d]=data(iChain,d,t+timeRange[0]);
                xTop[d]=data(iChain,d,l+timeRange[0]);






                mean[d] = 
                (
                xBottom[d]*(l-t-1) 
                + xTop[d] )
                /( l - t) 
                ;

                Real variance = (l-t-1) * 1. /(l-t) *timeStep;

                data(iChain,d, timeRange[0] + t + 1 ) = mean[d] + eta *sqrt(variance) ;
                }
        }
        

        //configurationsNew.copyData(  {configurationsNew.nBeads() , timeRange[1]-1} , iChain , 0  , iChainNext);
        //configurationsNew.fillHead(iChain);

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

    auto timeRanges = splitPeriodicTimeSlice(timeRange,confs.nBeads());
    const auto & currentChain = confs.getChain(iChain);

    if (timeRange[0] <= currentChain.tail)
    {
        return false;
    }
    int iChainNext = currentChain.next;

    if (
        (timeRange[1] > currentChain.head and  (currentChain.next == -1 ) )
        or (   ( currentChain.next != -1 ) and
            ( confs.getChain(iChainNext).head < timeRanges[1][1]   )
    )
        )
    {
        return false; // do not accept if the time range is not included inside a chain (including time PBC)
    }
    
    auto & data = confs.dataTensor();

    
    auto sOld = S.evaluate(confs, timeRanges[0], iChain);
    sOld += S.evaluate(confs, timeRanges[1], iChainNext);


    // copy to internal buffer beads to move
     confs.copyDataToBuffer(buffer,{timeRanges[0][0] , timeRanges[0][1]+ 1},iChain);
     if (iChainNext!=-1)
     {
        confs.copyDataToBuffer(buffer,timeRanges[1],iChainNext, timeRanges[0][1]  - timeRanges[0][0] + 2);
     }


      if (timeRange[1] > confs.nBeads() )//copy the end bead in second chain to the first chain
    {
        confs.copyData( { timeRanges[1][1] , timeRanges[1][1] } , iChainNext, timeRange[1],iChain );

        for(int d=0;d<getDimensions();d++)
        {
            data(iChain,d,timeRange[1])-= data(iChainNext,d,0) - data(iChain,d,confs.nBeads() );
        }// ensures to reconstruct along a continuos path

        
    }
    _levy.apply(confs,timeRange,iChain,S,randG);
    if (iChainNext!=-1)
     {
         confs.copyData( { timeRanges[0][1]+1 , timeRange[1]-1 } , iChain, 0,iChainNext ); // time periodic boundary conditions

        for(int t=0;t<timeRanges[1][1];t++)
            for(int d=0;d<getDimensions();d++)
            {
                data(iChainNext,d,t)+=
                data(iChainNext,d,timeRanges[1][1]) - data(iChain,d,timeRange[1] );


                //std::cout << data(iChainNext,d,timeRanges[1][1]) - data(iChain,d,timeRange[1] )<<std::endl;
            }// ensures to reconstruct along a continuos path
     }




    auto  sNew= S.evaluate(confs,timeRanges[0], iChain) ;
    sNew+=S.evaluate(confs,timeRanges[1], iChainNext) ;

    const auto actionDifference = sNew - sOld;

    bool accepted = sampler.acceptLog(-actionDifference,randG);

    if (! accepted)
    {
        // copy back old beads
        confs.copyDataFromBuffer(buffer,{timeRanges[0][0],timeRanges[0][1]+1},iChain);
        if (iChainNext!=-1)
         {
        confs.copyDataFromBuffer(buffer,timeRanges[1],iChainNext, timeRanges[0][1]  - timeRanges[0][0] + 2);
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

    std::array<Real, 3> distance;
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


    std::array<Real, 3> distance;
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

    int iPartner=particleSampler.sample(randG) + group.iStart;
    

    if ( confs.getChain(iPartner).hasHead()   )
    {
        return false;
    }


    

    int iPartnerNext = confs.getChain(iPartner).next;
    
    assert(iPartnerNext != -1);

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


     if ( confs.getChain(iNewChainHead).tail >= tHead)
        {
            return false;
        }

    assert(iNewChainHead != -1);

    
     if ( confs.getChain(iNewChainHead).hasHead() or 
     confs.getChain(iNewChainHead).hasTail()
     )
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


    // performs levy reconstruction between the head and the bead
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
        if (tHead + l >=M )
        {
            confs.join(iChainHead,iPartner );
        }
        else
        {
            confs.copyData({tJoin+1,M},iPartner,iChainHead);
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
        TESTCHAIN(confs,iPartner,M,-1);
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

void sectorTableMoves::push_back(move * move_,Real weight,const std::string & name)
{
    sampler.accumulateWeight(weight);
    _moves.push_back(move_);
    _names.push_back(name);
    _nTrials.push_back(0);
    _nSuccess.push_back(0);

};

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


    //evaluate the action in the new configuration

     for (auto iCurrentChain : chainsInThePolimer)
     {
        deltaS+= Spot.evaluate(confs,timeRange,iCurrentChain);   
     }



    bool accept = sampler.acceptLog(-deltaS,randG);

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

    std::array<Real,3> difference;
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
                std::cout <<    data(iChain,d,iHead)-data(iChain,d,t0) << std::endl;
                //return false;
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
sigma2(sigma_*sigma_),C(C_)
{

}


bool createWorm::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    Real timeStep = S.getTimeStep();
    int M = confs.nBeads();
    const auto & geo = S.getGeometry();

    int l= std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) ) + 1 ;

    int lMax= _maxReconstructedLength -2;


    
    #ifndef NDEBUG
    int oldSize = confs.getGroups()[getSet()].size();
    #endif

    auto & data = confs.dataTensor();
    
    int tTail=std::floor(uniformRealNumber(randG) * M ) - 1 ;

    int t0 = tTail + 1;

    //std::cout << t0 << " -> " << t0+l << std::endl;

    int iChainTail=confs.pushChain(getSet());
    int iChainHead=iChainTail;
    int tHead=t0+l;
    confs.setHeadTail(iChainTail,std::min(t0 + l,M),tTail);

    if (t0 + l  > M)
    {
        iChainHead=confs.pushChain(getSet() );
        tHead=t0 + l - M; 
        confs.setHeadTail(iChainHead,tHead,-1);

    }

    // generates the tail according to a gaussian distribution
    std::array<Real,getDimensions() > x1;
    for(int d=0;d<getDimensions();d++)
    {
        x1[d]=gauss(randG);
    } 

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

    if (t0 + l > M)
    {
        confs.copyData( { M , t0 + l } , iChainTail, 0,iChainHead ); 
    }


    deltaS+=sPot.evaluate(confs,timeRange,iChainTail);
    deltaS+=sPot.evaluate(confs,timeRange2,iChainHead);

    auto propRatio = -deltaS + confs.getChemicalPotential()*l*timeStep - probabilityInitialPosition(geo,x1) + log(C*M*lMax);


    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        
        if (t0+l > M)
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
        
        if (t0+l > M)
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
sigma2(sigma_*sigma_),C(C_)
{

}

bool deleteWorm::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
    Real timeStep = S.getTimeStep();
    int M = confs.nBeads();
    const auto & geo = S.getGeometry();

    int lMax=_maxReconstructedLength -2;


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

    auto propRatio = -deltaS - confs.getChemicalPotential()*l*timeStep + probabilityInitialPosition(geo,x1) - log(C*M*lMax);

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
        Real coeff=N*C;

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

    int iChain = confsSampler.sampleChain(confs,getSet(),randG);
    

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
    confs.deleteBeads(timeRange,iChain);

    deltaS+=sPot.evaluate(confs,timeRange,iChain);

    // comput deltaS on the next chain ( if tTail overflows to the next chain)

    int t0_2=0;
    int t1_2=(tHead + l - M  );
    std::array<int,2> timeRange2={t0_2,t1_2-1};
    int iChainNext=confs.getChain(iChain).next;
    deltaS-=S.evaluate(confs,timeRange2,iChainNext);
    confs.deleteBeads(timeRange2,iChainNext);
    deltaS+=S.evaluate(confs,timeRange2,iChainNext);

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
        difference[d]=
               geo.difference( 
                data(iChainTail,d,tTail+1)-data(iChain,d,tHead),d
            );
    }


    auto propRatio = -deltaS - freeParticleLogProbability(difference,S.getTimeStep()*l,mass)  -  confs.getChemicalPotential()*l*timeStep + log( openCloseRatioCoefficient(N,M) );


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
                testChainsJoined(confs,iChainTail,iChainNext);
#endif 

        }
    }
    else
    {
        confs.createBeads(timeRange,iChain);
        confs.createBeads(timeRange2,iChainNext);
    }

    /* std::cout << "After open" << std::endl;
    confs >> std::cout;
    std::cout << std::endl;
 */
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

    std::array<Real,getDimensions()> difference; 
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=data(0,d,tHead+1)-data(0,d,tHead);

        if (order == 2)
        {
            deltaS-= 0.5 * 0.5 *( std::pow(data(0,d,tHead+1),2 ) + std::pow(data(0,d,tHead),2 )  )*timeStep;
        }
        if (order == 1)
        {
            deltaS-= 0.5*( std::pow(data(0,d,tHead),2 )  )*timeStep;
        }
        
    }
    


    auto propRatio = -deltaS - freeParticleLogProbability(difference,S.getTimeStep()*l,mass)   + log( C );


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
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]=data(0,d,0)-data(0,d,tHead);

        if (order == 2)
        {
            deltaS+= 0.5 * 0.5 *( std::pow(data(0,d,0),2 ) + std::pow(data(0,d,tHead),2 )  )*timeStep;
        }


        if (order == 1)
        {
            deltaS+= 0.5*( std::pow(data(0,d,tHead),2 )  )*timeStep;
        }
    }
    

    auto propRatio = -deltaS + freeParticleLogProbability(difference,S.getTimeStep()*l,mass)   - log( C );


    bool accept = sampler.acceptLog(propRatio,randG);


    if ( accept)
    {  

        confs.setHead(0,M);
        //confs.setHeadTail(0,M,-1);
        //confs.copyData({tHead+1,M}  , 1, 0  );
        confs.join(0,0);
        confs.fillHead(0);
        //confs.removeChain(1);
    }

    return accept;
};



closeMove::closeMove(Real C_ , int set,int maxReconstructionLength) : C(C_),_levy(2*(maxReconstructionLength+2)),_maxLength(maxReconstructionLength+2),buffer((maxReconstructionLength+2)*2,getDimensions()),gauss(0,1),uniformRealNumber(0,1),
setStartingBeadRandom(true),
setLengthCutRandom(true),
startingBead(-1),
lengthCut(-1),
singleSetMove(set)
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


    std::array<Real,3> difference;

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
        Real coeff=N*C;

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
bool closeMove::attemptGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG )
{
  /*   std::cout << "Before close" << std::endl;
    confs >> std::cout;
    std::cout << std::endl; */



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


    deltaS-=sPot.evaluate(confs,timeRange,iChainHead);

    t0=0;
    t1=tHead + l - M  ; 
    std::array<int,2> timeRange2={t0,t1-1 };


    deltaS-=sPot.evaluate(confs,timeRange2,iChainTail);
    configurations_t::copyData(confs, {tTail+1,tTail+1}, {iChainTail,iChainTail}, confs, tHead + l, iChainHead   );


    _levy.apply(confs,{tHead,tHead + l},iChainHead,S,randG);

    configurations_t::copyData(confs, {M,tHead + l}, {iChainHead,iChainHead}, confs, 0 , iChainTail   );


    std::array<Real,3> difference;


    confs.setHead(iChainHead,std::min(tHead+l,M));

    if (tHead + l >= M)
    {
        confs.setTail(iChainTail,-1);
              
        confs.join(iChainHead,iChainTail);
    }

    deltaS+=sPot.evaluate(confs,timeRange,iChainHead);
    deltaS+=sPot.evaluate(confs,timeRange2,iChainTail);

    for (int d=0;d<getDimensions();d++)
    {   
        difference[d]= geo.difference( -  data(iChainHead,d,tHead)  + data(iChainTail,d,tTail+1),d);
    }


    Real mass = confs.getGroupByChain(iChainHead).mass;

    auto propRatio = -deltaS + freeParticleLogProbability(difference,S.getTimeStep()*l,mass) -log( openCloseRatioCoefficient(N+1,M)   ) + confs.getChemicalPotential()*l*timeStep ;

    bool accept = sampler.acceptLog(propRatio,randG);
    
    if ( accept)
    {
        
        if (( tHead + l) < M)
        {
            confs.copyData({tTail + 1,M},iChainTail,iChainHead);
            int iChainNext=confs.getChain(iChainTail).next;
            confs.setHead(iChainHead,M);
            confs.join(iChainHead,iChainNext);
            #ifndef NDEBUG
            TESTCHAIN(confs,iChainHead,M,-1);
            testChainsJoined(confs,iChainHead,iChainNext);
            #endif
            confs.removeChain(iChainTail);
            
        }
        else
        {
#ifndef NDEBUG
        TESTCHAIN(confs,iChainHead,M,-1);
        TESTCHAIN(confs,iChainTail,M,-1);
        testChainsJoined(confs,iChainHead,iChainTail);
#endif
        }
    }
    else
    {
        confs.setHead(iChainHead,tHead);
        confs.setTail(iChainTail,tTail);
    }

   /*  std::cout << "After close" << std::endl;
    confs >> std::cout;
    std::cout << std::endl; */

    return accept;

};

advanceHead::advanceHead(int maxAdvanceLength_,int set) :
_maxReconstructedLength(maxAdvanceLength_+2) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(set)
{
    setRandomLength();
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

    deltaS-=sPot.evaluate(confs,timeRange,iChain);


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

    
    deltaS+=sPot.evaluate(confs,timeRange,iChain);

    int iChainNew=-1100;

    if (tHead + l > M)
    {
        // creates a new chain
        iChainNew=confs.pushChain(getSet());
        confs.setHeadTail(iChainNew,timeRange2[1] + 1,-1);
        confs.copyData( { M , tHead + l } , iChain, 0,iChainNew ); 
        deltaS+=sPot.evaluate(confs,timeRange2,iChainNew);

    }


    auto propRatio = -deltaS + confs.getChemicalPotential()*l*timeStep;


    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        if (tHead + l > M)
        {
            confs.join(iChain,iChainNew);
        }
    }
    else
    {
        confs.setHead(iChain,tHead);
        if (tHead + l > M)
        {
            confs.removeChain(iChainNew);
        }
    }

    return accept;
}

recedeHead::recedeHead(int maxAdvanceLength_,int set) :
_maxReconstructedLength(maxAdvanceLength_+2) , _levy((maxAdvanceLength_+2)*2),gauss(0,1),uniformRealNumber(0,1),singleSetMove(set)
{
    setRandomLength();
    
}


bool recedeHead::attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)
{
    Real timeStep = S.getTimeStep();
    const auto & heads = confs.getGroups()[getSet()].heads;

    int iChain=heads[std::floor(uniformRealNumber(randG) * heads.size() )];    

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
    deltaS+=sPot.evaluate(confs,timeRange,iChain);

    if (tHead - l < 0 )
    {
        confs.setHead(iChainTail,tHead - l + M);
        deltaS+=sPot.evaluate(confs,timeRange2,iChainTail);
    }

    auto propRatio = -deltaS -  confs.getChemicalPotential()*l*timeStep;


    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {
        if (tHead - l < 0 )
        {
            confs.removeChain(iChain);
        }
    }
    else
    {
        confs.setHead(iChain,tHead);
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

    // evaluates the action
    deltaS+=sPot.evaluate(confs,timeRange,iChain);


    auto propRatio = -deltaS;

    bool accept = sampler.acceptLog(propRatio,randG);

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

    
    int iChain=tails[std::floor(uniformRealNumber(randG) * tails.size() )];



    if (not confs.getChain(iChain).hasTail() )
    {
        throw invalidState("No tail to move.");
    }

    int iTail=confs.getChain(iChain).tail + 1;

    
    int l= std::floor( uniformRealNumber(randG) * (_maxReconstructedLength -2) )  + 1; // distance from itime where the head is formed
    //int l=_maxReconstructedLength - 2;


    int t1= iTail + l;

    if (t1 >= confs.getChain(iChain).head )
    {
        return false;
    }
    

    Real deltaS=0;
    std::array<int,2> timeRange={iTail  , t1 - 1 };


    auto & data = confs.dataTensor();

    std::array<Real,3> difference;
    
    auto & sPot = S.getPotentialAction();

    deltaS-=sPot.evaluate(confs,timeRange,iChain);

    confs.copyDataToBuffer(buffer,{iTail,t1},iChain,0);

    // generates the tail
    Real distanceSquared=0;
    const auto & geo = S.getGeometry();
    std::array<Real,getDimensions()> tailPosition;
    std::array<Real,getDimensions()> startPosition;

    for (int d=0;d<getDimensions();d++)
    {
        startPosition[d]=data(iChain,d,t1);
    }

    auto confs_const = const_cast<pimcConfigurations_t*>(&confs);

    Real mass = (*confs_const).getGroupByChain(iChain).mass;

    Real var=2*D*timeStep/ mass;

    confsSampler.sampleFreeParticlePosition(tailPosition,startPosition,timeStep*l,randG,mass);

    for (int d=0;d<getDimensions();d++)
    {
        data(iChain,d,iTail)=tailPosition[d];
    }

    
    // perform levy reconstruction on l beads
    _levy.apply(confs,{iTail,t1},iChain,S,randG);

    // evaluates the action
    deltaS+=sPot.evaluate(confs,timeRange,iChain);

    auto propRatio = -deltaS;

    bool accept = sampler.acceptLog(propRatio,randG);

    if ( accept)
    {

    }
    else
    {
        confs.copyDataFromBuffer(buffer,{iTail,t1},iChain,0);
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

}


