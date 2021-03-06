#include "pimcObservables.h"
#include <fstream>
namespace pimc
{

Real thermodynamicEnergyEstimator::operator()(configurations_t & confs, firstOrderAction & S)
{
    auto & geo = S.getGeometry();

    auto & kA = S.getKineticAction();
    auto & potA = S.getPotentialAction();

    auto sA=kA.evaluate(confs);
    auto sV=potA.evaluateTimeDerivative(confs);

    auto beta = confs.nBeads() * kA.getTimeStep(); 
    sA/=beta;
    sV/=confs.nBeads();

    Real e= sV - sA +  getDimensions()/(2.*kA.getTimeStep())*confs.nParticles();

    return e;
}

Real virialEnergyEstimator::operator()(configurations_t & confs, firstOrderAction & S)
{
    auto & geo = S.getGeometry();
    auto & Spot = S.getPotentialAction();
    Real e=0;
    Real e2=0 , e3 = 0 , e4=0;

    int N=0;
    buffer.setConstant(0.);

    Spot.addGradient(confs,buffer);
    
    for ( const auto & group : confs.getGroups() )
    {
        N+=group.iEnd - group.iStart + 1;
        int iStart=group.iStart;
        int iEnd=group.iEnd;

        

      for(int i= iStart;i<=iEnd ; i++)
        {
            int iPrev=confs.getChain(i).prev;
            assert(iPrev >=0 );
            for(int d=0;d<getDimensions();d++)
            {
                buffer(i,d,0)+=buffer(iPrev,d,confs.nBeads() );
            }
            
        }  
    }

    const auto & data = confs.dataTensor();     

    // compute rC
    rC.setConstant(0);
    for (const auto & group : confs.getGroups() )
    {

        for (int alpha=0; alpha< confs.nBeads();alpha++)
        {
            // sum on time slices greater than alpha, same chain
            for (int t=alpha;t<confs.nBeads();t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,alpha)+=0.5*data(i,d,t)/ confs.nBeads();
                    }
                }
                
            }
            // sum on time slices greater than alpha, next chain
            for (int t=0;t<alpha;t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    int ii=confs.getChain(i).next;

                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,alpha)+=0.5*
                        (data(ii,d,t) + data(i,d,confs.nBeads()) - data(ii,d,0) )/ confs.nBeads();
                    }
                }
            }

            // sum on time slices lesser than alpha, prev chain
            for (int t=alpha + 1;t< confs.nBeads();t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    int ii=confs.getChain(i).prev;

                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,alpha)+=
                        0.5*(
                            
                            data(ii,d,t) + data(i,d,0) - data(ii,d,confs.nBeads())
                        
                        )/ confs.nBeads();
                    }
                }
            }

            // sum on time slices lesser then alpha , same chain
            for (int t=0;t<=alpha;t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        rC(i,d,alpha)+=0.5*data(i,d,t)/ confs.nBeads();
                    }
                }
                
            }

           /*  std::cout << alpha<<" ";
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        std::cout << rC(i,d,alpha) << " ";
                    }
                }
            std::cout <<std::endl;
             */


        }

        // second term in the virial estimator
        // Does not give any contribution for a classical system

        
        for (int t=1;t<confs.nBeads();t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                   
                    int inext = confs.getChain(i).next;

                     assert( (inext>= group.iStart) and (inext<= group.iEnd) );
                     
                    for(int d=0;d<getDimensions();d++)
                    {
                        e2+= 
                        /*
                        ( geo.difference( data(inext,d,t ) - data(i,d,t) ,d )  ) *
                        ( geo.difference( data(inext,d,t ) - data(i,d,t) ,d )  );
                        
                        geo.difference( -data(i,d,t) + data(inext,d,t ) ,d )*
                        geo.difference( -data(i,d,t+1) + data(i,d,t ) ,d );
                        */
                        (-data(i,d,t) + data(inext,d,t ) + data(i,d,confs.nBeads())  - data(inext,d,0)   )*
                        (  data(inext,d,t-1)   - data(inext,d,t)   );


                    }
                    
                    
                }
            }
        //std::cout << e2 << std::endl;
        {
            int t=0;
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    int inext = confs.getChain(i).next;

                    for(int d=0;d<getDimensions();d++)
                    {
                        e2+= (-data(i,d,t) + data(i,d,confs.nBeads() )   )*( data(i,d,confs.nBeads()-1) - data(i,d,confs.nBeads()  ) );
                        
                    }
                    
                    
                }
            } 


        
        // third term in the virial estimator
        for (int t=0;t<confs.nBeads();t++)
            {
                for (int i = group.iStart ; i<=group.iEnd ; i++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        
                        auto delta=( data(i,d,t) - rC(i,d,t) ) ;
                        //auto delta=geo.difference( data(i,d,t) - rC(i,d,t) ,d) ;

                        
                        e3+= delta *
                        buffer(i,d,t);
                    }

                }
             
            }
    }

    Real beta = S.getTimeStep() * confs.nBeads();
    e4= S.getPotentialAction().evaluateTimeDerivative(confs);
    e4/=( confs.nBeads() );

    e3/=(2 * beta );
    e2/=(2*beta*beta);
    
    // classical gas free contribution
    Real e1 = confs.nParticles()*   getDimensions()/(2*beta);
    return e1 + e2  + e3 + e4;

}


pairCorrelation::pairCorrelation(int setA_, int setB_) :
setA(setA_),setB(setB_)
{
    
}

angleEstimator::angleEstimator(int setA_, int setB_) :
setA(setA_),setB(setB_)
{
    
}


Real pairCorrelation::getNormalizationFactor(const configurations_t & configurations, const firstOrderAction & S , const accumulator_t & acc) const
{
    auto geo = S.getGeometry();
    auto dx = acc.stepSize();

    const auto & groupA = configurations.getGroups()[setA];
    const auto & groupB = configurations.getGroups()[setB];
    
    auto  NA = groupA.iEnd - groupA.iStart + 1;
    auto  NB = groupB.iEnd - groupB.iStart + 1;


#if DIMENSIONS == 3
    //auto R = 0.5*std::sqrt( geo.getLBox(0)*geo.getLBox(0) + geo.getLBox(1)*geo.getLBox(1) + geo.getLBox(2)*geo.getLBox(2)  );
    //auto R = 0.5*geo.getLBox(0);
    //auto V = 4./3.*M_PI*R*R*R;
    auto V = geo.getLBox(0)*geo.getLBox(1)*geo.getLBox(2);
    auto _normalizationFactor=V/(dx*4*M_PI*NA*NB);
#endif

#if DIMENSIONS == 1
  auto V = geo.getLBox(0)*geo.getLBox(0);
  auto _normalizationFactor=std::pow(V,1)/(dx*NA*NB*2);   
#endif
  
  _normalizationFactor = (setA == setB) ? 2*_normalizationFactor : _normalizationFactor;
  
  return _normalizationFactor;
}


void pairCorrelation::operator()(configurations_t & configurations, firstOrderAction & S,  accumulator_t & histAcc)
{

    if (setA == setB)
    {
        accumulateUnDistinguishable(configurations,S,histAcc);
    }
    else
    {
        accumulateDistinguishable(configurations,S,histAcc);
    }
}

void angleEstimator::operator()(configurations_t & confs, firstOrderAction & S,  accumulator_t & acc)
{
    if (setA == setB)
    {
        const auto & groupA = confs.getGroups()[setA];
        const auto & data = confs.dataTensor();
        const auto & geo = S.getGeometry();

        auto N = groupA.iEnd - groupA.iStart + 1;
        
        for(int t=0;t<confs.nBeads();t++)
        {
            for(int i=groupA.iStart;i<=groupA.iEnd;i++)
                for(int j=groupA.iStart;j<i;j++)
                {
                    std::array<Real,DIMENSIONS> r,rNext;
                    Real cosTeta=0;
                    Real r2=0 , r2Next=0;

                    for (int d=0;d<getDimensions();d++)
                        {
                            r[d]=geo.difference( data(i,d,t) - data(j,d,t)  ,d );
                            rNext[d]=geo.difference( data(i,d,t+1) - data(j,d,t+1)  ,d );

                            r2+=r[d]*r[d];
                            r2Next+=rNext[d]*rNext[d];
                            cosTeta+=(r[d]*rNext[d] );
                        }
                    
                    r2= std::sqrt(r2);
                    r2Next= std::sqrt(r2Next);
                    if (r2*r2Next != 0 )
                    {
                        cosTeta/=(r2*r2Next);
                        //std::cout << cosTeta << std::endl;

                        if( 
                            ( cosTeta > acc.minx() ) and  
                             ( cosTeta < acc.maxx() )
                            )
                            {
                                acc.accumulate(1,cosTeta);
                            }
                    }
                    
                }
        }

        acc.weight()+=confs.nBeads()*N*(N-1)/2.;

    }
    else
    {
        throw std::runtime_error("angle estimator not supported for different sets A and B");
    }
}


void oneBodyEstimator::operator()(configurations_t & confs, firstOrderAction & S,  accumulator_t & acc)
{    
    const auto & groupA = confs.getGroups()[_set];
    const auto & geo = S.getGeometry();

    if (confs.isOpen(_set) )
    {
    auto iChainHead= groupA.heads[0];
    auto iChainTail= groupA.tails[0];

    auto tHead=confs.getChain(iChainHead).head;
    auto tTail=confs.getChain(iChainTail).tail ;
    const auto & data = confs.dataTensor();

  
        if (tHead == tTail + 1)
        {
            Real r2=0;
            for (int d=0;d<getDimensions();d++)
            {
                auto tmp=geo.difference( data(iChainHead,d,tHead) - data(iChainTail,d,tTail + 1)  ,d );
                r2+=tmp*tmp;
            }

            Real r= std::sqrt(r2);
            if( 
                    ( r > acc.minx() ) and  
                        ( r < acc.maxx() )
                )
                {
                    acc.accumulate( 1/(4*M_PI * r2*acc.stepSize()),r);
                } 
            acc.weight()+=1;
        }

    }


}



void pairCorrelation::accumulateUnDistinguishable(configurations_t & confs, firstOrderAction & S,  accumulator_t & acc)
{    
    const auto & groupA = confs.getGroups()[setA];
    
    const auto & data = confs.dataTensor();
    const auto & geo = S.getGeometry();

    auto norm = getNormalizationFactor(confs,S,acc);

    for(int t=0;t<confs.nBeads();t++)
    {
        for(int i=groupA.iStart;i<=groupA.iEnd;i++)
            for(int j=groupA.iStart;j<i;j++)
            {
                Real r2=0;
                for (int d=0;d<getDimensions();d++)
                    {
                        auto tmp=geo.difference( data(i,d,t) - data(j,d,t)  ,d );
                        r2+=tmp*tmp;
                    }
                Real r= std::sqrt(r2);

                if( 
                    ( r > acc.minx() ) and  
                    ( r < acc.maxx() )
                  )
                {
                    acc.accumulate(norm,r);
                } 

            }
    }

    acc.weight()+=confs.nBeads();
}

void pairCorrelation::accumulateDistinguishable(configurations_t & confs, firstOrderAction & S,  accumulator_t & acc)
{    
    const auto & groupA = confs.getGroups()[setA];
    const auto & groupB = confs.getGroups()[setB];

    const auto & data = confs.dataTensor();
    const auto & geo = S.getGeometry();

    auto norm = getNormalizationFactor(confs,S,acc);


    for(int t=0;t<confs.nBeads();t++)
    {
        for(int i=groupA.iStart;i<=groupA.iEnd;i++)
            for(int j=groupB.iStart;j<=groupB.iEnd;j++)
            {
                Real r2=0;
                for (int d=0;d<getDimensions();d++)
                    {
                        auto tmp=geo.difference( data(i,d,t) - data(j,d,t)  ,d );
                        r2+=tmp*tmp;
                    }
                Real r= std::sqrt(r2);

                if( 
                    ( r > acc.minx() ) and  
                    ( r < acc.maxx() )
                  )
                {
                    acc.accumulate(norm,r);
                } 

            }
    }
    
    acc.weight()+=confs.nBeads();
}


lengthEstimator::lengthEstimator( int l) : lengthEstimator::lengthEstimator()
{
    _setMaxLength=true;
    _l=l;

}

lengthEstimator::lengthEstimator( ) : 
_l(0),_startFromHeadOrTail(false),_t0(-1),_iChain(0),_setMaxLength(false) 
{

}

void lengthEstimator::setStartFromHeadOrTail(int iGroup)
{
    _iGroup=iGroup;
    _startFromHeadOrTail=true;
}

int lengthEstimator::initialChain( const configurations_t & configurations) const
    {
        if (not _startFromHeadOrTail)
        {
            return _iChain;
        }
        else
        {
            if (_l >= 0)
            {
                return configurations.getGroups()[_iGroup ].tails[0];
            }
            else
            {
                return configurations.getGroups()[_iGroup ].heads[0];
            }

        }
    } 

int lengthEstimator::initialTime( const configurations_t & configurations) const
    {
        if (not _startFromHeadOrTail)
        {
            return _t0;
        }
        else
        {
            int iChain=initialChain(configurations);
            if (_l >= 0)
            {
                return configurations.getChain(iChain).tail + 1;
            }
            else
            {
                return configurations.getChain(iChain).head;
            }

        }
    } 


Real lengthEstimator::operator()(configurations_t & configurations, firstOrderAction & S) {

    const auto & data = configurations.dataTensor();
    auto M=configurations.nBeads();

    int currentChain=initialChain(configurations);
    int l= _l;
    int t0= initialTime(configurations);

    if ( not _setMaxLength)
    {
        l=configurations.nChains()* configurations.nBeads();
    }




    std::array<int,2> timeRange;    
    

    Real sum=0;
    while(currentChain !=-1 )
    {
        

        auto tHead=configurations.getChain( currentChain ).head;
        auto tTail=configurations.getChain( currentChain ).tail ;


        if (l>0)
        {
            timeRange[0]=std::max(t0,tTail+1);
            timeRange[1]=std::min( t0+l,tHead );    
        }

        {
            timeRange[1]=std::min( t0 , tHead);
            timeRange[0]=std::max(t0+l,tTail + 1);
        }


        for(int t=timeRange[0];t<timeRange[1];t++ )
        {
            for(int d=0;d<getDimensions(); d++)
            {
                Real tmp=data(currentChain,d,t+1) - data(currentChain,d,t);
                sum+=tmp*tmp;
            }
        }


        if( ( timeRange[0]+ l) > M )
        {
            currentChain=configurations.getChain(currentChain).next;
            l-=( timeRange[1] - timeRange[0]  );
            t0=0;
        }
        else if ( timeRange[1] + l < 0 )
        {
            l+=timeRange[1]-timeRange[0];
            t0=M;
            currentChain=configurations.getChain(currentChain).prev;            
        } 
        else
        {
            currentChain=-1;
        }

        if ( ( not _setMaxLength) and ( currentChain == _iChain  ) )
        {
            currentChain=-1;
        }

    }

    return sum;
}

Real closedLengthEstimator::operator()(configurations_t & configurations, firstOrderAction & S) {

    const auto & group = configurations.getGroups()[_set];

    const auto & data=configurations.dataTensor();

    std::array<int,2> timeRange { 0 , configurations.nBeads() - 1};

    Real sum=0;

    for(int t=timeRange[0];t<=timeRange[1];t++ )
    {
        for(int i=group.iStart; i<=group.iEnd;i++ )
        {
            for(int d=0;d<getDimensions(); d++)
            {
                Real tmp=data(i,d,t+1) - data(i,d,t);
                sum+=tmp*tmp;
            }
        }
    }


        

    return sum;
}


Real nBeadsInWormEstimator::operator()(configurations_t & configurations, firstOrderAction & S) {

    auto currentChain=configurations.getGroups()[_iGroup].tails[0];
    int l=0;

    while(currentChain !=-1 )
    {
        

        auto tHead=configurations.getChain( currentChain ).head;
        auto tTail=configurations.getChain( currentChain ).tail ;
        l+=tHead  - (tTail + 1);

        currentChain=configurations.getChain(currentChain).next;
        
    }
    
    return l;
}


Real nConnectedChainsEstimator::operator()(configurations_t & configurations, firstOrderAction & S) {


    return chainCounter.count(configurations,_set);
}

void openRatio::accumulate(configurations_t & configurations)
{
    const auto & groups = configurations.getGroups();

    int j=0;
    for(int i=0;i<groups.size() ;i++)
    {
        j+= (int) (groups[i].isOpen()) *  std::pow(2,i); 
    }
    assert( j< nOpen.size());
    nOpen[j]+=1;
    n+=1;
}

void openRatio::out(int t)
{
    f << t << "\t";

    if (n==0)
    {
        return ;
    }

    for(int i=0;i<nOpen.size();i++)
    {
        f << "\t" << average( i ) ;
    }

    f<<std::endl;
}

openRatio::openRatio( int nSets ) :
n(0) , filename("ratio.dat")
{
    nOpen.resize(  1<< nSets , 0 );
    assert(nSets>=1);

    if (fs::exists(fs::path(filename)) )
    {
        f.open(filename,std::fstream::app);
    }
    else
    {
        f.open(filename);
    
        f << "iteration"<< "\t";

        if (nSets==1)
        {
            f << "cRatio" << "\t" << "oRatio" << std::endl;
        }
        else if (nSets==2) 
        {
            f << "ccRatio" << "\t" << "ocRatio" << "\t" << "coRatio" << "\t" << "ooRatio" << std::endl;
        }
        else{
            throw std::runtime_error("nSets should be 1 or 2");
        }
    }
        
}

openRatio::~openRatio()
{
    f.close();
}


magnetizationDistribution::magnetizationDistribution( const json_t & j)
{
    _mMin=j["min"].get<int>();
    _mMax=j["max"].get<int>();
    setA=j["sets"][0].get<int>();
    setB=j["sets"][1].get<int>();

    _Ms.resize( _mMax - _mMin + 1,0);
    _label=j["label"].get<std::string>();
    std::string filename = _label + ".dat";

    if (fs::exists(fs::path(filename)))
            {
                f.open(filename,std::fstream::app);
            }
            else
            {
                f.open(filename);
                f << "time" << "\t" << "M" << "\t" << "P" << std::endl;
            }
    
    n=0;

    if (_mMin < 0 )
    {
        recordAbsoluteValue=false;
    }
    else
    {
        recordAbsoluteValue=true;
    }


}


void magnetizationDistribution::accumulate(configurations_t & configurations, firstOrderAction & S)
{
    int M=configurations.nParticles(setA) - configurations.nParticles(setB) ;

    if ( recordAbsoluteValue)
    {
        M=std::abs(M);
    }


    if ( (M>=_mMin) and (M<=_mMax) )
    {
        _Ms[ M - _mMin ]+=1;
    }
    n+=1;
}

void magnetizationDistribution::out( size_t t)
{
    if (n !=0)
    {
        for(int i=0;i<_Ms.size();i++)
        {
            f << t << "\t" << _mMin + i << "\t" << _Ms[i]/n << std::endl;
        }
    }

    f << std::flush;
}

void magnetizationDistribution::clear( )
{
    
    n=0;
    std::fill(_Ms.begin(),_Ms.end(),0);
}

magnetizationDistribution::~magnetizationDistribution( )
{
    f.close();   
}


void thermodynamicEnergyEstimatorMagnetization::operator()(configurations_t & configurations, firstOrderAction & S, accumulator_t & acc)
{
    int M=std::abs(configurations.nParticles(setA) - configurations.nParticles(setB) );
    acc.accumulate( energyEst(configurations,S) , M  );

}

void virialEnergyEstimatorMagnetization::operator()(configurations_t & configurations, firstOrderAction & S, accumulator_t & acc)
{
    int M=std::abs(configurations.nParticles(setA) - configurations.nParticles(setB) );
    acc.accumulate( energyEst(configurations,S) , M  );
}

void superfluidFractionEstimator::operator()(configurations_t & configurations, firstOrderAction & S,accumulator_t & acc)
    {
        const auto & geo = S.getGeometry();

        const auto & groups= configurations.getGroups();
        const auto & data = configurations.dataTensor();
        auto nBeads = configurations.nBeads();

        auto beta= S.getTimeStep() * nBeads;

        for(int t=0;t<nBeads;t++)
            {
            std::vector< std::array<Real,DIMENSIONS>  > cm;
            cm.resize( groups.size() , {0,0,0}  );
            
            for (int iGroup=0;iGroup<groups.size();iGroup++)
                {
                const auto & group = groups[iGroup];
                for(int i=group.iStart;i<=group.iEnd;i++)
                    {
                        auto iNext=configurations.getChain(i).next;
                        for(int d=0;d<getDimensions();d++)
                        {
                            Real diff=0;
                            diff+=data(i,d,t) - data(i,d,nBeads);
                            diff+= data(iNext,d,0) - data(iNext,d,t);
                            cm[iGroup][d]+=diff;
                        }
                    }
                }
                int k=0;
                for(int iGroup=0;iGroup < groups.size();iGroup++)
                {
                    for(int jGroup=0;jGroup<=iGroup; jGroup++)
                    {

                        Real rho=0;
                        for(int d=0;d<DIMENSIONS;d++)
                        {
                            rho+=cm[iGroup][d] * cm[jGroup][d];
                        }

                        const auto & groupA = groups[iGroup];
                        const auto & groupB = groups[jGroup];

                        auto NA = groupA.iEnd - groupA.iStart + 1;
                        auto NB = groupB.iEnd - groupB.iStart + 1;

                        acc.accumulate(rho/( sqrt(NA*NB) * beta * DIMENSIONS), k );
                        k++;
                    }
                }

            }
    }

bool observable::isValidSector(const configurations_t & confs)
{
    auto nGroups = confs.getGroups().size();

    bool valid=true;

    for(int i=0;i<nGroups;i++)
    {
        valid = valid and (not confs.isOpen(i));
    }

    return valid;
}

bool oneBodyObservable::isValidSector(const configurations_t & confs)
{
    auto nGroups = confs.getGroups().size();
    bool valid=true;

    for(int i=0;i<nGroups;i++)
    {
        if (i== _set)
        {
            valid = valid and ( confs.isOpen(i));
        }
        else
        {
            valid = valid and (not confs.isOpen(i));
        }

    }



    return valid;
}
    






};