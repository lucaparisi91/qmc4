#include "pimcDriver.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "action.h"
#include "pimcConfigurations.h"
#include "moves.h"
#include "pimcObservablesFactory.h"
#include <filesystem>
#include "pimcPotentials.h"
#include "actionsConstructor.h"
#include "actionOneBodyConstructor.h"
#include "actionTwoBodyConstructor.h"
#include "propagators.h"
#include <csignal>

namespace fs = std::filesystem;





namespace pimc
{

    bool pimc_main_is_interrupted=false;

    void interceptSignal(int signal)
    {
        if (signal=2)
        {
            pimc_main_is_interrupted=true;
        }
    }




Real getTimeStep(json_t & j)
{
    int nBeads= j["nBeads"].get<Real>();
    Real beta= j["inverseTemperature"].get<Real>();

    return beta/nBeads;
}

pimcDriver::pimcDriver(const json_t & j_) : j(j_),
chemicalPotential({}),
currentEnsamble(ensamble_t::canonical)
{
    
    std::vector<Real> lBox;
    lBox=j["lBox"].get<std::vector<Real> >();

    if ( lBox.size() != getDimensions() )
    {
        throw invalidInput("This executable is compiled for " + std::to_string( getDimensions() ) + ". Input file implies " + std::to_string( lBox.size() ) + " dimensions" );
    }



    #if DIMENSIONS == 1
    geo=pimc::geometryPBC_PIMC(lBox[0],lBox[0],lBox[0]);
    #endif

    #if DIMENSIONS == 2
    geo=pimc::geometryPBC_PIMC(lBox[0],lBox[1],lBox[0]);
    #endif

    #if DIMENSIONS == 3
    geo=pimc::geometryPBC_PIMC(lBox[0],lBox[1],lBox[2]);
    #endif


    Real beta = j["inverseTemperature"].get<Real>();
    nParticles = j["particles"].get<std::vector<int> >();
    nMaxParticles = nParticles;

    nBeads= j["nBeads"].get<int>();

    timeStep = beta/nBeads;
    
    seed = j["seed"].get<int>();


    #if DIMENSIONS == 1
    pimc::geometryPBC_PIMC geo(lBox[0]);
    #endif

    #if DIMENSIONS == 2
    pimc::geometryPBC_PIMC geo(lBox[0],lBox[1]);
    #endif
    
    #if DIMENSIONS == 3
    pimc::geometryPBC_PIMC geo(lBox[0],lBox[1],lBox[2]);
    #endif

    //-------------- loading ensamble ---------------------
    

    if ( j.find("ensamble") != j.end() )
    {
        if ( j["ensamble"].get<std::string>() == "canonical" )
        {
            currentEnsamble= pimc::ensamble_t::canonical;
        }   
        else if ( j["ensamble"].get<std::string>() == "grandCanonical" )
        {
            if (j.find("chemicalPotential") == j.end() )
            {
                throw std::runtime_error("Chemical potential is required for the grandCanonical algorithm");

            }
            chemicalPotential = j["chemicalPotential"].get<std::vector<Real> >();
            currentEnsamble= pimc::ensamble_t::grandCanonical;

            if (j.find("maxParticles") != j.end() )
            {
                nMaxParticles=j["maxParticles"].get<std::vector<int> >();
            }
            else
            {
                for(int i=0;i<nParticles.size();i++)
                {
                    nMaxParticles[i]=10*nParticles[i];
                }
            }

        }
        else
        {
            throw std::runtime_error("Unkown ensamble");
        }
         
    }

    // ------------------ build moves ------------------------------


    pimc::moveConstructor pimcMoveConstructor(nMaxParticles,nBeads);
    pimcMoveConstructor.registerMove<pimc::levyMove>("levy");
    pimcMoveConstructor.registerMove<pimc::translateMove>("translate");
    pimcMoveConstructor.registerMove<pimc::openMove>("open");
    pimcMoveConstructor.registerMove<pimc::closeMove>("close");
    pimcMoveConstructor.registerMove<pimc::moveHead>("moveHead");
    pimcMoveConstructor.registerMove<pimc::moveTail>("moveTail");
    pimcMoveConstructor.registerMove<pimc::swapMove>("swap");

    if (currentEnsamble == ensamble_t::grandCanonical )
    {
        pimcMoveConstructor.registerMove<pimc::advanceHead>("advanceHead");
        pimcMoveConstructor.registerMove<pimc::recedeHead>("recedeHead");
        pimcMoveConstructor.registerMove<pimc::createWorm>("createWorm");
        pimcMoveConstructor.registerMove<pimc::deleteWorm>("deleteWorm");
    }

    tab = pimcMoveConstructor.createTable( j["movesTable"] );

    //------ loading information on I/O  --------------
    stepsPerBlock = j["stepsPerBlock"].get<int>();
    nBlocks = j["nBlocks"].get<int>();
    correlationSteps = j["correlationSteps"].get<int>();
    loadCheckPoint=false;
    doCheckPoint=false;

    if ( j.find("checkPointFile") != j.end() )
    {
        doCheckPoint=true;
        checkPointFile=j["checkPointFile"].get<std::string>();

        if (std::filesystem::exists(checkPointFile) )
        {
             loadCheckPoint=true;
        }
        ;
       
    }

}



void pimcDriver::run()
{
    // build action 
    std::shared_ptr<action> sT= std::make_shared<kineticAction>(timeStep, nParticles[0] , nBeads  , geo);

    

     #if DIMENSIONS == 3
     auto V = pimc::makePotentialFunctor(
         [](Real x,Real y , Real z) {return 0.5*(x*x + y*y + z*z) ;} ,
         [](Real x,Real y, Real z) {return x  ;},
         [](Real x,Real y,Real z) {return y ;},
         [](Real x,Real y,Real z) {return z ;}
         );
    #endif

    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return x ;} 
         );
    #endif

    int nChains = std::accumulate(nMaxParticles.begin(),nMaxParticles.end() , 0);

    auto sOneBodyC = std::make_shared<potentialActionOneBodyConstructor>();


    sOneBodyC->setTimeStep(timeStep);
    sOneBodyC->setGeometry( geo );
    sOneBodyC->setTimeStep(timeStep);

    sOneBodyC->registerPotential< isotropicHarmonicPotential>("harmonic");
    sOneBodyC->registerPotential<gaussianPotential>("gaussian");

    actionsConstructor sC;
    sC.addConstructor("oneBody",sOneBodyC);

    auto sTwoBodyCreator = std::make_shared< actionTwoBodyConstructor >() ;

    sTwoBodyCreator->setTimeStep(timeStep);
    sTwoBodyCreator->setGeometry(geo);
    sTwoBodyCreator->setNBeads(nBeads);
    sTwoBodyCreator->setNMaxParticles(nChains);


    
    
    

    sTwoBodyCreator->registerPotential<pimc::gaussianPotential>("gaussian");
    sTwoBodyCreator->registerPotential<pimc::isotropicHarmonicPotential>("harmonic");


    sTwoBodyCreator->registerGreenFunction<pimc::caoBernePropagator>("caoBerne");



    sC.addConstructor("twoBody",sTwoBodyCreator);


    //sC.registerPotential<isotropicHarmonicPotential>();
    //sC.registerPotential<gaussianPotential>();
    

    std::shared_ptr<action> sV=
    std::make_shared<sumAction>(sC.createActions(j["action"])); 

    S = pimc::firstOrderAction(sT, sV);
    
    randomGenerator_t randG(seed);

    std::vector<pimc::particleGroup> groups;
    int nStart=0;

    for ( int i=0;i<nParticles.size();i++)
    {
        pimc::particleGroup currentGroup{ nStart , nStart + nParticles[i] - 1, nStart +  nMaxParticles[i] - 1 , 1.0};
        groups.push_back(currentGroup);
        nStart+=nMaxParticles[i];     
    }

    pimc::pimcConfigurations configurations(nBeads, getDimensions() , groups );


    if (currentEnsamble == ensamble_t::grandCanonical)
    {
        configurations.setChemicalPotential(chemicalPotential);
    }

    // sets a random initial condition
    std::cout << "Generating initial configurations" << std::endl;

    
    configurations.setRandom( { geo.getLBox(0) ,geo.getLBox(1),geo.getLBox(2) } , randG );

    int iAttemptInitialCondition;

    while (not S.checkConstraints(configurations) )
    {
        configurations.setRandom( { geo.getLBox(0) ,geo.getLBox(1),geo.getLBox(2) } , randG );
        iAttemptInitialCondition++;

        if (iAttemptInitialCondition > 100000)
        {
            throw std::runtime_error("Max iteration reached in generating the initial condition");
        }
    }

    
    configurations.fillHeads();

    if (loadCheckPoint )
    {
        configurations=pimc::pimcConfigurations::loadHDF5(checkPointFile);
    }


    std::vector<std::shared_ptr<observable> > observables;
    
    pimcObservablesFactory obFactory(nBeads, nMaxParticles);

    obFactory.registerObservable<virialEnergyEstimator>("virialEnergy");
    obFactory.registerObservable<thermodynamicEnergyEstimator>("thermalEnergy");
    obFactory.registerObservable<particleNumberEstimator>("nParticles");
    obFactory.registerObservable<magnetizationSquaredEstimator>("magnetizationSquared");
    obFactory.registerObservable<magnetizationEstimator>("magnetization");

    obFactory.registerObservable<pairCorrelation>("pairCorrelation");

    if (j.find("observables") == j.end())
    {
        throw invalidInput("No abservables have been defined");
    }
    else
    {
       observables=obFactory.createObservables(j["observables"]) ;
    
    }
    
    auto eO=obFactory.getEnergyObservable();


    


   /*  pimc::levyMove freeMoves(5);

    Real delta=0.1;

    pimc::translateMove translMove(delta,(nBeads+1)*nParticles[0]);


   
     Real C = 1e-1;
    int l = 5;
    
    pimc::openMove openMove(C,l);
    pimc::closeMove closeMove(C,l);

    pimc::moveHead moveHeadMove(l);
    pimc::moveTail moveTailMove(l);

    pimc::swapMove swapMove(4,nParticles[0]);


     tab.push_back(& freeMoves,0.8,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(& freeMoves,0.8,pimc::sector_t::diagonal,"levy");

    tab.push_back(& translMove,0.2,pimc::sector_t::diagonal,"translate");
    tab.push_back(& translMove,0.2,pimc::sector_t::offDiagonal,"translate");

    tab.push_back(& openMove,0.2,pimc::sector_t::diagonal,"open");
    tab.push_back(& closeMove,0.2,pimc::sector_t::offDiagonal,"close");

    tab.push_back(& moveHeadMove,0.4,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(& moveTailMove,0.4,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(& swapMove,1.9,pimc::sector_t::offDiagonal,"swap");
 */
    
     // build initial  configuration


   


    pimc::thermodynamicEnergyEstimator energyEstimator;
    Real e=0;
    Real e2=0;

    std::ofstream f;

    std::ofstream ratioOut;

    ratioOut.open("nOpenClosed.dat");
    
    if ( ! fs::exists("configurations") ) 
    { 
        fs::create_directory("configurations"); // create src folder
    }

    
    //configurations.save("configurations/sample"+std::to_string(0));

    //configurations.save("configurations/sample"+std::to_string(0));
    int success = 0;

    // print the moves distribution

    int nGroups=configurations.getGroups().size();

    for (int i=0;i<nGroups;i++)
    {
        std::cout << "-----Set: " << i << ", open sector------"<<std::endl;

        std::cout << reportMoves(tab.getOpenSectorTable(i) );
        std::cout << "-----Set: " << i << ", closed sector------"<<std::endl;

        std::cout << reportMoves(tab.getClosedSectorTable(i)) ;

    }

    // start the simulation

    std::cout << "Start." << std::endl << std::flush;

    for (int i=0;(i< nBlocks) & (!pimc_main_is_interrupted); i++)
    {
        Real eStep=0;
        Real nClosed=0;
        Real nOpen=0;

        
        while ( nClosed < stepsPerBlock & (!pimc_main_is_interrupted) )
        {
            
            for (int j=0;(j<correlationSteps) & (!pimc_main_is_interrupted) ;j++)
            {
                bool accepted=tab.attemptMove(configurations, S, randG);

                if (accepted)
                {success+=1;}
            }
            
            if (!configurations.isOpen() )
            {
                
                for (auto & O : observables)
                {
                    O->accumulate(configurations,S);
                }
                nClosed+=1;
                
                
            }
            else
            {
               nOpen+=1;
                
            }
            
        }

        

        //std::cout << e << std::endl;
        if (eO != nullptr )
        {
            if (eO->weight() != 0 )
            {
                std::cout << "Energy: " << eO->average() << std::endl;
            }
        }

        for (auto & O : observables)
        {
            O->out(i);
            O->clear();
        }


        std::cout << "Acceptance ratio: " << success*1./((i+1)*stepsPerBlock*correlationSteps) << std::endl;

        
        tab >> std::cout;

        ratioOut << i << " " << nOpen<< " " << nClosed <<std::endl;
        

        configurations.saveHDF5("configurations/sample"+std::to_string(i+1) + ".hdf5" );

        if (doCheckPoint)
        {
            configurations.saveHDF5(checkPointFile);
        }


    }


    f.close();
    ratioOut.close();
    
    std::cout << "END." << std::endl;

}


}