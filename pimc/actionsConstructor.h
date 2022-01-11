#include "actionConstructor.h"
namespace pimc
{

    class actionsConstructor
    {
        public:
        using geometry_t=geometryPBC_PIMC;

        std::vector<std::shared_ptr<pimc::action> > createActions(const json_t & js)
        {
            std::vector< std::shared_ptr<pimc::action> > actions;
            
            for (const auto & j : js)
            {
                std::string key = j["kind"].get<std::string>();

                actions.push_back( mapKeyToConstructor.at(key)->create(j) );
            }

            if ( actions.size() == 0 )
            {
                json_t jEmpty;
                actions.push_back( mapKeyToConstructor.at("nullPotential")->create( jEmpty));
            }


            return actions;

        }


        virtual void addConstructor(std::string key, const std::shared_ptr<actionConstructor> & constructor  )
        {
            mapKeyToConstructor[key]=constructor;
        }

        private:

        std::map<std::string, std::shared_ptr<actionConstructor>  > mapKeyToConstructor;


    };


};