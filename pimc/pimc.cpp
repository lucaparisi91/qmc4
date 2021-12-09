#include "pimcDriver.h"
#include <stdexcept>  

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        throw std::invalid_argument("Executable expects one argument: the input file name.");
    }

    std::string filename = argv[1];
    std::ifstream ifs(filename);
    pimc::json_t j = pimc::json_t::parse(ifs);
    
    pimc::pimcDriver driver(j);
    driver.run();

}
