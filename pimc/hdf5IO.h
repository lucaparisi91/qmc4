#include <string>
#include "hdf5.h"
#include <iostream>
#include <fstream>
#include <map>
#include <vector>


class dataIO 
{
    public:
    virtual  void annotate(const std::string & name, const int * values,const  size_t * dimensions, int rank, const std::string & label )=0;


    virtual  void annotate(const std::string & name, int value , const std::string & label )
    {
        size_t dims[1] {1};

        annotate(name, &value,dims,1,label);
    };


    virtual void readNote(const std::string & name,  int * values,  const std::string & label ) = 0 ;

    virtual void readNote(const std::string & name,  int & value ,  const std::string & label ) 
    {
        readNote(name,&value,label);
    }



    virtual void write( const std::vector<int> & data, const std::string & label)
    {
        size_t dims[1] = {data.size()};
        write(data.data(),label,dims,1);
    }

    virtual void write( const std::vector<double> & data, const std::string & label)
    {
        size_t dims[1] = {data.size()};
        write(data.data(),label,dims,1);
    }

    virtual void write( const double * data , const std::string & label, const size_t * dimensions, int rank )=0;


    virtual void write( const int * data , const std::string & label, const size_t * dimensions, int rank )=0;

    virtual void read(  int * data , const std::string & label )=0;
    virtual void read(  double * data , const std::string & label )=0;


    template<class T>
     T get(const std::string & label);

    virtual std::vector<size_t> getDimensions(const std::string & label)=0;


};




class hdf5IO : public dataIO
{
public:
    using dataIO::write;
    using dataIO::annotate;
    using dataIO::read;
    
    

    hdf5IO(std::string filename,std::ios_base::openmode mode);

    void open(std::string filename , std::ios_base::openmode mode);
    void close();


    void write( const double * data , const std::string & label, const size_t * dimensions, int rank ); 

    void write( const int * data , const std::string & label, const size_t * dimensions, int rank );

    


    void read(  double * data , const std::string & label );
    void read(  int * data , const std::string & label );

    void annotate(const std::string & name, const int * values,const  size_t * dimensions, int rank, const std::string & label );

    void readNote(const std::string & name,  int * values,  const std::string & label );

    std::vector<size_t> getDimensions(const std::string & label);


    ~hdf5IO();


    private:

    hid_t file;
    std::map<std::string,hid_t> hObjects;
    std::ios_base::openmode mode;



};