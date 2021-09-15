#include "hdf5IO.h"
#include "hdf5.h"
#include <fstream>
#include "qmcExceptions.h"

hdf5IO::hdf5IO(std::string filename,std::ios_base::openmode mode_) : mode(mode_)
{

    open(filename,mode);
}

void hdf5IO::open(std::string filename,std::ios_base::openmode mode_)
{
    mode=mode_;
    
    if ( mode == std::ios::out)
    {
        file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    else if (mode ==  std::ios::in | std::ios::out)
    {
         file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);   
    }
    else
    {
        throw invalidInput("Mode not supported");
    }
    
}

void hdf5IO::close()
{
    
    H5Fclose(file);   
}


hdf5IO::~hdf5IO()
{
}

void hdf5IO::write( const double * data , const std::string & label, const size_t * dimensions , int rank)
{
    herr_t  status;
    hsize_t hDims[rank];

    for(int i=0;i<rank;i++)
    {
        hDims[i]=dimensions[i];
    }

    hid_t dataspace_id = H5Screate_simple(rank, hDims, NULL);

    hObjects[label + "_dataSpace"] = dataspace_id;


    hid_t dataset_id =
        H5Dcreate2(file, ("/" +label).c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    hObjects[label + "_dataSet"] = dataset_id;

    
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
}


void hdf5IO::write( const int * data , const std::string & label, const size_t * dimensions , int rank)
{
    herr_t  status;
    hsize_t hDims[rank];

    for(int i=0;i<rank;i++)
    {
        hDims[i]=dimensions[i];
    }

    

    hid_t dataspace_id = H5Screate_simple(rank, hDims, NULL);

    hObjects[label + "_dataSpace"] = dataspace_id;

    hid_t dataset_id =
        H5Dcreate2(file, ("/" +label).c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    hObjects[label + "_dataSet"] = dataset_id;

    
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
}


void hdf5IO::read( double * data , const std::string & label )
{
    if ( mode != ( std::ios::in | std::ios::out ) ) throw invalidState("Not open in read mode");

    hid_t dataset_id = H5Dopen2(file, ("/"+label).c_str(), H5P_DEFAULT);
    herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    status = H5Dclose(dataset_id);
}

void hdf5IO::read( int * data , const std::string & label )
{
    if ( mode != ( std::ios::in | std::ios::out ) ) throw invalidState("Not open in read mode");

    hid_t dataset_id = H5Dopen2(file, ("/"+label).c_str(), H5P_DEFAULT);

    herr_t status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    status = H5Dclose(dataset_id);
}

void hdf5IO::annotate(const std::string & name, const int * values, const size_t * dimensions,int rank,  const std::string & label )
{

   

    hsize_t hDims[rank];
    for(int i=0;i<rank;i++)
    {
        hDims[i]=dimensions[i];
    }
    hid_t dataspace_id = H5Screate_simple(rank, hDims, NULL);

    hid_t dataset_id = H5Dopen2(file,( "/" + label).c_str(), H5P_DEFAULT);


    hid_t attribute_id = H5Acreate2 (dataset_id, name.c_str(), H5T_NATIVE_INT,
                           dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    hObjects[label + "_attribute"]=attribute_id;

    herr_t status = H5Awrite(attribute_id, H5T_NATIVE_INT, values);    

    
    status = H5Aclose(attribute_id);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);  

}


void hdf5IO::readNote(const std::string & name,  int * values,  const std::string & label )
{

    hid_t attr = H5Aopen_by_name(file, ("/"+label).c_str(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT);

    herr_t status = H5Aread(attr, H5T_NATIVE_INT, values);
    status = H5Aclose(attr);



}

std::vector<size_t> hdf5IO::getDimensions(const std::string & label)
    {
        hid_t dataset_id = H5Dopen2(file, ("/"+label).c_str(), H5P_DEFAULT);
        hid_t dspace = H5Dget_space(dataset_id);

        int ndims = H5Sget_simple_extent_ndims(dspace);

        hsize_t dims[ndims];
        H5Sget_simple_extent_dims(dspace, dims, NULL);

        std::vector<size_t> outDims;

        outDims.resize(ndims);

        for(int i=0;i<ndims;i++)
        {
            outDims[i]=dims[i];
        }

        H5Dclose(dataset_id);
        H5Sclose(dspace);
        

        return outDims;

    }


    template<>
    std::vector<int> dataIO::get<std::vector<int> >(const std::string & label)
    {
        std::vector<int> data;
        auto dims = getDimensions(label);

        size_t totSize=0;
        
        for (auto dim : dims)
        {
            totSize+=dim;
        }
        data.resize(totSize);
        
        read(data.data(),label);


        return data;
    }

    template<>
    std::vector<double> dataIO::get<std::vector<double> >(const std::string & label)
    {
        std::vector<double> data;
        auto dims = getDimensions(label);

        size_t totSize=0;
        
        for (auto dim : dims)
        {
            totSize+=dim;
        }
        data.resize(totSize);
        
        read(data.data(),label);


        return data;
    }