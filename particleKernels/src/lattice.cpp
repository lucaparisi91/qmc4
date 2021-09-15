#include "lattice.h"

using namespace particleKernels;

lattice::lattice(std::array<size_t,3>  shape, std::array<Real,3> lowerEdge, std::array<Real,3 > higherEdge,std::array<size_t,3> nGhosts) :
_shape(shape),
_nGhosts(nGhosts),
_lowerEdge(lowerEdge),
_higherEdge(higherEdge),
 dims(3),
 periodic(true),
 _nNeighboursPerCell(27)
{
    initialize();
};

void lattice::initialize()
{
    _size=1;
    _extendedSize=1;

    for (int d=0;d<dims;d++)
        {
            _lBox[d]=_higherEdge[d]-_lowerEdge[d];
            _lBoxInverse[d]=1./_lBox[d];
            _extendedShape[d]=_shape[d] + 2*_nGhosts[d];
            _size*=_shape[d];
            _extendedSize*=_extendedShape[d];
        }
    buildBC();
    buildIndexOffsets();
}

void lattice::buildBC()
{
    for(int d=0;d<dims;d++)
    {
        _wrap[d].resize(extendedSize() , 0 );
    }

    if ( dimensions() == 3)
    {
        for (int i=0;i<_extendedShape[0];i++)
            for (int j=0;j<_extendedShape[1];j++)
                for (int k=0;k<_extendedShape[2];k++)
                {
                    int currentIndex=index(i,j,k);
                    std::array<int,3> index3d{i,j,k};

                    for(int d=0;d<dimensions();d++)
                    {
                        if (index3d[d]<lowIndex(d))
                        {
                            _wrap[d][currentIndex]=-_lBox[d];
                        }
                        else if (index3d[d]>highIndex(d))
                        {
                            _wrap[d][currentIndex]=+_lBox[d];
                        }
                        else
                        {
                            _wrap[d][currentIndex]=0;
                        }
                        
                    }
                
                }
    }
}

void lattice::buildIndexOffsets()
{
    if (dims == 3)
    {
        int t=0;
        for (int k=0;k<=2;k++)
            for (int j=0;j<=2;j++)
                for (int i=0;i<=2;i++)
                {
                    _neighbourOffsets[t]=index(i,j,k) - index(1,1,1);
                    t++;
                }
    }

}

void lattice::checkNeighbourIndexing()
{
     if ( dimensions() == 3)
    {
        for (int i=lowIndex(0);i<=highIndex(0);i++)
            for (int j=lowIndex(1);j<=highIndex(1);j++)
                for (int k=lowIndex(2);k<=highIndex(2);k++)
                {

                    int iCell = index(i,j,k);

                    for (int ii=-1;ii<=1;ii++)
                        for (int jj=-1;jj<=1;jj++)
                            for (int kk=-1;kk<=1;kk++)
                            {
                                int iCell2=index(i+ii,j+jj,k+kk);
                                bool found=false;

                                for(int t=0;t<nCellsNeighbourhood();t++)
                                {
                                    int iCell2Check=getNeighbour(iCell,t);

                                    if (iCell2Check==iCell2)
                                    {
                                        found=true;
                                    }
                                }

                                assert(found==true);
                                

                            }
                }
    }
}