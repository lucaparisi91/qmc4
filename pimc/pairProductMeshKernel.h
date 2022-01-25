
namespace pimc
{

    template<class greenFunction_t>
    class pairProductMeshKernel 
    {
        public:

        pairProductMeshKernel( std::shared_ptr<greenFunction_t> greenFunction_ , size_t N, size_t M) :
        _greenFunction(greenFunction_),
        disBuffer(N,M)
        {

        }

        void setGeometry(const geometry_t & geo) {_geo=geo;}

        Real evaluate(const Eigen::Tensor<Real,getDimensions()> & tn, const linkedCellParticles & particles,const std::array<int,2> & timeRange, int iParticle )
        {
            disBuffer.buildDistanceList(tn,particles,_geo,{timeRange[0],timeRange[1]-1},iParticle);

            const auto & buffer = disBuffer.getBufferDistances();
            const auto & bufferOff =  disBuffer.getBufferDistancesOffset() ;
            const auto & nDistances = disBuffer.getNDistances();



            Real sum=0;

            for (int t=timeRange[0];t<timeRange[1];t++)
            {
                for(int k=0;k<nDistances[k];k++)
                    {
                        sum+=_greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{bufferOff(k,0,t),bufferOff(k,1,t),bufferOff(k,2,t)}) ;
                    }     
            }

            return sum;

        }


        private:

        std::shared_ptr<greenFunction_t> _greenFunction;
        twoBodyDistancesBuffer disBuffer;
        geometry_t _geo;

    };

}

