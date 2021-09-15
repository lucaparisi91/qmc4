#include "geometry.h"
#include "tools.h"

real_t norm(real_t x,real_t y, real_t z)
{
	return std::sqrt(x*x + y*y + z*z);
}

geometryPBC::diff_t geometryPBC::differencesOneBody(const geometryPBC::particles_t & particleData, const std::array<real_t,3> & x) const
	{
	  int N = getN(particleData);
	  constexpr int D=getDimensions();

		diff_t diffs(N,D);

		for (int i=0;i<N;i++)
		{
			for(int d=0;d<D;d++)
				{
				diffs(i,d)=difference( particleData(i,d) - x[d] ,d);
				}
		}

		return diffs;
	}

distance_t norm( const difference_t & diffs)
{
  const int N = diffs.rows();
  constexpr int D=difference_t::ColsAtCompileTime;
	
  distance_t norms(N);
  for(int i=0;i<N;i++)
    {
      real_t tmp=0;
      for(int d=0;d<D;d++)
	{
	  tmp+= diffs(i,d)*diffs(i,d);	
	}
      norms(i)=std::sqrt(tmp);
    }

  return norms;
}


geometryPBC::diff_t geometryPBC::differencesTwoBody(const geometryPBC::particles_t & particleData) const 
	{
	         const int N = particleData.rows();
                 constexpr int D=difference_t::ColsAtCompileTime;

		diff_t diffs( (N*(N-1))/2,D);

		int k=0;
 		for(int i=0;i<N;i++)
 		{
 			for(int j=0;j<i;j++)
 			{
 				for (int d=0;d<D;d++)
 					diffs(k,d)=difference( particleData(i,d) - particleData(j,d) , d);
 				k++;
 			}	
		}
		return diffs;
	}


void geometryPBC::differencesTwoBody( diff_t &diffs, const geometryPBC::particles_t & particleData , int i) const 
	{
	    const int N = particleData.rows();

        constexpr int D=difference_t::ColsAtCompileTime;
		int k=(i*(i-1) )/2;

		// loop on a row of differences matrix
 		for(int j=0;j<i;j++)
 		{
 			for (int d=0;d<D;d++)
			 {
 				diffs(k,d)=difference( particleData(i,d) - particleData(j,d) , d);
 				k++;
 			}	
		 }

		 // loop on a column of the differences matrix
		 for (int j=i+1;j<N ;j++ )
		 {
			int k = (j*(j-1) )/2 + i;
			for (int d=0;d<D;d++)
			 {
			 diffs(k,d)=difference( particleData(j,d) - particleData(i,d) , d); // i and j inverted
			 }

		 }

	}



 geometryPBC::diff_t geometryPBC::differencesTwoBody(const geometryPBC::particles_t & data1,const geometryPBC::particles_t & data2) const 
	{
	  int N1 = data1.rows();
	  int N2 = data2.rows();
	  
	  constexpr int D=difference_t::ColsAtCompileTime;
	  

      
		diff_t diffs(N1*N2,D);

		int k=0;
 		for(int i=0;i<N1;i++)
 		{
 			for(int j=0;j<N2;j++)
 			{
 				for (int d=0;d<D;d++)
 					diffs(k,d)=difference( data1(i,d) - data2(j,d) , d);
 				k++;
 			}	
		}
		return diffs;
	}



geometry::diff_t geometry::differencesOneBody(const geometry::particles_t & data, std::initializer_list<real_t> l) const
	{
		std::array<real_t, 3> arr;
		int d=0;
		for (real_t x : l)
		 {
		 	arr[d]=x;d++;
		 }
		return differencesOneBody(data,arr);

	 }
