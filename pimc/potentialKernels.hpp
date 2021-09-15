
namespace pimc
{

    namespace kernels
    {
        namespace potentials
        {   

            template<class functor_t>
            Real reduceOneBodyPotentials2OrderIsotropic(const functor_t & V,const Eigen::Tensor<Real,3> & data, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange)
            {
                Real sum=0;

                for(int t=timeRange[0]+1;t<=timeRange[1] ; t++ )
                    {
                    for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                        {
                            Real r2=0;
                            for (int d=0;d<getDimensions();d++)
                            {
                            r2+=data(i,d,t)*data(i,d,t); 
                            }

                            sum+=V( std::sqrt(r2) );
                        }
                    }
                
                    {
                        int t = timeRange[0];
                        for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                            {
                                Real r2=0;
                                for (int d=0;d<getDimensions();d++)
                                {
                                r2+=data(i,d,t)*data(i,d,t); 
                                }

                                sum+=0.5*V( std::sqrt(r2) );
                            }
                    }

                    {
                        int t = timeRange[1] + 1;
                        for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                            {
                                Real r2=0;
                                for (int d=0;d<getDimensions();d++)
                                {
                                r2+=data(i,d,t)*data(i,d,t); 
                                }

                                sum+=0.5*V( std::sqrt(r2) );
                            }
                    }


 
                return sum;
            };


            template<class functor_t>
            Real reduceOneBodyPotentialChin(const functor_t & V,const Eigen::Tensor<Real,3> & data, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange, const std::array<Real,2> & coefficients)
            {
                Real sum=0;
                int tFirstA=((timeRange[0] + 1)/3 )*3;
                int tLastA=(timeRange[1] /3 )*3;
                        
                for(int t=timeRange[0];t<=timeRange[1] ; t++ )
                    {

                        Real alpha=coefficients[t%2];

                        alpha=  ( (t == tFirstA) or (t == tLastA) ) ? alpha/2 : alpha;


                    for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                        {
                            Real r2=0;
                            for (int d=0;d<getDimensions();d++)
                            {
                            r2+=data(i,d,t)*data(i,d,t); 
                            }

                            sum+=V( std::sqrt(r2) )*alpha;
                        }
                    }
                
                    

                   


 
                return sum;
            };


            template<class functor_t>
            Real reduceOneBodyPotentialChin(const functor_t & V,const Eigen::Tensor<Real,3> & data, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange, const std::array<Real,2> & coefficients,const pimcConfigurations::tags_t & mask)
            {

                // beads of type A
                Real sum=0;
                int tFirstA=((timeRange[0] + 1)/3 )*3;
                Real alpha=coefficients[0];

                for(int t=tFirstA;t<=timeRange[1]+1 ; t+=3 )
                    {

                        for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                        {
                            Real r2=0;
                            for (int d=0;d<getDimensions();d++)
                            {
                            r2+=data(i,d,t)*data(i,d,t); 
                            }

                            sum+=V( std::sqrt(r2) )*0.5*(mask(i,t) + mask(i,t-1));
                        }
                    }
                
                // beads of type B

                int tFirstB= ( (tFirstA -2) < timeRange[0] ) ? tFirstA + 1 : tFirstA - 2;
                alpha=coefficients[1];
                for(int t=tFirstB;t<=timeRange[1]+1 ; t+=3 )
                    {

                        for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                        {
                            Real r2=0;
                            for (int d=0;d<getDimensions();d++)
                            {
                            r2+=data(i,d,t)*data(i,d,t); 
                            }

                            sum+=V( std::sqrt(r2) )*mask(i,t) ;
                        }
                    }

                // beads of type C

                alpha=coefficients[0];
                int tFirstC= (tFirstA -1 < timeRange[0] ) ? tFirstA + 2 : tFirstA - 1;

                for(int t=tFirstC;t<=timeRange[1]+1 ; t+=3 )
                    {

                        for (size_t i=particleRange[0];i<=particleRange[1];i++ )
                        {
                            Real r2=0;
                            for (int d=0;d<getDimensions();d++)
                            {
                            r2+=data(i,d,t)*data(i,d,t); 
                            }

                            sum+=V( std::sqrt(r2) )*mask(i,t) ;
                        }
                    }

                    

                   


 
                return sum;
            };






            template<class functor_t>
            Real reduceOneBodyPotentialIsotropic(const functor_t & V,const Eigen::Tensor<Real,3> & data, std::array<int ,2 > timeRange, std::array<int, 2>  particleRange, const pimcConfigurations::tags_t & mask)
            {
                Real sum=0;

                // tail part
                {
                int t=timeRange[0];
                for (int i=particleRange[0];i<=particleRange[1];i++ )
                        {
                            Real r2=0;
                            for (int d=0;d<getDimensions();d++)
                            {
                            r2+=data(i,d,t)*data(i,d,t); 
                            }
                            sum+=
                            0.5*(mask(i,t) )* V(std::sqrt(r2) );
                            ;
                        }
                }

                // central part
                for(int t=timeRange[0]+1;t<=timeRange[1] ; t++ )
                    {
                    for (int i=particleRange[0];i<=particleRange[1];i++ )
                        {
                            Real r2=0;
                            for (int d=0;d<getDimensions();d++)
                            {
                            r2+=data(i,d,t)*data(i,d,t); 
                            }


                            sum+=
                            0.5*(mask(i,t) + mask(i,t-1) )*V(std::sqrt(r2));
                        }
                    }

                // end part
                {
                int t = timeRange[1] + 1;
                for (int i=particleRange[0];i<=particleRange[1];i++ )
                        {
                            Real r2=0;
                            for (int d=0;d<getDimensions();d++)
                            {
                            r2+=data(i,d,t)*data(i,d,t); 
                            }

                            sum+=
                            0.5*(mask(i,t-1) )* V(std::sqrt(r2));
                            ;
                        }
                }
                return sum;
            };

        }
    }

}
