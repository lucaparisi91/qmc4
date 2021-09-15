#include <vector>
#include<string>
#include <iostream>
#include <type_traits>


struct arrayData
{
    double * old_positions;
    double * new_positions;
    double * dx;

};


namespace operators
{

struct operatorBase {};

struct laplacian : public operatorBase
{
    public:

    inline double operator()( int i, const double * positions , const double * dx ) const
    {
        return (  -2 * positions[i] +   positions[i+1] + positions[i-1])/dx[0]*dx[0];
    };

    private:


};

struct cubicTerm : public operatorBase
{
    inline double operator()( int i, const double * positions , const double * dx ) const
    {
        return  positions[i]*positions[i]*positions[i];
    };

};


class reductionOperator : public operatorBase{};


struct norm : public reductionOperator
{

    inline double operator()(int i,const double * data ) const
    {
    
        return data[i]*data[i];
    };

};


template<class T1,class T2>
class sumOperators : public operatorBase
{
    public:
    sumOperators(T1 left,T2 right) : _left(left),_right(right) {}


    template<class ... Args>
    auto operator()(Args&& ... args) {

        
        return _left(std::forward<Args>(args)...)   + _right(std::forward<Args>(args)...)  ; 
        }
    
    private:

    T1 _left;
    T2 _right;
}; 





template<class T1,class T2>
class productOperators : public operatorBase
{
    public:
    productOperators(T1 left,T2 right) : _left(left),_right(right) {}

    
    template<class ... Args>
    auto operator()(Args&& ... args) {

        
        return _left(std::forward<Args>(args)...)   * _right(std::forward<Args>(args)...)  ; 
        }
    
    private:

    T1 _left;
    T2 _right;
}; 

template<class T1,class T2>
class productScalarOperator : public operatorBase
{
    public:
    productScalarOperator(T1 left,T2 right) : _left(left),_right(right) {}

    
    template<class ... Args>
    auto operator()(Args&& ... args) {
        return _left   * _right(std::forward<Args>(args)...)  ; 
        }
    
    private:

    T1 _left;
    T2 _right;
}; 



template<class T1,class T2, std::enable_if_t< std::is_base_of<operatorBase,T1>::value  >* =nullptr   ,std::enable_if_t< std::is_base_of<operatorBase,T2>::value  >* =nullptr  >
sumOperators<T1,T2> operator+(T1 left, T2 right)
{
    return sumOperators<T1,T2>(left,right);
}

template<class T1,class T2, std::enable_if_t< !std::is_base_of<operatorBase,T1>::value  >* =nullptr   ,std::enable_if_t< std::is_base_of<operatorBase,T2>::value  >* =nullptr  >
auto operator*(T1 left, T2 right)
{
    return productScalarOperator<T1,T2>(left,right);
}

}


class loop
{

public:
    loop(){}

    template<class op_t,class ... Args>
    auto apply(int iStart, int iEnd, op_t &&  currentOp, Args&&...  args)
    {
        if constexpr (std::is_base_of<operators::reductionOperator,op_t>() )
        {
            return reduceSum(iStart,iEnd,currentOp,std::forward<Args>(args)  ...);
        }
        else
        {
            for ( int i=iStart;i<=iEnd;i++)
            {
                currentOp(i,std::forward<Args>(args) ... );
            }

        
        }
      
    };


    template<class op_t,class ... Args>
    auto reduceSum(int iStart, int iEnd, op_t &&  currentOp, Args&&... args)
    {
        using return_t=decltype(currentOp(0,std::forward<Args>(args)...) );
        return_t sum{0};

        for (int i=iStart;i<=iEnd;i++)
        {
            sum+=currentOp(i,std::forward<Args>(args)...);
        }
        return sum;
    }

};


int main(int argc,char** argv)
{

    const int N=0;

    arrayData data;
    data.old_positions = (double *)aligned_alloc( 128, N*sizeof(double) );
    data.new_positions = (double *)aligned_alloc( 128, N*sizeof(double) );
    data.dx = (double *)aligned_alloc( 128, 1*sizeof(double) );
    data.dx[0]=0.1;


    loop loop1d;


    auto f =   - 0.5 * operators::laplacian() + 3* operators::cubicTerm();

    

    loop1d.apply(1,N-2 , data.new_positions , data.old_positions );



    //auto sum= loop1d.reduce(data,0,N-1,norm());


}