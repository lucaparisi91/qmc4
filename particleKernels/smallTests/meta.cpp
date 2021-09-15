#include <vector>
#include<string>
#include <iostream>


template <class... Ts> struct cell {

    cell(){nParticles=0;}
    void addParticle(){};

    protected:

    int nParticles;
};


template <class T, class... Ts>
struct cell<T, Ts...> :  cell<Ts...> 
{
    using cell<Ts... >::nParticles;


  cell()   {}


    template<int d,std::enable_if_t<d==0>* = nullptr   >
        const auto &  get(){return field;}
    template<int d,std::enable_if_t<d!=0>* = nullptr   >
        const auto &  get(){return cell<Ts...>::template get<d-1>()       ;}
    
    void addParticle( T t, Ts... args)
    {
        field.push_back(t);
        return cell<Ts...>::addParticle(args...);
        nParticles++;
    }

    private:
    
    std::vector<T> field ;
};








int main(int argc,char** argv)
{

    struct field
    {
        int a;
        int b;
        int c;
    };

    cell<int,int> indexCell;


    const auto & indices = indexCell.get<0>();
    const auto & leavingIndices = indexCell.get<1>();

    indexCell.addParticle(1,1);
    
}