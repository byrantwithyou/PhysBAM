#include <Tools/Arrays/ARRAY.h>
using namespace PhysBAM;

//ARRAY, the array can be used as a cycle.

template<class T>
class CYCLIC_ARRAY:public ARRAY <T>
{
public:
    T& operator()(const int i)
    {
        int ii= (i-1+m)%m+1;
        assert(ii<=m); 
        return array(ii);}

    const T& operator()(const int i) const
    {
        int ii= (i-1+m)%m+1;
        assert(ii<=m); 
        return array(ii);
    }
};
