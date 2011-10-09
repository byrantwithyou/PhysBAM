#ifndef __DATA_EXCHANGE_FIXED_MATRIX__
#define __DATA_EXCHANGE_FIXED_MATRIX__

#include <vector>

namespace data_exchange{
template<class T, int m,int n>
struct fixed_matrix
{
    T data[m][n];
};

typedef fixed_matrix<float,1,1> mf1;
typedef fixed_matrix<float,2,2> mf2;
typedef fixed_matrix<float,3,3> mf3;
typedef fixed_matrix<float,4,4> mf4;
}

#endif
