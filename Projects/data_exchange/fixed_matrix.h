#ifndef __DATA_EXCHANGE_FIXED_MATRIX__
#define __DATA_EXCHANGE_FIXED_MATRIX__

#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace data_exchange{
template<class T, int m,int n>
struct fixed_matrix
{
    T data[m][n];

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        for(int i=0; i<m; i++) for(int j=0; j<n; j++) ar & data[i][j];
    }
};

typedef fixed_matrix<float,1,1> mf1;
typedef fixed_matrix<float,2,2> mf2;
typedef fixed_matrix<float,3,3> mf3;
typedef fixed_matrix<float,4,4> mf4;
}

#endif
