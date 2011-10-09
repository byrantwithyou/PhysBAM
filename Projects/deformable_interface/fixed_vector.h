#ifndef __DATA_EXCHANGE_FIXED_VECTOR__
#define __DATA_EXCHANGE_FIXED_VECTOR__

#include <boost/static_assert.hpp>

namespace data_exchange{
template<class T, int n>
struct fixed_vector
{
    T data[n];

    fixed_vector(){}

    fixed_vector(const T& a) {BOOST_STATIC_ASSERT(n==1); data[0]=a;}

    fixed_vector(const T& a, const T& b) {BOOST_STATIC_ASSERT(n==2); data[0]=a; data[1]=b;}

    fixed_vector(const T& a, const T& b, const T& c) {BOOST_STATIC_ASSERT(n==3); data[0]=a; data[1]=b; data[2]=c;}

    fixed_vector(const T& a, const T& b, const T& c, const T& d) {BOOST_STATIC_ASSERT(n==4); data[0]=a; data[1]=b; data[2]=c; data[3]=d;}
};

typedef fixed_vector<float,1> vf1;
typedef fixed_vector<float,2> vf2;
typedef fixed_vector<float,3> vf3;
typedef fixed_vector<float,4> vf4;
typedef fixed_vector<int,1> vi1;
typedef fixed_vector<int,2> vi2;
typedef fixed_vector<int,3> vi3;
typedef fixed_vector<int,4> vi4;
}

#endif
