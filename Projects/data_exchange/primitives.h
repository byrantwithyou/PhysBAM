#ifndef __DATA_EXCHANGE_PRIMITIVES__
#define __DATA_EXCHANGE_PRIMITIVES__

#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

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

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        for(int i=0; i<n; i++) ar & data[i];
    }
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

namespace boost {
namespace serialization {
template<class T, int n>
void serialize(boost::archive::binary_iarchive & ar, std::vector<data_exchange::fixed_vector<T,n> > & g, const unsigned int version)
{
    size_t s;
    ar >> s;
    g.resize(s);
    ar.load_binary((unsigned char*)&g[0].data, g.size()*n*sizeof(T));
}
template<class T, int n>
void serialize(boost::archive::binary_oarchive & ar, std::vector<data_exchange::fixed_vector<T,n> > & g, const unsigned int version)
{
    size_t s = g.size();
    ar << s;
    ar.save_binary((unsigned char*)&g[0].data, g.size()*n*sizeof(T));
}
}
}

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
