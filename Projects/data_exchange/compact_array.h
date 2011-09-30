#ifndef __DATA_EXCHANGE_PRIMITIVES__
#define __DATA_EXCHANGE_PRIMITIVES__

#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace data_exchange{
template<class T>
struct compact_array
{
    size_t size;
    T constant_value;
    vector<T> all_values;

    fixed_vector(): size(0), constant_value(0) {}

    void set_constant_value(size_t array_size, const T& value)
    {
        size = array_size;
        constant_value = value;
    }

    void populate_full_array()
    {
        if(!all_values.size())
        {
            all_values.resize(size);
            std::fill(all_values.begin(), all_values.end(), constant_value);
        }
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        for(int i=0; i<n; i++) ar & data[i];
    }
};

typedef fixed_vector<float,3> vf3;
typedef fixed_vector<int,3> vi3;
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
} // namespace serialization
} // namespace boost

#endif
