#ifndef __DATA_EXCHANGE_SIMULATION_OBJECT__
#define __DATA_EXCHANGE_SIMULATION_OBJECT__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>

namespace data_exchange{

struct simulation_object
{
    virtual ~simulation_object() {}
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    }
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(simulation_object)
}

#endif
