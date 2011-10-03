#ifndef __DATA_EXCHANGE_FORCE__
#define __DATA_EXCHANGE_FORCE__

#include <vector>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>

namespace data_exchange{
struct force
{
    virtual ~force(){}
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    }
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(force)
}

#endif
