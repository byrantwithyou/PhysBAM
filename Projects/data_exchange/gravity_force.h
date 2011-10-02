#ifndef __DATA_EXCHANGE_GRAVITY_FORCE__
#define __DATA_EXCHANGE_GRAVITY_FORCE__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "force.h"

namespace data_exchange{
struct gravity_force : public force
{
    float magnitude;
    vf3 direction;
    std::vector<int> bodies_affected;

    gravity_force(): magnitude(9.8f), direction(0,-1,0) {}

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<force>(*this);
        ar & magnitude & direction & bodies_affected;
    }
};
}

#endif
