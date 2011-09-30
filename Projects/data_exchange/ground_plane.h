#ifndef __DATA_EXCHANGE_GROUND_PLANE__
#define __DATA_EXCHANGE_GROUND_PLANE__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "simulation_object.h"
#include "polygon_mesh.h"

namespace data_exchange{
struct ground_plane : public simulation_object
{
    vf3 position;
    vf3 normal;

    ground_plane(): position(0,0,0), normal(0,1,0) {}

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<simulation_object>(*this);
        ar & position & normal;
    }
};
}

#endif
