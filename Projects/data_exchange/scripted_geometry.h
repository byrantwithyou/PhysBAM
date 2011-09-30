#ifndef __DATA_EXCHANGE_SCRIPTED_GEOMETRY__
#define __DATA_EXCHANGE_SCRIPTED_GEOMETRY__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "simulation_object.h"
#include "polygon_mesh.h"

namespace data_exchange{
struct scripted_geometry : public simulation_object
{
    // Reference geometry
    polygon_mesh mesh;
    std::vector<vf3> position;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<simulation_object>(*this);
        ar & mesh & position;
    }
};
}

#endif
