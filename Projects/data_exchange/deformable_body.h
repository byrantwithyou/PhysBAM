#ifndef __DATA_EXCHANGE_DEFORMABLE_BODY__
#define __DATA_EXCHANGE_DEFORMABLE_BODY__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "simulation_object.h"
#include "polygon_mesh.h"

namespace data_exchange{
struct deformable_body : public simulation_object
{
    polygon_mesh mesh;

    std::vector<vf3> position;
    std::vector<vf3> velocity;
    float mass;

    deformable_body(): mass(1) {}

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<simulation_object>(*this);
        ar & mesh & position & velocity & mass;
    }
};

struct deformable_body_output : public simulation_object_output
{
    std::vector<vf3> position;
    std::vector<vf3> velocity;

    deformable_body_output() {}

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<simulation_object_output>(*this);
        ar & position & velocity;
    }
};
}

#endif
