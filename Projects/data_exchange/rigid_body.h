#ifndef __DATA_EXCHANGE_RIGID_BODY__
#define __DATA_EXCHANGE_RIGID_BODY__

#include "geometry.h"

namespace data_exchange{

struct rigid_body: public simulation_object
{
    polygon_mesh mesh;

    vf3 pos;
    mf3 orientation;
    vf3 vel;
    vf3 angular_vel;
    float mass;
    mf3 inertia;
    float density;
    bool scripted;

    vector<vf3> vel; // velocity
    vector<float> mass; // mass

    rigid_body(): constant_mass(0) {}

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & mesh & vel & mass;
    }
};

}

#endif
