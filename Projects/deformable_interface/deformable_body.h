#ifndef __DATA_EXCHANGE_DEFORMABLE_BODY__
#define __DATA_EXCHANGE_DEFORMABLE_BODY__

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
};
}

#endif
