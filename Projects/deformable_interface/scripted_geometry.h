#ifndef __DATA_EXCHANGE_SCRIPTED_GEOMETRY__
#define __DATA_EXCHANGE_SCRIPTED_GEOMETRY__

#include "simulation_object.h"
#include "polygon_mesh.h"

namespace data_exchange{
struct scripted_geometry : public simulation_object
{
    // Reference geometry
    polygon_mesh mesh;
    std::vector<vf3> position;
};
}

#endif
