#ifndef __DATA_EXCHANGE_SIMULATION_OBJECT__
#define __DATA_EXCHANGE_SIMULATION_OBJECT__

/*
  Base class for objects in the simulation.  Derived classes are
  deformable_object, scripted_geometry, and ground_plane.
 */

namespace data_exchange{
struct simulation_object
{
    virtual ~simulation_object() {}
};
}

#endif
