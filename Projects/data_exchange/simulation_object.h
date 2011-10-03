#ifndef __DATA_EXCHANGE_SIMULATION_OBJECT__
#define __DATA_EXCHANGE_SIMULATION_OBJECT__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>

/*
  Base class for objects in the simulation.  Derived classes are
  deformable_object, scripted_geometry, and ground_plane.
 */

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

/*
  Base class for holding the results of a frame of simulation.  The only derived
  class is deformable_body_output.
 */

struct simulation_object_output
{
    virtual ~simulation_object_output() {}
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    }
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(simulation_object_output)
}

#endif
