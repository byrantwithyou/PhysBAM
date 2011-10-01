#ifndef __DATA_EXCHANGE_DEFORMABLE_BODY_SIMULATION__
#define __DATA_EXCHANGE_DEFORMABLE_BODY_SIMULATION__

#include "force.h"
#include "simulation_object.h"
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include "deformable_body.h"
#include "gravity_force.h"
#include "ground_plane.h"
#include "scripted_geometry.h"
#include "volumetric_force.h"

namespace data_exchange{
struct deformable_body_simulation
{
    std::vector<simulation_object*> simulation_objects;
    std::vector<force*> simulation_forces;

    float dt;
    int number_frames;
    std::string output_directory;
    std::string filename_base;
    // If output_directory = "output" and filename_base = "frame_data"
    // and number_frames = 100, then output will be found in:
    // output/frame_data.0  (Initial configuration)
    // output/frame_data.1
    // output/frame_data.2
    // ...
    // output/frame_data.100

    deformable_body_simulation():
        dt(1./24), number_frames(120),
       output_directory("output"), filename_base("frame_data")
    {}

    virtual ~deformable_body_simulation()
    {
        for(size_t i=0; i<simulation_objects.size(); i++) delete simulation_objects[i];
        for(size_t i=0; i<simulation_forces.size(); i++) delete simulation_forces[i];
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar.register_type(static_cast<deformable_body*>(NULL));
        ar.register_type(static_cast<ground_plane*>(NULL));
        ar.register_type(static_cast<scripted_geometry*>(NULL));
        ar & simulation_objects;
        ar.register_type(static_cast<volumetric_force*>(NULL));
        ar.register_type(static_cast<gravity_force*>(NULL));
        ar & simulation_forces;
    }
};
}

#endif
