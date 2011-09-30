#ifndef __DATA_EXCHANGE_DEFORMABLE_BODY__
#define __DATA_EXCHANGE_DEFORMABLE_BODY__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "geometry.h"

namespace data_exchange{

struct simulation_object
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & position & velocity;
    }
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(simulation_object)

struct deformable_body : public simulation_object
{
    polygon_mesh mesh;

    vector<vf3> position;
    vector<vf3> velocity;
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

struct ground_plane : public simulation_object
{
    vf3 position;
    vf3 normal;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<simulation_object>(*this);
        ar & position & normal;
    }
};

struct scripted_geometry : public simulation_object
{
    // Reference geometry
    polygon_mesh mesh;
    vector<vf3> position;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<simulation_object>(*this);
        ar & mesh & position;
    }
};

struct deformable_body_simulation
{
    vector<simulation_object*> simulation_objects;

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

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar.register_type(static_cast<deformable_body*>(NULL));
        ar.register_type(static_cast<ground_plane*>(NULL));
        ar.register_type(static_cast<scripted_geometry*>(NULL));
        ar & simulation_objects;
    }
};

struct deformable_body_output_frame
{
    vector<vf3> position;
    vector<vf3> velocity;

    deformable_body_output_frame() {}

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & position & velocity;
    }
};

}

#endif
