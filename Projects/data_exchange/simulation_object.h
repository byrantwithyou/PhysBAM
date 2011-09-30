#ifndef __DATA_EXCHANGE_SIMULATION_OBJECT__
#define __DATA_EXCHANGE_SIMULATION_OBJECT__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "geometry.h"

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

struct deformable_body : public simulation_object
{
    polygon_mesh mesh;

    std::vector<vf3> position;
    std::vector<vf3> velocity;
    float mass;

    deformable_body(): mass(1) {}
    virtual ~deformable_body() {}

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

    virtual ~ground_plane() {}

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
    std::vector<vf3> position;

    virtual ~scripted_geometry() {}

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
