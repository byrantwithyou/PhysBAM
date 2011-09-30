#ifndef __DATA_EXCHANGE_FORCES__
#define __DATA_EXCHANGE_FORCES__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "geometry.h"

namespace data_exchange{

struct force
{
    std::vector<int> affected_bodies;

    virtual ~force(){}
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & affected_bodies;
    }
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(force)

struct volumetric_force : public force
{
    float stiffness;
    float poissons_ratio;
    float damping;

    volumetric_force(): stiffness(1e3f), poissons_ratio(.45f), damping(1.f) {}
    virtual ~volumetric_force() {}

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<force>(*this);
        ar & stiffness & poissons_ratio & damping;
    }
};

struct gravity_force : public force
{
    float magnitude;
    vf3 direction;

    gravity_force(): magnitude(9.8f), direction(0,-1,0) {}

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<force>(*this);
        ar & magnitude & direction;
    }
};
}

#endif
