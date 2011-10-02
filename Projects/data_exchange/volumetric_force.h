#ifndef __DATA_EXCHANGE_VOLUMETRIC_FORCE__
#define __DATA_EXCHANGE_VOLUMETRIC_FORCE__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "force.h"

namespace data_exchange{
struct volumetric_force : public force
{
    float stiffness;
    float poissons_ratio;
    float damping;
    std::vector<int> bodies_affected;

    volumetric_force(): stiffness(1e3f), poissons_ratio(.45f), damping(1.f) {}
    virtual ~volumetric_force() {}

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<force>(*this);
        ar & stiffness & poissons_ratio & damping & bodies_affected;
    }
};
}

#endif
