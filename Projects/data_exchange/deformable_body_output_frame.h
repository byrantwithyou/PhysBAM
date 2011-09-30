#ifndef __DATA_EXCHANGE_DEFORMABLE_BODY_OUTPUT_FRAME__
#define __DATA_EXCHANGE_DEFORMABLE_BODY_OUTPUT_FRAME__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "primitivies.h"

namespace data_exchange{
struct deformable_body_output_frame
{
    std::vector<vf3> position;
    std::vector<vf3> velocity;

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
