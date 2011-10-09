#ifndef __DATA_EXCHANGE_GRAVITY_FORCE__
#define __DATA_EXCHANGE_GRAVITY_FORCE__

#include "force.h"

namespace data_exchange{
struct gravity_force : public force
{
    float magnitude;
    vf3 direction;

    gravity_force(): magnitude(9.8f), direction(0,-1,0) {}
};
}

#endif
