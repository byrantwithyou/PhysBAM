#ifndef __DATA_EXCHANGE_VOLUMETRIC_FORCE__
#define __DATA_EXCHANGE_VOLUMETRIC_FORCE__

#include "force.h"

/*
  Members:

  * stiffness

  Indicates how resistant the object is to changing its shape.  Any positive
  number is valid, but values will normally be rather large, such as the default
  of 1000.

  * poissons_ratio

  This indicates the degree to which an object will resist changes in volume.
  Reasonable values are perhaps 0.2 to 0.45.

  * damping

  Simulators typically introduce some damping into a simulation automatically.
  This parameter specifies that additional damping should be included.  It is
  generally good to include a little extra.  Any positive or zero number is
  fine.

  * bodies_affected

  This is a list of body indices to whom this force should be applied.  The
  indexing follows the indexing of the simulated_

 */

namespace data_exchange{
struct volumetric_force : public force
{
    float stiffness;
    float poissons_ratio;
    float damping;

    volumetric_force(): stiffness(1e3f), poissons_ratio(.45f), damping(1.f) {}
    virtual ~volumetric_force() {}
};
}

#endif
