//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STRAWMAN_DRIVER__
#define __STRAWMAN_DRIVER__
#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
namespace PhysBAM{

template<class TV> class STRAWMAN_EXAMPLE;

template<class TV>
class STRAWMAN_DRIVER : public DRIVER<TV>
{
    typedef DRIVER<TV> BASE;
    typedef typename TV::SCALAR T;
    using BASE::time;
    using BASE::current_frame;
    using BASE::output_number;

    STRAWMAN_EXAMPLE<TV>& example;
  public:
    STRAWMAN_DRIVER(STRAWMAN_EXAMPLE<TV>& example);
    virtual ~STRAWMAN_DRIVER();

//#####################################################################
    void Initialize() PHYSBAM_OVERRIDE;
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
//#####################################################################
};

}
#endif
