//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAGGED_FORCE
//#####################################################################
#ifndef __LAGGED_FORCE__
#define __LAGGED_FORCE__

#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class LAGGED_FORCE:public DEFORMABLES_FORCES<TV>
{
private:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:

    LAGGED_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input)
        :DEFORMABLES_FORCES<TV>(particles_input)
    {}

    virtual ~LAGGED_FORCE()
    {}

    virtual void Lagged_Update_Position_Based_State(const T time)=0;
//#####################################################################
};
}
#endif
