//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __MPLE_DRIVER__
#define __MPLE_DRIVER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
class MPLE_DRIVER: public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    
    enum WORKAROUND{ghost=3};

public:
    
    ARRAY<TV> w;             // data
    GRID<TV> grid;           // grid

    ARRAY<T,TV_INT>* u;      // segmentation function
    ARRAY<T,TV_INT>* u_new;  // new segmentations funtion
    
    MPLE_DRIVER()
    {
        u=new ARRAY<T,TV_INT>;
        u_new=new ARRAY<T,TV_INT>;
    }

    ~MPLE_DRIVER()
    {
        delete u;
        delete u_new;
    }

    void Initialize()
    {
        u->Resize(grid.Node_Indices(ghost),false);
        u_new->Resize(grid.Node_Indices(ghost),false);
    }
};
}

#endif
