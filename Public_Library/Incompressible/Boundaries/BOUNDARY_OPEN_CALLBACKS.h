//#####################################################################
// Copyright 2006, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_OPEN_CALLBACKS
//##################################################################### 
#ifndef __BOUNDARY_OPEN_CALLBACKS__
#define __BOUNDARY_OPEN_CALLBACKS__

#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV>
class BOUNDARY_OPEN_CALLBACKS
{
    typedef typename TV::SCALAR T;
public:

    BOUNDARY_OPEN_CALLBACKS()
    {}

    virtual ~BOUNDARY_OPEN_CALLBACKS()
    {}

//#####################################################################
    virtual T Open_Water_Phi(const TV& location,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;}
//#####################################################################
};
}
#endif
