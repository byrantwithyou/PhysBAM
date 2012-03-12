//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_STRUCTURE_LIST
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_STRUCTURE_LIST__
#define __READ_WRITE_STRUCTURE_LIST__

#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <string>
namespace PhysBAM{

template<class RW,class TV,class ID>
class Read_Write<STRUCTURE_LIST<TV,ID>,RW>
{
    typedef DYNAMIC_LIST<STRUCTURE<TV>,ID> OBJECT_BASE;
public:
};
}

#endif
#endif
