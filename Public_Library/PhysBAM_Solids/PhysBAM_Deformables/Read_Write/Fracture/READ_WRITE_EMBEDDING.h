//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_EMBEDDING
//#####################################################################
#ifndef __READ_WRITE_EMBEDDING__
#define __READ_WRITE_EMBEDDING__

#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDING.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<EMBEDDING<TV>,RW>:public Read_Write<STRUCTURE<TV>,RW>
{
public:
};
}
#endif
