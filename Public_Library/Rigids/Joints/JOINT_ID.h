//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class JOINT_ID
//#####################################################################
#ifndef __JOINT_ID__
#define __JOINT_ID__

#include <Core/Data_Structures/ELEMENT_ID.h>
namespace PhysBAM{

#ifndef COMPILE_ID_TYPES_AS_INT
PHYSBAM_DECLARE_ELEMENT_ID(JOINT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
#else
typedef int JOINT_ID;
#endif
}
#endif
