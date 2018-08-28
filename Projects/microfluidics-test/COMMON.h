//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMMON__
#define __COMMON__
#include <Core/Data_Structures/ELEMENT_ID.h>

namespace PhysBAM{
PHYSBAM_DECLARE_ELEMENT_ID(PARTICLE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(EDGE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical);
PHYSBAM_DECLARE_ELEMENT_ID(TRIANGLE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(PIPE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical);
PHYSBAM_DECLARE_ELEMENT_ID(BLOCK_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(VERTEX_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical);
PHYSBAM_DECLARE_ELEMENT_ID(BC_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical);
PHYSBAM_DECLARE_ELEMENT_ID(DOF_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(CODE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical);
PHYSBAM_DECLARE_ELEMENT_ID(LOCAL_V_CODE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical);
PHYSBAM_DECLARE_ELEMENT_ID(LOCAL_P_CODE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical);
}
#endif
