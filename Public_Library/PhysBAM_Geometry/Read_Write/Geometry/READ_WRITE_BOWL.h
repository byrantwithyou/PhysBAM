//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BOWL
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_BOWL__
#define __READ_WRITE_BOWL__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOWL.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<BOWL<T>,RW>
{
public:
    static void Read(std::istream& input,BOWL<T>& object)
    {Read_Binary<RW>(input,object.hole_radius,object.depth,object.thickness,object.height,object.inner_radius,object.outer_radius);}

    static void Write(std::ostream& output,const BOWL<T>& object)
    {Write_Binary<RW>(output,object.hole_radius,object.depth,object.thickness,object.height,object.inner_radius,object.outer_radius);}
};
}
#endif
#endif
