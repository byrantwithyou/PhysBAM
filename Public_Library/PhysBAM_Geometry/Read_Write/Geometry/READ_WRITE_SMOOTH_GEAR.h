//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SMOOTH_GEAR
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SMOOTH_GEAR__
#define __READ_WRITE_SMOOTH_GEAR__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SMOOTH_GEAR.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SMOOTH_GEAR<VECTOR<T,2> >,RW>
{
public:
    static void Read(std::istream& input,SMOOTH_GEAR<VECTOR<T,2> >& object)
    {Read_Binary<RW>(input,object.r,object.s,object.n);object.Compute_Centers();}

    static void Write(std::ostream& output,const SMOOTH_GEAR<VECTOR<T,2> >& object)
    {Write_Binary<RW>(output,object.r,object.s,object.n);}
};

template<class RW,class T>
class Read_Write<SMOOTH_GEAR<VECTOR<T,3> >,RW>
{
public:
    static void Read(std::istream& input,SMOOTH_GEAR<VECTOR<T,3> >& object)
    {Read_Binary<RW>(input,object.w,object.g);}

    static void Write(std::ostream& output,const SMOOTH_GEAR<VECTOR<T,3> >& object)
    {Write_Binary<RW>(output,object.w,object.g);}
};
}
#endif
#endif
