//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BONE__
#define __BONE__
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>

namespace PhysBAM{

template<class T>
class BONE
{
    typedef VECTOR<T,3> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    FRAME<TV> transform,targeted_transform;
    FRAME<TV> rotation,targeted_rotation;
    FRAME<TV> translation,targeted_translation;
    T length;

    BONE()
        :length(0)
    {}

    ~BONE()
    {}

    bool operator==(const BONE<T>& bone_input) const
    {return (length==bone_input.length&&transform==bone_input.transform&&targeted_transform==bone_input.targeted_transform);}

    bool operator!=(const BONE<T>& bone_input) const
    {return !(*this==bone_input);}

    BONE<T>& operator=(const BONE<T>& bone_input)
    {transform=bone_input.transform;targeted_transform=bone_input.targeted_transform;length=bone_input.length;return *this;}

    void Rescale(T scaling_factor)
    {transform.t*=scaling_factor;targeted_transform.t*=scaling_factor;length*=scaling_factor;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,transform,targeted_transform,length);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,transform,targeted_transform,length);}
//#####################################################################    
};
template<class T> inline std::istream& operator>>(std::istream& input,BONE<T>& bone)
{PHYSBAM_NOT_IMPLEMENTED();}

template<class T> inline std::ostream& operator<<(std::ostream& output,const BONE<T>& bone)
{output<<bone.targeted_transform<<" "<<bone.transform<<" "<<bone.length;return output;}
}
#endif
