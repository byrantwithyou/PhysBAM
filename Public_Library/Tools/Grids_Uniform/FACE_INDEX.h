//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_INDEX
//#####################################################################
#ifndef __FACE_INDEX__
#define __FACE_INDEX__

#include <Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<int d>
class FACE_INDEX
{
    typedef VECTOR<int,d> TV_INT;
public:
    typedef int HAS_UNTYPED_READ_WRITE;

    FACE_INDEX()
        :axis(0)
    {}

    FACE_INDEX(int axis_input,const TV_INT& index_input)
        :axis(axis_input),index(index_input)
    {}

    int axis;
    TV_INT index;

    bool operator==(const FACE_INDEX& fi) const
    {return axis==fi.axis && index==fi.index;}

    TV_INT First_Cell_Index() const
    {TV_INT i(index);i(axis)--;return i;}

    TV_INT Second_Cell_Index() const
    {return index;}

    TV_INT Cell_Index(int j) const
    {TV_INT i(index);i(axis)+=j-2;return i;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,axis,index);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,axis,index);}
};
template<int d>
inline std::ostream& operator<<(std::ostream& output,const FACE_INDEX<d>& fi)
{output<<"("<<fi.axis<<" "<<fi.index<<")";return output;}
template<int d> struct HASH_REDUCE<FACE_INDEX<d> >
{static int H(const FACE_INDEX<d>& key){return int_hash(key.axis,HASH_REDUCE<VECTOR<int,d> >::H(key.index));}};
}
#endif
