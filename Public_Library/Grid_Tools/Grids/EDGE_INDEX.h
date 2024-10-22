//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EDGE_INDEX
//#####################################################################
#ifndef __EDGE_INDEX__
#define __EDGE_INDEX__

#include <Core/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<int d>
class EDGE_INDEX
{
    typedef VECTOR<int,d> TV_INT;
public:
    EDGE_INDEX()
    {}

    EDGE_INDEX(int axis_input,const TV_INT& index_input)
        :axis(axis_input),index(index_input)
    {}

    int axis;
    TV_INT index;

    bool operator==(const EDGE_INDEX& fi) const
    {return axis==fi.axis && index==fi.index;}

    bool operator!=(const EDGE_INDEX& fi) const
    {return !(*this==fi);}

    TV_INT First_Node_Index() const
    {return index;}

    TV_INT Second_Node_Index() const
    {TV_INT i(index);i(axis)++;return i;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,axis,index);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,axis,index);}
};
template<int d> struct HASH_REDUCE<EDGE_INDEX<d> >
{static int H(const EDGE_INDEX<d>& key){return int_hash(key.axis,HASH_REDUCE<VECTOR<int,d> >::H(key.index));}};
}
#endif
