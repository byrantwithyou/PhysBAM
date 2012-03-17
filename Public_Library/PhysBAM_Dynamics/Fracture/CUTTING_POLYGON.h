//#####################################################################
// Copyright 2006, Kevin Der.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_POLYGON
//##################################################################### 
#ifndef __CUTTING_POLYGON__
#define __CUTTING_POLYGON__

#include <PhysBAM_Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
namespace PhysBAM{

class CUTTING_POLYGON
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    enum POLYGON_TYPE {FACE_BOUNDARY,FACE_INTERIOR,TRIANGLE_CLIPPED};
    int polygon_index; // indexes into polygon_mesh
    int simplex_owner; // indexes into cutting_simplices
    bool flipped;
    POLYGON_TYPE polygon_type;

    CUTTING_POLYGON()
    {}

    CUTTING_POLYGON(const int polygon_index,const int simplex_owner,const bool flipped,const POLYGON_TYPE polygon_type)
        :polygon_index(polygon_index),simplex_owner(simplex_owner),flipped(flipped),polygon_type(polygon_type)
    {}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,polygon_index,simplex_owner,flipped,polygon_type);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,polygon_index,simplex_owner,flipped,polygon_type);}
//#####################################################################    
};
}
#endif
