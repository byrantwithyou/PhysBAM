//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_CALLBACKS
//##################################################################### 
#ifndef __CONSERVATION_CALLBACKS__
#define __CONSERVATION_CALLBACKS__

#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/GRID.h>
namespace PhysBAM{

template<class T>
class CONSERVATION_CALLBACKS
{
public:
    CONSERVATION_CALLBACKS()
    {}

    virtual ~CONSERVATION_CALLBACKS()
    {}

    virtual void Get_Neumann_Face_Location(const GRID<VECTOR<T,1> >& grid_1d,const int face_index,T& location) const
    {location=grid_1d.Face(FACE_INDEX<1>(1,VECTOR<int,1>(face_index))).x;}
};
}
#endif
