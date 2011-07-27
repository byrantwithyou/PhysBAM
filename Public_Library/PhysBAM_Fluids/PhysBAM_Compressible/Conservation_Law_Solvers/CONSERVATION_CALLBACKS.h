//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_CALLBACKS
//##################################################################### 
#ifndef __CONSERVATION_CALLBACKS__
#define __CONSERVATION_CALLBACKS__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
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
    {location=grid_1d.Face(1,VECTOR<int,1>(face_index)).x;}
};
}
#endif
