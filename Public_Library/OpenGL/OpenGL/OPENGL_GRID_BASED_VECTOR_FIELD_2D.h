//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_BASED_VECTOR_FIELD_2D
//#####################################################################
#ifndef __OPENGL_GRID_BASED_VECTOR_FIELD_2D__
#define __OPENGL_GRID_BASED_VECTOR_FIELD_2D__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
namespace PhysBAM{

template<class T_input>
class OPENGL_GRID_BASED_VECTOR_FIELD_2D:public OPENGL_VECTOR_FIELD_2D<ARRAY<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    using OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >::size;using OPENGL_OBJECT<T>::World_Space_Box;

    ARRAY<TV> vector_field,vector_locations;
    GRID<TV>& grid;
    ARRAY<TV,TV_INT>& V;
    TV_INT selected_cell;
    TV_INT selected_node;

    OPENGL_GRID_BASED_VECTOR_FIELD_2D(STREAM_TYPE stream_type,GRID<TV>& grid,ARRAY<VECTOR<T,2>,VECTOR<int,2> >& V);
    virtual ~OPENGL_GRID_BASED_VECTOR_FIELD_2D();

    void Update();  // Call when grid/V change

    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    void Print_Selection_Info(std::ostream& stream) const override;

};
}
#endif
