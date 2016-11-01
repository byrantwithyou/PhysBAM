//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_BASED_VECTOR_FIELD_3D
//#####################################################################
#ifndef __OPENGL_GRID_BASED_VECTOR_FIELD_3D__
#define __OPENGL_GRID_BASED_VECTOR_FIELD_3D__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
namespace PhysBAM{

template<class TV> class RANGE;

template<class T_input>
class OPENGL_GRID_BASED_VECTOR_FIELD_3D:public OPENGL_VECTOR_FIELD_3D<T_input>
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    using OPENGL_VECTOR_FIELD_3D<T>::slice;
    using OPENGL_VECTOR_FIELD_3D<T>::vector_field;using OPENGL_VECTOR_FIELD_3D<T>::vector_locations;
    using OPENGL_VECTOR_FIELD_3D<T>::size;
    using OPENGL_VECTOR_FIELD_3D<T>::World_Space_Box;

    GRID<TV>& grid;
    ARRAY<TV,TV_INT>& V;
    int max_vectors_3d;
    TV_INT selected_cell;
    TV_INT selected_node;
    
    OPENGL_GRID_BASED_VECTOR_FIELD_3D(STREAM_TYPE stream_type,GRID<TV>& grid,ARRAY<TV,TV_INT>& V);
    virtual ~OPENGL_GRID_BASED_VECTOR_FIELD_3D();

    void Update();  // Call when grid/V change

    virtual RANGE<TV> Bounding_Box() const override;
    virtual void Slice_Has_Changed() override { Update(); }
    void Print_Selection_Info(std::ostream& stream) const override;

};
}
#endif
