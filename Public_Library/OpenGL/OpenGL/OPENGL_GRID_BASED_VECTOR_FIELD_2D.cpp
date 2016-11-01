//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/RANGE.h>
#include <OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <OpenGL/OpenGL/OPENGL_GRID_BASED_VECTOR_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// OPENGL_GRID_BASED_VECTOR_FIELD_2D
//#####################################################################
template<class T> OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
OPENGL_GRID_BASED_VECTOR_FIELD_2D(STREAM_TYPE stream_type,GRID<TV>& grid,ARRAY<TV,TV_INT>& V)
    :OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >(stream_type,vector_field,vector_locations),grid(grid),V(V),
    selected_cell(-1,-1),selected_node(-1,-1)
{}
//#####################################################################
// ~OPENGL_GRID_BASED_VECTOR_FIELD_2D
//#####################################################################
template<class T> OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
~OPENGL_GRID_BASED_VECTOR_FIELD_2D()
{}
//#####################################################################
// Update
//#####################################################################
template<class T> void OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
Update()
{
    int idx=0;
    vector_field.Resize(V.Size().Product());
    vector_locations.Resize(V.Size().Product());
    for(RANGE_ITERATOR<TV::m> it(V.domain);it.Valid();it.Next()){
        vector_field(idx)=V(it.index);
        vector_locations(idx)=grid.X(it.index);
        idx++;}
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
Bounding_Box() const
{
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
Print_Selection_Info(std::ostream& stream) const
{
    // TODO: interpolate to particles
    if(selected_node.x>=0 && !grid.Is_MAC_Grid()){
        if(V.Valid_Index(selected_node)) stream<<V(selected_node);}
    if(selected_cell.x>=0 && grid.Is_MAC_Grid()){
        if(V.Valid_Index(selected_cell)) stream<<V(selected_cell);}
    stream<<std::endl;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_GRID_BASED_VECTOR_FIELD_2D<float>;
template class OPENGL_GRID_BASED_VECTOR_FIELD_2D<double>;
}
