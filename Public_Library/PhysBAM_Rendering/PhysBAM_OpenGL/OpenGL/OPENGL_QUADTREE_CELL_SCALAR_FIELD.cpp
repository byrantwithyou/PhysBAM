#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2009, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_QUADTREE_CELL_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_QUADTREE_GRID.h>
using namespace PhysBAM;

//#####################################################################
// Display
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_CELL_SCALAR_FIELD<T,T2>::
Display(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(point_size);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
    for(DYADIC_GRID_ITERATOR_CELL<QUADTREE_GRID<T> > iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next())
        if(value(iterator.Cell_Index())){
            color_map->Lookup(value(iterator.Cell_Index())).Send_To_GL_Pipeline();
            OpenGL_Vertex(iterator.Location(),vertices);}
    OpenGL_Draw_Arrays(GL_POINTS,2,vertices);
    glPopAttrib();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_QUADTREE_CELL_SCALAR_FIELD<T,T2>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.uniform_grid.domain.min_corner.x,grid.uniform_grid.domain.max_corner.x,grid.uniform_grid.domain.min_corner.y,grid.uniform_grid.domain.max_corner.y,0,0);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_CELL_SCALAR_FIELD<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::QUADTREE_CELL){
        int index=((OPENGL_SELECTION_QUADTREE_CELL<T>*)current_selection)->index;
        output_stream<<value(index);}
    output_stream<<std::endl;
}
//#####################################################################
template class OPENGL_QUADTREE_CELL_SCALAR_FIELD<float>;
template class OPENGL_QUADTREE_CELL_SCALAR_FIELD<float,int>;
template class OPENGL_QUADTREE_CELL_SCALAR_FIELD<float,bool>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_QUADTREE_CELL_SCALAR_FIELD<double>;
template class OPENGL_QUADTREE_CELL_SCALAR_FIELD<double,int>;
template class OPENGL_QUADTREE_CELL_SCALAR_FIELD<double,bool>;
#endif
#endif
