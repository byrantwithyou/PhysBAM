//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_GRID_1D.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_GRID_1D<T>::
Display() const
{
    if(!draw)return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    color.Send_To_GL_Pipeline();

    int ghost_cells=draw_ghost_values?3:0;

    // Draw masks
    VECTOR<T,1> min_corner(grid.Node(TV_INT(1-ghost_cells)));
    VECTOR<T,1> max_corner(grid.Node(TV_INT(grid.numbers_of_cells.x+ghost_cells+1)));

    // Outline boundary of real domain
    if(ghost_cells>0){
        glPushAttrib(GL_ALL_ATTRIB_BITS);OPENGL_COLOR(.07,.12,.42).Send_To_GL_Pipeline();
        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);
        OpenGL_Begin(GL_TRIANGLE_STRIP);
        OpenGL_Triangle_Strip_2D(VECTOR<T,2>(grid.domain.min_corner.x,(T)-1000.),VECTOR<T,2>(grid.domain.max_corner.x,(T)1000.));
        OpenGL_End();
        glPopAttrib();}

    // Draw grid
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glPointSize(3.0);
    OpenGL_Begin(GL_POINTS);
    for(int i=-ghost_cells;i<grid.counts.x+ghost_cells;i++) OpenGL_Vertex(VECTOR<T,3>(grid.X(TV_INT(i)).x,(T)0,(T)0));
    OpenGL_End();
    glPopAttrib();

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_GRID_1D<T>::
Set_Frame()
{
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_GRID_1D<T>::
Bounding_Box() const
{
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Function Toggle_Draw_Ghost_Values
//#####################################################################
template<class T> void OPENGL_GRID_1D<T>::
Toggle_Draw_Ghost_Values()
{
    draw_ghost_values=!draw_ghost_values;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_GRID_1D<float>;
template class OPENGL_GRID_1D<double>;
}
