//#####################################################################
// Copyright 2002-2005, Robert Bridson, Eran Guendelman, Eilene Hao, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_AXES
//##################################################################### 
#include <OpenGL/OpenGL/OPENGL_AXES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_AXES<T>::
OPENGL_AXES(const FRAME<TV>& frame_input,const RANGE<TV>& box_input,bool draw_box_input,
    bool draw_xz_grid_input,bool draw_xy_grid_input,bool draw_yz_grid_input,T grid_spacing_input)
    :box(box_input),draw_box(draw_box_input),draw_xz_grid(draw_xz_grid_input),draw_xy_grid(draw_xy_grid_input),
    draw_yz_grid(draw_yz_grid_input),grid_spacing(grid_spacing_input)
{
    *frame=frame_input;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_AXES<T>::
~OPENGL_AXES()
{
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_AXES<T>::
Display() const
{
    GLboolean lighting_enabled;
    glGetBooleanv(GL_LIGHTING,&lighting_enabled);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    glDisable(GL_LIGHTING);

    OpenGL_Begin(GL_LINES);

    glColor3f(1,.25,.25);
    // draw all lines parallel to x axis
    if(draw_box){
        OpenGL_Line(box.min_corner,TV(box.max_corner.x,box.min_corner.y,box.min_corner.z));
        OpenGL_Line(TV(box.min_corner.x,box.max_corner.y,box.min_corner.z),TV(box.max_corner.x,box.max_corner.y,box.min_corner.z));
        OpenGL_Line(TV(box.min_corner.x,box.max_corner.y,box.max_corner.z),TV(box.max_corner.x,box.max_corner.y,box.max_corner.z));
        OpenGL_Line(TV(box.min_corner.x,box.min_corner.y,box.max_corner.z),TV(box.max_corner.x,box.min_corner.y,box.max_corner.z));}
    else{
        OpenGL_Line(TV(box.min_corner.x,0,0),TV(box.max_corner.x,0,0));}
    if(draw_xz_grid)
        for(T z=box.min_corner.z;z<box.max_corner.z;z+=grid_spacing){
            OpenGL_Line(TV(box.min_corner.x,box.min_corner.y,z),TV(box.max_corner.x,box.min_corner.y,z));}
    if(draw_xy_grid)
        for(T y=box.min_corner.y;y<box.max_corner.y;y+=grid_spacing){
            OpenGL_Line(TV(box.min_corner.x,y,box.min_corner.z),TV(box.max_corner.x,y,box.min_corner.z));}

    OpenGL_End();
    OpenGL_Begin(GL_LINES);

    glColor3f(.25f,1,.25f);
    // draw all lines parallel to y axis
    if(draw_box){
        OpenGL_Line(TV(box.min_corner.x,box.min_corner.y,box.min_corner.z),TV(box.min_corner.x,box.max_corner.y,box.min_corner.z));
        OpenGL_Line(TV(box.max_corner.x,box.min_corner.y,box.min_corner.z),TV(box.max_corner.x,box.max_corner.y,box.min_corner.z));
        OpenGL_Line(TV(box.max_corner.x,box.min_corner.y,box.max_corner.z),TV(box.max_corner.x,box.max_corner.y,box.max_corner.z));
        OpenGL_Line(TV(box.min_corner.x,box.min_corner.y,box.max_corner.z),TV(box.min_corner.x,box.max_corner.y,box.max_corner.z));}
    else{
        OpenGL_Line(TV(0,box.min_corner.y,0),TV(0,box.max_corner.y,0));}
    if(draw_yz_grid)
        for(T z=box.min_corner.z;z<box.max_corner.z;z+=grid_spacing){
            OpenGL_Line(TV(box.min_corner.x,box.min_corner.y,z),TV(box.min_corner.x,box.max_corner.y,z));}
    if(draw_xy_grid)
        for(T x=box.min_corner.x;x<box.max_corner.x;x+=grid_spacing){
            OpenGL_Line(TV(x,box.min_corner.y,box.min_corner.z),TV(x,box.max_corner.y,box.min_corner.z));}
    
    OpenGL_End();
    OpenGL_Begin(GL_LINES);

    glColor3f(.25f,.25f,1);
    // draw all lines parallel to z axis
    if(draw_box){
        OpenGL_Line(TV(box.min_corner.x,box.min_corner.y,box.min_corner.z),TV(box.min_corner.x,box.min_corner.y,box.max_corner.z));
        OpenGL_Line(TV(box.max_corner.x,box.min_corner.y,box.min_corner.z),TV(box.max_corner.x,box.min_corner.y,box.max_corner.z));
        OpenGL_Line(TV(box.max_corner.x,box.max_corner.y,box.min_corner.z),TV(box.max_corner.x,box.max_corner.y,box.max_corner.z));
        OpenGL_Line(TV(box.min_corner.x,box.max_corner.y,box.min_corner.z),TV(box.min_corner.x,box.max_corner.y,box.max_corner.z));}
    else{
        OpenGL_Line(TV(0,0,box.min_corner.z),TV(0,0,box.max_corner.z));}
    if(draw_yz_grid)
        for(T y=box.min_corner.y;y<box.max_corner.y;y+=grid_spacing){
            OpenGL_Line(TV(box.min_corner.x,y,box.min_corner.z),TV(box.min_corner.x,y,box.max_corner.z));}
    if(draw_xz_grid)
        for(T x=box.min_corner.x;x<box.max_corner.x;x+=grid_spacing){
            OpenGL_Line(TV(x,box.min_corner.y,box.min_corner.z),TV(x,box.min_corner.y,box.max_corner.z));}

    OpenGL_End();
    glPopMatrix();
    if(lighting_enabled) glEnable(GL_LIGHTING);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_AXES<T>::
Bounding_Box() const
{
    return World_Space_Box(box);
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_AXES<float>;
template class OPENGL_AXES<double>;
}
