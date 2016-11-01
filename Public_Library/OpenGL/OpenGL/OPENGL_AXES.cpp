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
OPENGL_AXES(STREAM_TYPE stream_type,const FRAME<TV>& frame_input,const RANGE<TV>& box_input)
    :OPENGL_OBJECT<T>(stream_type),box(box_input),draw_box(false),number_grid_spaces(10)
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

    TV draw_grid_dX=box.Edge_Lengths()/number_grid_spaces;
    for(int a=0;a<TV::m;a++){
        glColor3f(a==0?1:.25,a==1?1:.25,a==2?1:.25);
        // draw all lines parallel to x axis
        if(draw_box){
            for(int i=0;i<(1<<TV::m);i++)
                if(!(i&(1<<a))){
                    TV start;
                    for(int j=0;j<TV::m;j++)
                        start(j)=(i&(1<<j))?box.min_corner(j):box.max_corner(j);
                    TV end=start;
                    end(a)=box.max_corner(a);
                    OpenGL_Line(start,end);}
        else{
            TV start=box.min_corner,end=box.min_corner;
            end(a)=box.max_corner(a);
            OpenGL_Line(start,end);}}
        for(int b=0;b<TV::m;b++){
            if(a==b) continue;
            int dg=TV::m*2-3-a-b;
            if(!draw_grid(dg)) continue;
            T d=draw_grid_dX(b);
            for(int i=0;i<=number_grid_spaces;i++){
                TV A=box.min_corner,B=box.min_corner;
                B(a)=box.max_corner(a);
                A(b)+=i*d;
                B(b)+=i*d;
                OpenGL_Line(A,B);}}}

    OpenGL_End();
    glPopMatrix();
    if(lighting_enabled) glEnable(GL_LIGHTING);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> auto OPENGL_AXES<T>::
Bounding_Box() const -> RANGE<TV>
{
    return World_Space_Box(box);
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_AXES<float>;
template class OPENGL_AXES<double>;
}
