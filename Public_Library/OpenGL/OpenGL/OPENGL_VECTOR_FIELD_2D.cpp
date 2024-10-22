//#####################################################################
// Copyright 2003, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_ARRAY> OPENGL_VECTOR_FIELD_2D<T_ARRAY>::
OPENGL_VECTOR_FIELD_2D(T_ARRAY& vector_field,T_ARRAY& vector_locations,const OPENGL_COLOR &color,double size,bool draw_arrowhead,bool draw_value)
    :vector_field(vector_field),vector_locations(vector_locations),vector_color(color),size(size),draw_arrowhead(draw_arrowhead),draw_value(draw_value),draw(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_ARRAY> OPENGL_VECTOR_FIELD_2D<T_ARRAY>::
~OPENGL_VECTOR_FIELD_2D()
{
}
//#####################################################################
// Function Display
//#####################################################################
template<class T_ARRAY> void OPENGL_VECTOR_FIELD_2D<T_ARRAY>::
Display() const
{
    if(!draw)return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    
    glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_TEXTURE_BIT);
    glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);glDisable(GL_DEPTH_TEST);

    vector_color.Send_To_GL_Pipeline();

    OpenGL_Begin(GL_LINES);
    if(draw_arrowhead) 
        for(int i=0;i<vector_locations.Size();i++){OPENGL_SHAPES::Draw_Arrow(vector_locations(i),vector_locations(i)+(T)size*vector_field(i));}
    else 
        for(int i=0;i<vector_locations.Size();i++){OpenGL_Line(vector_locations(i),vector_locations(i)+(T)size*vector_field(i));}
    OpenGL_End();

    if(draw_value){
        (vector_color+OPENGL_COLOR(.8f,.8f,.8f)).Send_To_GL_Pipeline();
        for(int i=0;i<vector_locations.Size();i++)
            OpenGL_String(vector_locations(i)+(T)1.1*(T)size*vector_field(i),LOG::sprintf("%.3f %.3f",vector_field(i).x,vector_field(i).y));}

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T_ARRAY> RANGE<VECTOR<typename T_ARRAY::SCALAR,3> > OPENGL_VECTOR_FIELD_2D<T_ARRAY>::
Bounding_Box() const
{
    RANGE<VECTOR<T,3> > box(World_Space_Point(vector_locations(0)));
    for(int i=1;i<vector_locations.Size();i++)
        box.Enlarge_Nonempty_Box_To_Include_Point(World_Space_Point(vector_locations(i)));
    return box;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T_ARRAY> void OPENGL_VECTOR_FIELD_2D<T_ARRAY>::
Scale_Vector_Size(const T scale)
{
    size*=scale;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_VECTOR_FIELD_2D<ARRAY<VECTOR<float,2> > >;
template class OPENGL_VECTOR_FIELD_2D<ARRAY_VIEW<VECTOR<float,2> > >;
template class OPENGL_VECTOR_FIELD_2D<ARRAY<VECTOR<double,2> > >;
template class OPENGL_VECTOR_FIELD_2D<ARRAY_VIEW<VECTOR<double,2> > >;
}
