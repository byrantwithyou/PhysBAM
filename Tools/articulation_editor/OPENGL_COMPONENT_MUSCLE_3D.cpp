//#####################################################################
// Copyright 2005, Jared Go, Ranjitha Kumar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <sstream>
#include <string>
#include "ARTICULATION_VISUALIZATION.h"
#include "OPENGL_COMPONENT_MUSCLE_3D.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_MUSCLE_3D<T,RW>::
OPENGL_COMPONENT_MUSCLE_3D(ARTICULATION_VISUALIZATION<T> *o)
    :OPENGL_COMPONENT("Muscle"),owner(o),selected(false),targpoint(-1),segment(-1)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_MUSCLE_3D<T,RW>::
~OPENGL_COMPONENT_MUSCLE_3D()
{
}
//#####################################################################
// Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_MUSCLE_3D<T,RW>::
AddPointInFrame(const VECTOR<T,3> &pt, const FRAME_3D<T> &frame, const std::string &bname, const std::string &a, const std::string &st)
{
    points.Append(pt);
    frames.Append(frame);
    bones.Append(bname);
    attach.Append(a);
    // first via_point has no corresponding segment type
    segment_type.Append(st);
}

//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_MUSCLE_3D<T,RW>::
Display(const int in_color) const
{

    glPushAttrib(GL_LINE_BIT|GL_ENABLE_BIT|GL_CURRENT_BIT);

    OPENGL_COLOR muscle_color;
    glLineWidth(selected?2.0f:0.5f);    
    glDisable(GL_LIGHTING);

    glPushName(0);
    for (int i = 2; i <= points.m; i++) {
        if (segment_type(i)=="linear")
            muscle_color = selected ?  OPENGL_COLOR(0.05f,0.05f,0.9f) :OPENGL_COLOR(0.8f,0.05f,0.05f);
        else if (segment_type(i)=="analytic")
            muscle_color = selected ?  OPENGL_COLOR(0.05f,0.05f,0.9f) :OPENGL_COLOR(0.05f,0.8f,0.05f);
        muscle_color.Send_To_GL_Pipeline();
        glPushName(i); glBegin(GL_LINE_STRIP); OpenGL_Vertex(frames(i-1) * points(i-1));OpenGL_Vertex(frames(i) * points(i));glEnd();
        glPopName();}
    glPopName();

    // highlight selected muscle segment
    if(selected) {
        glPushAttrib(GL_DEPTH_BUFFER_BIT);glDepthFunc(GL_LEQUAL);
        glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);OPENGL_PREFERENCES::selection_highlight_color.Send_To_GL_Pipeline();
        glBegin(GL_LINE_STRIP); OpenGL_Vertex(frames(segment-1)*points(segment-1));OpenGL_Vertex(frames(segment)*points(segment));glEnd();
        glPopAttrib();
    }

    
    glPushName(1);
    glPushAttrib(GL_POINT_BIT);
    glPointSize(OPENGL_PREFERENCES::selection_point_size);    
    for (int i = 1; i <= points.m; i++) 
    {
        glPushName(i);
        OPENGL_SHAPES::Draw_Dot(frames(i) * points(i),OPENGL_COLOR((T)0.25,(T)0.33,1), selected?5:2);
        glPopName();
    }
    glPopAttrib();
    glPopName();

    
    glPopAttrib();

    if (selected)
    {
        glPushName(2);
        for (int i = 1; i <= points.m; i++) 
        {
            glPushName(i);    
            T len = (T)0.02;
            BOX_3D<T> axes_box(0,len,0,len,0,len);
            VECTOR<T,3> worldpos = frames(i) * points(i);

            (OPENGL_AXES<T>(FRAME_3D<T>(worldpos),axes_box)).Display();

            glPopName();    
        }
        glPopName();    

        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glColor3f(1,1,1);     
        for (int i = 1; i <= points.m; i++) 
        {
            VECTOR<T,3> worldpos = frames(i) * points(i);
            char buf[16];
            sprintf(buf, "%d",  i);
            OpenGL_String(worldpos, this->name + buf + "(" + bones(i) + ")");
        }
        glPopAttrib();
    
    }
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> BOX_3D<float> OPENGL_COMPONENT_MUSCLE_3D<T,RW>::
Bounding_Box() const
{
    return BOX_3D<float>();
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_MUSCLE_3D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    
    switch (buffer[0])
    {
    case 0:
        //muscles
        std::cout<<"Clicked muscle "<<this->name<<".\n";
        targpoint=-1;
        segment=buffer[1];
        owner->Handle_Object_Selected(this, 1);
        return new OPENGL_SELECTION((OPENGL_SELECTION::TYPE)0x1000, this);
        break;
    case 1:
        //via points
        std::cout<<"Clicked via point: "<<buffer[1]<<" on muscle "<<this->name<<".\n";
        targpoint=buffer[1];
        segment=targpoint;
        if(segment==1)
            segment=2;
        owner->Handle_Object_Selected(this, buffer[1]);
        return new OPENGL_SELECTION((OPENGL_SELECTION::TYPE)0x1001, this);
        break;
    case 2:
        //via points
        std::cout<<"Clicked axes for point: "<<buffer[1]<<" on muscle "<<this->name<<".\n";
        targpoint=buffer[1];
        segment=targpoint;
        if(segment==1)
            segment=2;
        owner->Handle_Object_Selected(this, buffer[1]);
        return new OPENGL_SELECTION((OPENGL_SELECTION::TYPE)0x1002, this);
        break;
    default:
        break;
    }    
    return NULL;
}
   
//#####################################################################
template class OPENGL_COMPONENT_MUSCLE_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_MUSCLE_3D<double,double>;
#endif
