//#####################################################################
// Copyright 2005, Jared Go, Ranjitha Kumar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <sstream>
#include <string>
#include "../../Public_Library/Articulated_Rigid_Bodies/BEND_JOINT.h"
#include "ARTICULATION_VISUALIZATION.h"
#include "OPENGL_COMPONENT_JOINT_3D.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_JOINT_3D<T,RW>::
OPENGL_COMPONENT_JOINT_3D(ARTICULATION_VISUALIZATION<T> *o, JOINT_3D<T>*j)
    :OPENGL_COMPONENT("Muscle"),owner(o),joint(j),selected(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_JOINT_3D<T,RW>::
~OPENGL_COMPONENT_JOINT_3D()
{
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_JOINT_3D<T,RW>::
Display(const int in_color) const
{
    //glPushName(i); 
    FRAME_3D<T> jf = FRAME_3D<T>(*(parent_surf->opengl_triangulated_surface.frame)) * joint->frame_pj;
    OPENGL_SHAPES::Draw_Dot(jf.t,OPENGL_COLOR((T)0.25,(T)0.33,1), 6);

    T len = (T)0.01;
    BOX_3D<T> axes_box(0,len,0,len,0,len);    
    (OPENGL_AXES<T>(jf,axes_box)).Display();

    //glPopName();        
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> BOX_3D<float> OPENGL_COMPONENT_JOINT_3D<T,RW>::
Bounding_Box() const
{
    return BOX_3D<float>();
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_JOINT_3D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{ 
    return NULL;
}

//#####################################################################
template class OPENGL_COMPONENT_JOINT_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_JOINT_3D<double,double>;
#endif
