//#####################################################################
// Copyright 2005, Jared Go, Ranjitha Kumar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_JOINT_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_JOINT_3D__
#define __OPENGL_COMPONENT_JOINT_3D__

#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_SURFACE.h>
#include <string>

namespace PhysBAM
{

template<class T> class ARTICULATION_VISUALIZATION;

template<class T,class RW=T>
class OPENGL_COMPONENT_JOINT_3D : public OPENGL_COMPONENT
{
public:
    OPENGL_COMPONENT_JOINT_3D(ARTICULATION_VISUALIZATION<T> *owner, JOINT_3D<T>*j);
    virtual ~OPENGL_COMPONENT_JOINT_3D();
    
    // Owning vis
    ARTICULATION_VISUALIZATION<T> *owner;

    // The joint 3d we represent
    JOINT_3D<T> *joint;

    // Should we be hilighed?
    bool selected;    

    // Calls
    virtual void Display(const int in_color=1) const;
    virtual bool Use_Bounding_Box() const { return false; }
    virtual BOX_3D<float> Bounding_Box() const;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    
private:
    
public:
    OPENGL_COMPONENT_TRIANGULATED_SURFACE<T> *parent_surf;
    OPENGL_COMPONENT_TRIANGULATED_SURFACE<T> *child_surf;

private:    
};

}

#endif
