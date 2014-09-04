//#####################################################################
// Copyright 2003, 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_BOX_3D
//##################################################################### 
#ifndef __OPENGL_BOX_3D__
#define __OPENGL_BOX_3D__

#include <Tools/Math_Tools/RANGE.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{

template<class T>
class OPENGL_BOX_3D:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    RANGE<TV>& box;
    OPENGL_COLOR color;

    OPENGL_BOX_3D(RANGE<TV>& box_input,const OPENGL_COLOR& color_input=OPENGL_COLOR::White()) 
        :box(box_input),color(color_input)
    {}

    void Display() const PHYSBAM_OVERRIDE;
    virtual RANGE<TV> Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
