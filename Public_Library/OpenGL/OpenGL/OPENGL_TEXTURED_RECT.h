//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TEXTURED_RECT
//#####################################################################
#ifndef __OPENGL_TEXTURED_RECT__
#define __OPENGL_TEXTURED_RECT__

#include <Core/Math_Tools/RANGE.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_TEXTURE.h>

namespace PhysBAM
{
template<class T>
class OPENGL_TEXTURED_RECT:public OPENGL_OBJECT<T>
{
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    double width, height;
    OPENGL_TEXTURE *texture;

    OPENGL_TEXTURED_RECT(STREAM_TYPE stream_type);

    void Set_Texture(OPENGL_TEXTURE *texture_input);
    void Display() const override;
    RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

}

#endif
