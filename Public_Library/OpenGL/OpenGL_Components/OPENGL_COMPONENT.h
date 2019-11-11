//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT
//#####################################################################
#ifndef __OPENGL_COMPONENT__
#define __OPENGL_COMPONENT__

#include <OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <string>

namespace PhysBAM
{
class VIEWER_DIR;

template <class T>
class OPENGL_COMPONENT:public OPENGL_OBJECT<T>
{
public:
    using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT(const VIEWER_DIR& viewer_dir,const std::string &name = "");
    virtual ~OPENGL_COMPONENT();

    void Set_Name(const std::string &name) { component_name = name; }

    bool Use_Bounding_Box() const override;

    virtual void Set_Frame()=0;
    virtual void Set_Draw(bool draw_input = true);
    virtual void Draw_All_Objects();

    const VIEWER_DIR& viewer_dir;
    bool draw;
    std::string component_name;
};

}

#endif
