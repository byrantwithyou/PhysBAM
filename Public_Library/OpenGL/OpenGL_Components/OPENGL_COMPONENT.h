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

template <class T>
class OPENGL_COMPONENT:public OPENGL_OBJECT<T>
{
public:
    using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT(STREAM_TYPE stream_type,const std::string &name = "");
    virtual ~OPENGL_COMPONENT();

    void Set_Name(const std::string &name) { component_name = name; }

    bool Use_Bounding_Box() const override;

    virtual bool Valid_Frame(int frame_input) const;
    virtual bool Is_Up_To_Date(int frame) const;

    virtual void Set_Frame(int frame_input);
    virtual void Set_Draw(bool draw_input = true);
    virtual void Draw_All_Objects();
    bool Is_Animated() { return is_animation; }

    int     frame;
    bool    draw;
    bool    is_animation;

    std::string component_name;
};

}

#endif
