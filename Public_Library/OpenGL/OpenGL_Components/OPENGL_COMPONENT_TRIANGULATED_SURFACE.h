//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_TRIANGULATED_SURFACE
//##################################################################### 
#ifndef __OPENGL_COMPONENT_TRIANGULATED_SURFACE__
#define __OPENGL_COMPONENT_TRIANGULATED_SURFACE__

#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T>
class OPENGL_COMPONENT_TRIANGULATED_SURFACE:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::is_animation;using OPENGL_COMPONENT<T>::stream_type;
    OPENGL_COMPONENT_TRIANGULATED_SURFACE(STREAM_TYPE stream_type,const std::string &filename, bool use_display_list = true);
    virtual ~OPENGL_COMPONENT_TRIANGULATED_SURFACE();
    
    bool Valid_Frame(int frame_input) const override;
    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<TV> Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override
    {return opengl_triangulated_surface.Get_Selection_Priority(indices);}

    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override
    {return opengl_triangulated_surface.Set_Selection(indices,modifiers);}

    virtual void Clear_Selection() override { opengl_triangulated_surface.Clear_Selection(); }

    virtual RANGE<TV> Selection_Bounding_Box() const override;
    
private:
    void Reinitialize();    // Needs to be called after some state changes

public:
    TRIANGULATED_SURFACE<T>& triangulated_surface;
    OPENGL_TRIANGULATED_SURFACE<T> opengl_triangulated_surface;

private:
    std::string filename;
    int frame_loaded;
    bool valid;
    bool use_display_list;
};

}

#endif
