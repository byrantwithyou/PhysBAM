//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_TRIANGULATED_AREA
//##################################################################### 
#ifndef __OPENGL_COMPONENT_TRIANGULATED_AREA__
#define __OPENGL_COMPONENT_TRIANGULATED_AREA__

#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_TRIANGULATED_AREA:public OPENGL_COMPONENT<T>
{
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::is_animation;
    OPENGL_COMPONENT_TRIANGULATED_AREA(const std::string &filename);
    OPENGL_COMPONENT_TRIANGULATED_AREA(const std::string &filename,const std::string &color_map_filename);
    virtual ~OPENGL_COMPONENT_TRIANGULATED_AREA();
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display() const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

private:
    void Reinitialize();    // Needs to be called after some state changes
    ARRAY<OPENGL_COLOR >* color_map;

public:
    TRIANGULATED_AREA<T>& triangulated_area;
    OPENGL_TRIANGULATED_AREA<T> opengl_triangulated_area;

private:
    std::string filename;
    const std::string *color_map_filename;
    int frame_loaded;
    bool valid;
};

}

#endif
