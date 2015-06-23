//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_1D
//#####################################################################
#ifndef __OPENGL_COMPONENT_LEVELSET_1D__
#define __OPENGL_COMPONENT_LEVELSET_1D__

#include <OpenGL/OpenGL/OPENGL_LEVELSET_1D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T>
class OPENGL_COMPONENT_LEVELSET_1D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,1> TV;typedef VECTOR<int,TV::m> TV_INT;
private:
    std::string levelset_filename;
    int frame_loaded;
    bool valid;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;
    OPENGL_LEVELSET_1D<T>* opengl_levelset;

//##################################################################### 
    OPENGL_COMPONENT_LEVELSET_1D(STREAM_TYPE stream_type,GRID<TV> &grid,const std::string& levelset_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color);
    virtual ~OPENGL_COMPONENT_LEVELSET_1D();
    bool Valid_Frame(int frame_input) const override;
    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input=true) override;
    void Display() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
private:
    void Reinitialize();
//##################################################################### 
};
}
#endif
