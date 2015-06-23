//#####################################################################
// Copyright 2006-2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SCALAR_FIELD_1D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_SCALAR_FIELD_1D__
#define __OPENGL_COMPONENT_SCALAR_FIELD_1D__

#include <OpenGL/OpenGL/OPENGL_SCALAR_FIELD_1D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>
namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_COMPONENT_SCALAR_FIELD_1D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,1> TV;
private:
    std::string scalar_field_filename;
    int frame_loaded;
    bool valid;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;using OPENGL_COMPONENT<T>::viewer_callbacks;
    OPENGL_SCALAR_FIELD_1D<T,T2>  opengl_scalar_field;

    bool Is_Up_To_Date(int frame) const override
    {return valid && frame_loaded==frame;}

    bool Use_Bounding_Box() const  override
    {return draw && valid;}

//##################################################################### 
    OPENGL_COMPONENT_SCALAR_FIELD_1D(STREAM_TYPE stream_type,const GRID<TV> &grid_input,const std::string &scalar_field_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color);
    virtual ~OPENGL_COMPONENT_SCALAR_FIELD_1D();
    bool Valid_Frame(int frame_input) const override;
    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input=true) override;
    void Display() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    void Scale(const T scale);
    void Increase_Scale();
    void Decrease_Scale();
private:
    void Reinitialize();
//##################################################################### 
};
}
#endif
