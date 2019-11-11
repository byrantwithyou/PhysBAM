//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D__
#define __OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D__
#include <OpenGL/OpenGL/OPENGL_FACE_SCALAR_FIELD_1D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>

namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,1> TV;
    typedef VECTOR<T,TV::m+(TV::m==1)> TV_BOX;
public:
    OPENGL_FACE_SCALAR_FIELD_1D<T,T2> opengl_face_scalar_field;
private:
    std::string values_filename;
    bool valid;

public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;
    using OPENGL_COMPONENT<T>::viewer_dir;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid_input,const std::string &values_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color);
    virtual ~OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D();

    bool Use_Bounding_Box() const override
    {return draw && valid;}

//#####################################################################
    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;
    void Display() const override;
    void Print_Selection_Info(std::ostream& stream) const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    virtual void Set_Slice(OPENGL_SLICE *slice_input) override {slice=slice_input;opengl_face_scalar_field.Set_Slice(slice_input);}
    virtual void Slice_Has_Changed() override {opengl_face_scalar_field.Slice_Has_Changed();}
    void Scale(const T scale);
    void Increase_Scale();
    void Decrease_Scale();
private:
    void Reinitialize();
//#####################################################################
};
}
#endif
