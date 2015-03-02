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
public:
    OPENGL_FACE_SCALAR_FIELD_1D<T,T2> opengl_face_scalar_field;
private:
    std::string values_filename;
    int frame_loaded;
    bool valid;

public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;
    OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D(STREAM_TYPE stream_type,const GRID<TV> &grid_input,const std::string &values_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color);
    virtual ~OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D();

    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE
    {return valid && frame_loaded==frame;}

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return draw && valid;}

//#####################################################################
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE {slice=slice_input;opengl_face_scalar_field.Set_Slice(slice_input);}
    virtual void Slice_Has_Changed() PHYSBAM_OVERRIDE {opengl_face_scalar_field.Slice_Has_Changed();}
    void Scale(const T scale);
    void Increase_Scale();
    void Decrease_Scale();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D,Increase_Scale,"Increase scale");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D,Decrease_Scale,"Decrease scale");
private:
    void Reinitialize();
//#####################################################################
};
}
#endif
