//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D__
#define __OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D__
#include <OpenGL/OpenGL/OPENGL_FACE_SCALAR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>
namespace PhysBAM{

template<class T,class T2=T>
class OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<T2,FACE_INDEX<3> > internal_scalar_field;
    OPENGL_FACE_SCALAR_FIELD_3D<T,T2>  opengl_scalar_field;
private:
    std::string values_filename;
    int frame_loaded;
    bool valid;

public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;
    OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D(STREAM_TYPE stream_type,const GRID<TV> &grid_input,
        const std::string &values_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input);
    virtual ~OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D();

    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE
    {slice=slice_input;opengl_scalar_field.Set_Slice(slice_input);}

    virtual void Slice_Has_Changed()
    {if(draw) opengl_scalar_field.Slice_Has_Changed();}

    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE
    {return valid && frame_loaded==frame;}

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return draw && valid;}

//##################################################################### 
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
private:
    void Reinitialize();
//##################################################################### 
};
}
#endif
