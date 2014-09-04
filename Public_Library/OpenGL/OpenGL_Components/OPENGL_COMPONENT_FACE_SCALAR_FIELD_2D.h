//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D__
#define __OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D__
#include <OpenGL/OpenGL/OPENGL_FACE_SCALAR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>

namespace PhysBAM{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D:public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
public:
    OPENGL_FACE_SCALAR_FIELD_2D<T,T2>  opengl_scalar_field;
private:
    std::string values_filename;
    std::string x_face_values_filename,y_face_values_filename;
    int frame_loaded;
    bool valid;

public:
    OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D(const GRID<TV> &grid_input,const std::string &values_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input);
    OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D(const GRID<TV> &grid_input,const std::string &x_face_values_filename_input,const std::string &y_face_values_filename_input,
        OPENGL_COLOR_MAP<T2>* color_map_input);
    virtual ~OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D();

    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE
    {return valid && frame_loaded==frame;}

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return draw && valid;}

//#####################################################################
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE {slice=slice_input;opengl_scalar_field.Set_Slice(slice_input);}
    virtual void Slice_Has_Changed() PHYSBAM_OVERRIDE {opengl_scalar_field.Slice_Has_Changed();}
private:
    void Reinitialize();
//#####################################################################
};
}
#endif
