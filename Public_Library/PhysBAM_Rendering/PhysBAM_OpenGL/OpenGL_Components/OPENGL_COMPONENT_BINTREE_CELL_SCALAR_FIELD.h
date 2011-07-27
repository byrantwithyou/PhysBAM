//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD__
#define __OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BINTREE_CELL_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD:public OPENGL_COMPONENT
{
public:
    OPENGL_BINTREE_CELL_SCALAR_FIELD<T,T2> opengl_bintree_cell_scalar_field;

private:
    std::string scalar_field_filename;
    int frame_loaded;
    bool valid;
public:
    OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD(BINTREE_GRID<T> &grid,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map,const bool draw_as_points);
    ~OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;

    void Increase_Point_Size();
    void Decrease_Point_Size();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD, Increase_Point_Size, "Increase point size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD, Decrease_Point_Size, "Decrease point size");

private:
    void Reinitialize();

};

}

#endif
#endif
