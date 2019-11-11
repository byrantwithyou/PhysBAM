//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D__
#define __OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D__

#include <OpenGL/OpenGL/OPENGL_GRID_BASED_VECTOR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T>
class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::viewer_dir;
    using OPENGL_COMPONENT<T>::World_Space_Box;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid,const std::string &vector_field_filename);
    virtual ~OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D();

    void Display() const override;
    virtual RANGE<TV> Bounding_Box() const override;
    virtual void Set_Slice(OPENGL_SLICE *slice_input) override { slice=slice_input; opengl_grid_based_vector_field.Set_Slice(slice_input); }
    virtual void Slice_Has_Changed() override { opengl_grid_based_vector_field.Slice_Has_Changed(); }
    void Print_Selection_Info(std::ostream& stream) const override;

    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();

private:
    void Reinitialize();

public:
    OPENGL_GRID_BASED_VECTOR_FIELD_3D<T> opengl_grid_based_vector_field;

private:
    std::string vector_field_filename;
    bool valid;
};

}

#endif
