//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D__
#define __OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D__

#include <OpenGL/OpenGL/OPENGL_MAC_VELOCITY_FIELD_3D.h>
#include <OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T>
class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,TV::m> TV_INT;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;
    using OPENGL_OBJECT<T>::viewer_callbacks;using OPENGL_COMPONENT<T>::viewer_dir;
    OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid,const std::string &velocity_filename_input);
    virtual ~OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D();

    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    virtual void Set_Slice(OPENGL_SLICE *slice_input) override {slice=slice_input;opengl_mac_velocity_field.Set_Slice(slice_input);opengl_vorticity_magnitude.Set_Slice(slice_input);}
    virtual void Slice_Has_Changed() override { opengl_mac_velocity_field.Slice_Has_Changed();opengl_vorticity_magnitude.Slice_Has_Changed(); }
    void Print_Selection_Info(std::ostream& stream) const override;

    void Toggle_Velocity_Mode();
    void Toggle_Velocity_Mode_And_Draw();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();
    void Toggle_Draw_Vorticity();
    void Normalize_Vorticity_Color_Map();

private:
    void Reinitialize();
    void Update_Vorticity();

public:
    OPENGL_MAC_VELOCITY_FIELD_3D<T> opengl_mac_velocity_field;
    ARRAY<T,VECTOR<int,3> > opengl_vorticity_magnitude_array;
    OPENGL_SCALAR_FIELD_3D<T> opengl_vorticity_magnitude;
    bool draw_vorticity;
    TV_INT selected_index;

private:
    std::string velocity_filename;
    bool valid;
    T min_vorticity,max_vorticity;
};

}

#endif
