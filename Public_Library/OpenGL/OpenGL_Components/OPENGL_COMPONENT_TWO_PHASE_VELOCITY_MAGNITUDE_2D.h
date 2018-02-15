//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D__
#define __OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D__

#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class T> class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D;
template<class T> class OPENGL_COMPONENT_LEVELSET_2D;
template<class T> class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D;
template<class T> class OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D;

template<class T>
class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D(OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>& V_minus_component,OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>& V_plus_component,OPENGL_COMPONENT_LEVELSET_2D<T>& levelset_component);
    virtual ~OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D();

    bool Valid_Frame(int frame_input) const override;
    bool Is_Up_To_Date(int frame) const override { return valid; }

    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    void Toggle_3D_Mode();
    void Increase_Point_Size();
    void Decrease_Point_Size();

private:
    void Reinitialize(const bool force_even_if_not_drawn=false);
    bool valid;

private:
    T magnitude_height_scale;
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>& V_minus_component;
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>& V_plus_component;
    OPENGL_COMPONENT_LEVELSET_2D<T>& levelset_component;
public:
    OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>& opengl_two_phase_velocity_magnitude;
};

}

#endif
