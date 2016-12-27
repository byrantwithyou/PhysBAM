//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D__
#define __OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D__

#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

template<class T> class OPENGL_GRID_BASED_VECTOR_FIELD_2D;

namespace PhysBAM
{
template<class TV> class GRID;

template<class T>
class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::component_name;
    using OPENGL_COMPONENT<T>::is_animation;using OPENGL_COMPONENT<T>::stream_type;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D(STREAM_TYPE stream_type,const GRID<TV> &grid,const std::string &vector_field_filename_input);
    virtual ~OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D();

    bool Valid_Frame(int frame_input) const override;
    bool Is_Up_To_Date(int frame) const override { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();

private:
    void Reinitialize(bool force_load_even_if_not_drawn=false);
    template<class> friend class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D;
public:
    OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>& opengl_grid_based_vector_field;

private:
    std::string vector_field_filename;
    int frame_loaded;
    bool valid;
};

}

#endif
