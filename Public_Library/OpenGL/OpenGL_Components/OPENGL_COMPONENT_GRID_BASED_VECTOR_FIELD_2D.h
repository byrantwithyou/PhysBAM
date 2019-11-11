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
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;
    using OPENGL_OBJECT<T>::viewer_callbacks;using OPENGL_COMPONENT<T>::viewer_dir;
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid,const std::string &vector_field_filename_input);
    virtual ~OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D();

    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();

private:
    void Reinitialize();
    template<class> friend class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D;
public:
    OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>& opengl_grid_based_vector_field;

private:
    std::string vector_field_filename;
    bool valid;
};

}

#endif
