//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_MPM_PARTICLES_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_MPM_PARTICLES_3D__
#define __OPENGL_COMPONENT_MPM_PARTICLES_3D__

#include <Core/Arrays/ARRAY.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class MPM_PARTICLES;
template<class T>
class OPENGL_COMPONENT_MPM_PARTICLES_3D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;
    using OPENGL_COMPONENT<T>::Slice_Has_Changed;using OPENGL_COMPONENT<T>::World_Space_Box;
    using OPENGL_OBJECT<T>::viewer_callbacks;using OPENGL_COMPONENT<T>::viewer_dir;
    using OPENGL_COMPONENT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::Set_Slice;
    OPENGL_COMPONENT_MPM_PARTICLES_3D(const VIEWER_DIR& viewer_dir,const std::string &filename);
    virtual ~OPENGL_COMPONENT_MPM_PARTICLES_3D();

    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override;
    virtual RANGE<TV> Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    void Print_Selection_Info(std::ostream &output_stream) const override;
    virtual RANGE<TV> Selection_Bounding_Box() const override;
    bool Destroy_Selection_After_Frame_Change() override;
    void Toggle_Draw_Velocities();
    void Toggle_Draw_Phases();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();
    void Set_Slice(OPENGL_SLICE *slice_input) override;

private:
    OPENGL_INDEXED_COLOR_MAP* color_map;
    void Reinitialize();

public:
    MPM_PARTICLES<TV>& particles;
    OPENGL_COLOR default_color;
    OPENGL_COLOR velocity_color;
    bool draw_velocities;
    bool draw_phases;
    bool draw_arrows;
    bool draw_B;
    bool draw_F;
    VECTOR<OPENGL_COLOR,TV::m> B_color;
    VECTOR<OPENGL_COLOR,TV::m> F_color;
    T scale_velocities;

private:
    std::string filename;
    bool valid;
    int selected_index;
};

}

#endif
