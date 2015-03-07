//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D__
#define __OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D__

#include <OpenGL/OpenGL/OPENGL_MAC_VELOCITY_FIELD_2D.h>
#include <OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T>
class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,TV::m> TV_INT;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::Toggle_Draw;using OPENGL_COMPONENT<T>::stream_type;
    typedef LINEAR_INTERPOLATION_UNIFORM<TV,T> T_LINEAR_INTERPOLATION_VECTOR;
    OPENGL_MAC_VELOCITY_FIELD_2D<T>* opengl_mac_velocity_field;
    OPENGL_SCALAR_FIELD_2D<T>* opengl_vorticity_magnitude;
    bool draw_vorticity;
private:
    std::string velocity_filename;
    int frame_loaded;
    bool valid;
    bool draw_divergence;
    bool draw_streamlines,use_seed_for_streamlines;
    ARRAY<T,VECTOR<int,2> > divergence;
    OPENGL_SCALAR_FIELD_2D<T>* opengl_divergence_field;
    SEGMENTED_CURVE_2D<T> streamlines;
    OPENGL_SEGMENTED_CURVE_2D<T> opengl_streamlines;
    std::string psi_N_psi_D_basedir;
    int number_of_steps;
    T min_vorticity,max_vorticity;
    unsigned int streamline_seed;

public:
    OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D(STREAM_TYPE stream_type,const GRID<TV> &grid,const std::string &velocity_filename_input);
    virtual ~OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Set_Psi_N_Psi_D_Basedir_For_Divergence(std::string psi_N_psi_D_basedir_input)
    {psi_N_psi_D_basedir=psi_N_psi_D_basedir_input;}

    void Toggle_Velocity_Mode();
    void Toggle_Velocity_Mode_And_Draw();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();
    void Toggle_Draw_Divergence();
    void Toggle_Draw_Streamlines();
    void Toggle_Use_Streamline_Seed();
    void Set_Streamline_Seed(const unsigned int seed=0);
    void Lengthen_Streamlines();
    void Shorten_Streamlines();
    void Toggle_Draw_Vorticity();
    void Normalize_Vorticity_Color_Map();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Toggle_Velocity_Mode, "Toggle velocity mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Toggle_Velocity_Mode_And_Draw, "Toggle velocity mode and draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Toggle_Arrowhead, "Toggle arrowhead");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Toggle_Draw_Divergence, "Toggle draw divergence");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Toggle_Draw_Streamlines, "Toggle draw streamlines");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Toggle_Use_Streamline_Seed, "Toggle draw consistent streamlines");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Lengthen_Streamlines, "Lengthen streamlines");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Shorten_Streamlines, "Shorten streamlines");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Toggle_Draw_Vorticity, "Toggle draw vorticity");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D, Normalize_Vorticity_Color_Map, "Normalize vorticity map based on current frame");

private:
    void Reinitialize();
    void Update_Divergence();
    void Update_Streamlines();
    void Update_Vorticity();
};

}

#endif
