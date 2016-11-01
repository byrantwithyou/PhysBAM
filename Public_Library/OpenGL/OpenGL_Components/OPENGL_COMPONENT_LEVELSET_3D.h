//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_LEVELSET_3D__
#define __OPENGL_COMPONENT_LEVELSET_3D__

#include <Core/Arrays/ARRAY.h>
#include <OpenGL/OpenGL/OPENGL_LEVELSET_MULTIVIEW.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM{

template<class T>
class OPENGL_COMPONENT_LEVELSET_3D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::Is_Up_To_Date;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_LEVELSET_3D(STREAM_TYPE stream_type,const std::string& levelset_filename,
                                 const std::string& triangulated_surface_filename = "",
                                 const std::string& filename_set_input = "",
                                 const std::string& filename_triangulated_surface_set_input = "",
                                 bool write_generated_triangulated_surface = false,
                                 bool check_triangulated_surface_file_time = true);

    virtual ~OPENGL_COMPONENT_LEVELSET_3D();

    void Set_Surface_Material(const OPENGL_MATERIAL &front_surface_mat,
                              const OPENGL_MATERIAL &back_surface_mat);
    void Set_Overlayed_Surface_Material(const OPENGL_MATERIAL &overlayed_surface_mat);
    void Set_Slice_Color(const OPENGL_COLOR &inside_slice_color,
                         const OPENGL_COLOR &outside_slice_color);

    bool Valid_Frame(int frame_input) const override;

    void Display() const override;
    virtual RANGE<TV> Bounding_Box() const override;
    void Print_Selection_Info(std::ostream& output_stream) const override;
    void Turn_Smooth_Shading_On() override;
    void Turn_Smooth_Shading_Off() override;
    virtual void Slice_Has_Changed() override;

    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;
    bool Use_Sets() const {return use_sets;}

    void Toggle_Display_Overlay();
    void Toggle_Slice_Color_Mode();
    void Toggle_Smooth_Slice();
    void Next_Set();
    void Previous_Set();
    void Toggle_Draw_Multiple_Levelsets();

private:
    void Reinitialize();
    void Reinitialize_Levelset(const std::string& levelset_filename, const std::string& triangulated_surface_filename, OPENGL_LEVELSET_MULTIVIEW<T>* levelset_multiview);

public:
    OPENGL_LEVELSET_MULTIVIEW<T>* opengl_levelset_multiview;
    ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>* > opengl_levelset_multiviews;

private:
    std::string levelset_filename;
    std::string triangulated_surface_filename;
    std::string filename_set;
    std::string filename_triangulated_surface_set;
    bool write_generated_triangulated_surface;
    int frame_loaded;
    bool check_triangulated_surface_file_time;
    int set;
    int set_loaded;
    bool use_sets;
    bool draw_multiple_levelsets;
public:
    int ghost_cells;
    TV_INT selected_cell;
    TV_INT selected_node;
};

}

#endif
