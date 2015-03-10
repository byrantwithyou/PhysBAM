//#####################################################################
// Copyright 2004, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_LEVELSET_2D__
#define __OPENGL_COMPONENT_LEVELSET_2D__

#include <OpenGL/OpenGL/OPENGL_LEVELSET_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T>
class OPENGL_COMPONENT_LEVELSET_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::stream_type;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::is_animation;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_LEVELSET_2D(STREAM_TYPE stream_type,const std::string& levelset_filename_input,const std::string filename_set_input="");
    virtual ~OPENGL_COMPONENT_LEVELSET_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display() const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* current_selection) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded == frame;}

    void Toggle_Color_Mode();
    void Toggle_Smooth();
    void Toggle_Normals();
    void Toggle_Draw_Mode();
    void Toggle_Draw_Sign();
    void Next_Set();
    void Previous_Set();
    void Toggle_Draw_Multiple_Levelsets();
    void Toggle_Draw_Ghost_Values();
    
private:
    void Reinitialize(const bool force_even_if_not_drawn=false);
    template<class> friend class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D;

public:
    OPENGL_LEVELSET_2D<T>* opengl_levelset;
    ARRAY<OPENGL_LEVELSET_2D<T>*> opengl_levelsets;

private:
    std::string levelset_filename;
    std::string filename_set;
    int frame_loaded;
    int set;
    bool use_sets;
    int set_loaded;
    bool valid;
    bool draw_multiple_levelsets;
};

}

#endif
