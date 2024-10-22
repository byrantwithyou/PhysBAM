//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANIMATED_VISUALIZATION
//#####################################################################
#ifndef __ANIMATED_VISUALIZATION__
#define __ANIMATED_VISUALIZATION__

#include <OpenGL/OpenGL/BASIC_VISUALIZATION.h>
#include <climits>

namespace PhysBAM
{

// Needed to deal with having multiple-question prompt in the event-loop-based glut framework
struct CAPTURE_FRAMES_PROMPT_STATE
{
    int step;
    std::string filename_pattern;
    int start_frame, end_frame;
    int jpeg_quality;
};

template<class T>
class ANIMATED_VISUALIZATION:public BASIC_VISUALIZATION<T>
{
public:
    using BASIC_VISUALIZATION<T>::opengl_world;using BASIC_VISUALIZATION<T>::component_list;
    using BASIC_VISUALIZATION<T>::Set_Current_Selection;using BASIC_VISUALIZATION<T>::selected_object;

    ANIMATED_VISUALIZATION(VIEWER_DIR& viewer_dir);
    virtual ~ANIMATED_VISUALIZATION(){}

private:
    void Capture_Frames(const std::string &filename_pattern, int capture_start_frame, int capture_end_frame=INT_MAX, int jpeg_quality=95, bool swap_buffers=true);

protected:
    void Add_Arguments(PARSE_ARGS &parse_args) override;
    void Parse_Arguments(PARSE_ARGS &parse_args) override;
    void Add_OpenGL_Initialization() override;
    void Initialize_Components_And_Key_Bindings() override;
    void Update_OpenGL_Strings() override;
    void Goto_Start_Frame() override;
    void Render_Offscreen() override;

    virtual void Set_Frame();
    virtual void Pre_Frame_Extra() {}
    virtual void Set_Frame_Extra() {}

public:
    void    Next_Frame();
    void    Prev_Frame();
    void    Goto_Frame();
    void    Reset();
    void    Toggle_Play();
    void    Toggle_Loop();
    void    Toggle_Fixed_Frame_Rate();
    void    Goto_Last_Frame();
    
    OPENGL_CALLBACK next_frame_cb,prev_frame_cb,goto_frame_cb,reset_cb,toggle_play_cb,toggle_loop_cb,null_cb;
    OPENGL_CALLBACK toggle_fixed_frame_rate_cb,goto_last_frame_cb,goto_frame_prompt_cb,capture_frames_cb,capture_frames_prompt_cb;
    VIEWER_DIR& viewer_dir;
    int substeps_level;
    
private:
    void    Goto_Frame_Prompt();
    void    Capture_Frames();
    void    Capture_Frames_Prompt();
    CAPTURE_FRAMES_PROMPT_STATE capture_frames_prompt_state;

protected:
    bool    play;
    bool    loop;
    bool    fixed_frame_rate;
    float   last_frame_time;
    int     start_frame, stop_frame, frame_rate;
    std::string saved_frame_filename;
    int     jpeg_quality;
    std::string frame_title;
};

}

#endif
