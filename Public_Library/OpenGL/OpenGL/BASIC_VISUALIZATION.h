//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIC_VISUALIZATION
//#####################################################################
#ifndef __BASIC_VISUALIZATION__
#define __BASIC_VISUALIZATION__

#include <Tools/Arrays/ARRAY.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
class PARSE_ARGS;
template<class T> class OPENGL_COMPONENT;
template<class T> class OPENGL_AXES;

template<class T>
class BASIC_VISUALIZATION
{
public:
    enum WORKAROUND{OWNED=1,SELECTABLE=2,START_HIDDEN=4};

    BASIC_VISUALIZATION();
    virtual ~BASIC_VISUALIZATION();

    void Initialize(PARSE_ARGS &parse_args);
    void Run();
    void Initialize_And_Run(PARSE_ARGS &parse_args);

    virtual void Process_Hits(GLint hits, GLuint buffer[]);

    DEFINE_CALLBACK_CREATOR(BASIC_VISUALIZATION, Draw_All_Objects);
protected:
    virtual void Add_Arguments(PARSE_ARGS &parse_args);
    virtual void Parse_Arguments(PARSE_ARGS &parse_args);
    virtual void Initialize_Components_And_Key_Bindings();
    virtual void Add_OpenGL_Initialization();
    virtual void Initialize_Scene();
    virtual void Update_OpenGL_Strings();
    virtual void Goto_Start_Frame() {}
    virtual void Render_Offscreen() {}
    virtual void Selection_Callback();

    void Add_Component(OPENGL_COMPONENT<T>* component,const std::string &name,const char toggle_draw_key,const int flags);
    const OPENGL_COMPONENT<T>* Find_Component(const std::string& name) const;
    OPENGL_COMPONENT<T>* Find_Component(const std::string& name);
    void Set_Current_Selection(OPENGL_SELECTION<T>* selection);
    int &Selection_Priority(typename OPENGL_SELECTION<T>::TYPE selection_type);

private:
    void Parse_Args(PARSE_ARGS &parse_args);
    void PreInitialize_OpenGL_World();
    void PostInitialize_OpenGL_World();
    void Reset_Objects_In_World();

private:
    void Reset_View();
    void Reset_Up();
    void Toggle_Axes();
    void Draw_All_Objects();

    DEFINE_CALLBACK_CREATOR(BASIC_VISUALIZATION,Reset_View);
    DEFINE_CALLBACK_CREATOR(BASIC_VISUALIZATION,Reset_Up);
    DEFINE_CALLBACK_CREATOR(BASIC_VISUALIZATION,Toggle_Axes);

public:
    ARRAY<OPENGL_COMPONENT<T>*> component_list;
    ARRAY<OPENGL_COMPONENT<T>*> owned_components;
    HASHTABLE<std::string,OPENGL_COMPONENT<T>*> component_by_name;
    OPENGL_AXES<T> * opengl_axes;

    OPENGL_WORLD<T> opengl_world;
    int width,height;
    bool set_window_position;
    VECTOR<int,2> window_position;
    T fovy;
    std::string opengl_window_title;
    bool add_axes;
    bool render_offscreen;
    std::string camera_script_filename;
    std::string initialization_key_sequence;
    bool opt_left_handed,opt_smooth;

    // Selection stuff
    bool selection_enabled;
    OPENGL_SELECTION<T>* current_selection;
    ARRAY<int> selection_priority;     // higher priority takes precedence; priority=0 is unselectable
};

}

#endif
