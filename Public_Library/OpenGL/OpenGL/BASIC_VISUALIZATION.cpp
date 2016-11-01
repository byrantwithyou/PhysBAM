//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <OpenGL/OpenGL/BASIC_VISUALIZATION.h>
#include <OpenGL/OpenGL/OPENGL_AXES.h>
#include <OpenGL/OpenGL/OPENGL_LIGHT.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_WINDOW.h>
#include <climits>
#include <sstream>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> BASIC_VISUALIZATION<T>::
BASIC_VISUALIZATION(STREAM_TYPE stream_type) 
    :stream_type(stream_type),opengl_axes(0),opengl_world(stream_type),set_window_position(false),opengl_window_title("OpenGL Visualization"),add_axes(true),render_offscreen(false),
    opt_left_handed(false),opt_smooth(false),selection_enabled(true),selected_object(0)
{
    reset_view_cb={[this](){Reset_View();},"Center view on selection (or reset if none)"};
    reset_up_cb={[this](){Reset_Up();},"Reset up vector for view"};
    toggle_axes_cb={[this](){Toggle_Axes();},"Toggle axes"};
    draw_all_objects_cb={[this](){Draw_All_Objects();},"draw_all_objects"};
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> BASIC_VISUALIZATION<T>::
~BASIC_VISUALIZATION() 
{
    for(int i=owned_components.m-1;i>=0;i--) delete owned_components(i);
    delete opengl_axes;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Initialize(PARSE_ARGS &parse_args)
{
    Parse_Args(parse_args);
    PreInitialize_OpenGL_World();
    Initialize_Components_And_Key_Bindings();
    Goto_Start_Frame();
    PostInitialize_OpenGL_World();
}
//#####################################################################
// Function Run
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Run()
{
    opengl_world.window->Main_Loop();
    if(render_offscreen) Render_Offscreen();
}
//#####################################################################
// Function Initialize_And_Run
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Initialize_And_Run(PARSE_ARGS &parse_args)
{
    Initialize(parse_args);
    Run();
}
//#####################################################################
// Function Add_Component
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Add_Component(OPENGL_COMPONENT<T>* component,const std::string &name,const char toggle_draw_key,const int flags)
{
    LOG::cout<<"Using Component '"<<name<<"'"<<std::endl;
    component->Set_Name(name);
    if(toggle_draw_key!='\0') opengl_world.Append_Bind_Key(toggle_draw_key,component->viewer_callbacks.Get("toggle_draw"));
    if(flags & SELECTABLE) component->selectable=true;
    if(flags & OWNED) owned_components.Append(component);
    if(flags & START_HIDDEN) component->viewer_callbacks.Get("toggle_draw").func();
    component_list.Append(component);
    component_by_name.Insert(name,component);
}
//#####################################################################
// Function Find_Component
//#####################################################################
template<class T> const OPENGL_COMPONENT<T>* BASIC_VISUALIZATION<T>::
Find_Component(const std::string &name) const
{
    if(const OPENGL_COMPONENT<T>* const* comp=component_by_name.Get_Pointer(name)) return *comp;
    return 0;
}
//#####################################################################
// Function Find_Component
//#####################################################################
template<class T> OPENGL_COMPONENT<T>* BASIC_VISUALIZATION<T>::
Find_Component(const std::string &name)
{
    if(OPENGL_COMPONENT<T>** comp=component_by_name.Get_Pointer(name)) return *comp;
    return 0;
}
//#####################################################################
// Function Add_Arguments
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS &parse_args)
{
    width=1024;
    height=768;
    fovy=20;
    parse_args.Add("-w",&width,"width","window width");
    parse_args.Add("-h",&height,"height","window height");
    parse_args.Add("-window_position",&window_position,"position","initial window corner position");
    parse_args.Add("-window_title",&opengl_window_title,"title","window title");
    parse_args.Add("-fov",&fovy,"angle","full field of view (y direction) in degrees");
    parse_args.Add("-offscreen",&render_offscreen,"render offscreen");
    parse_args.Add("-smooth",&opt_smooth,"smooth defaults");
    parse_args.Add("-camera_script",&camera_script_filename,"script","camera script filename");
    parse_args.Add("-keys",&initialization_key_sequence,"keys","initialization key sequence");
    parse_args.Add("-left_handed",&opt_left_handed,"treat coordinate system as left handed");
}
//#####################################################################
// Function Parse_Arguments
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS& parse_args)
{
    if(opt_smooth) OPENGL_PREFERENCES::Set_Smooth_Defaults();
    opengl_world.Set_Left_Handed(opt_left_handed);
}
//#####################################################################
// Function Parse_Args
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Parse_Args(PARSE_ARGS &parse_args)
{
    Add_Arguments(parse_args);
    parse_args.Parse();
    Parse_Arguments(parse_args);
}
//#####################################################################
// Function Initialize_Components_And_Key_Bindings
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Initialize_Components_And_Key_Bindings()
{
    opengl_world.Set_Key_Binding_Category("Default Keys (BASIC_VISUALIZATION<T>)");

    // Remove some silly key bindings
    opengl_world.Unbind_Keys("^m^n^o^j^k^h^l^v^d");

    opengl_world.Bind_Key("^r",reset_view_cb);
    opengl_world.Bind_Key("^u",reset_up_cb);
    opengl_world.Bind_Key("^a",toggle_axes_cb);

    if(!camera_script_filename.empty()){
        opengl_world.Bind_Key("^c",{[this](){opengl_world.Save_View(camera_script_filename,true);},"Save view"});
        opengl_world.Bind_Key('c',{[this](){opengl_world.Load_View(camera_script_filename,true);},"Load view"});}

    opengl_world.Set_Key_Binding_Category("User-Defined Keys");
}
//#####################################################################
// Function Add_OpenGL_Initialization
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Add_OpenGL_Initialization()
{
    opengl_world.Set_Ambient_Light(0.3);
    opengl_world.Add_Light(new OPENGL_LIGHT(VECTOR<double,3>(0.3,0.3,1),.3));
    opengl_world.Add_Light(new OPENGL_LIGHT(VECTOR<double,3>(-0.3,0.3,1),.3));
    opengl_world.Set_Zoom_Direction(false);
    opengl_world.Set_Lighting_For_Wireframe(false);
    opengl_world.fovy=fovy;

    if(OPENGL_PREFERENCES::smooth_points) glEnable(GL_POINT_SMOOTH);
    if(OPENGL_PREFERENCES::smooth_lines) glEnable(GL_LINE_SMOOTH);
    glLineWidth(OPENGL_PREFERENCES::line_width);
    glPointSize(OPENGL_PREFERENCES::point_size);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
}
//#####################################################################
// Function PreInitialize_OpenGL_World
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
PreInitialize_OpenGL_World()
{
    opengl_world.Initialize(opengl_window_title,width,height,render_offscreen);
    if(set_window_position) opengl_world.window->Request_Move(window_position.x,window_position.y);
}
//#####################################################################
// Function PostInitialize_OpenGL_World
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
PostInitialize_OpenGL_World()
{
    Add_OpenGL_Initialization();
    Reset_Objects_In_World();
    Initialize_Scene();
    
    if(selection_enabled){
        opengl_world.load_names_for_selection=true;
        opengl_world.process_hits_cb=[=](GLint hits,GLuint buffer[],int modifiers){Process_Hits(hits,buffer,modifiers);};}
}
//#####################################################################
// Function Initialize_Scene
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Initialize_Scene()
{
    if(!initialization_key_sequence.empty()){
        LOG::cout<<"Initialization key sequence: '"<<initialization_key_sequence<<"'"<<std::endl;
        ARRAY<OPENGL_KEY> key_list;OPENGL_KEY::Parse_Key_Sequence(initialization_key_sequence,key_list);
        for(int i=0;i<key_list.m;i++)opengl_world.Handle_Keypress_Main(key_list(i),0,0);}

    if(camera_script_filename.empty() || !opengl_world.Load_View(camera_script_filename))
        opengl_world.Center_Camera_On_Scene();
}
//#####################################################################
// Function Update_OpenGL_Strings
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Update_OpenGL_Strings()
{
    std::ostringstream output_stream;
    if(selected_object) selected_object->Print_Selection_Info(output_stream);
    opengl_world.Add_String(output_stream.str());
}
//#####################################################################
// Function Reset_Objects_In_World
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Reset_Objects_In_World()
{
    opengl_world.Clear_All_Objects();
    // Add components
    for(int i=0;i<component_list.m;i++) opengl_world.Add_Object(component_list(i),true,true);
    if(add_axes){
        if(!opengl_axes) opengl_axes=new OPENGL_AXES<T>(stream_type);
        opengl_world.Add_Object(opengl_axes,false);}
}
//#####################################################################
// Function Reset_View
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Reset_View()
{
    if(selected_object) opengl_world.Center_Camera_On_Bounding_Box(selected_object->Selection_Bounding_Box(),false);
    else{
        opengl_world.Center_Camera_On_Scene();
        opengl_world.Reset_Camera_Orientation();}
}
//#####################################################################
// Function Reset_Up
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Reset_Up()
{
    opengl_world.Reset_Camera_Orientation(true);
}
//#####################################################################
// Function Toggle_Axes
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Toggle_Axes()
{
    add_axes=!add_axes;
    Reset_Objects_In_World();
}
//#####################################################################
// Function Draw_All_Objects
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Draw_All_Objects()
{
    for(int i=0;i<component_list.m;i++) component_list(i)->Draw_All_Objects();
}
//#####################################################################
// Function Process_Hits
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Process_Hits(GLint hits,GLuint buffer[],int modifiers)
{
#ifndef NDEBUG
    opengl_world.Print_Hits(hits,buffer);
#endif

    int current_priority=INT_MIN;
    int current_idx=-1;
    float current_min_depth=FLT_MAX;
    int idx=0;
    for(int i=0;i<(int)hits;i++){
        GLint names=buffer[idx];
        int object_id=buffer[idx+3];
        ARRAY_VIEW<GLuint> name_buffer(names-1,&buffer[idx+4]);
        PHYSBAM_ASSERT(0<=object_id && object_id<opengl_world.object_list.m && opengl_world.object_list(object_id)->selectable);
        unsigned int denom=0xffffffff;
        T min_depth=(T)buffer[idx+1]/denom;
        int this_priority=opengl_world.object_list(object_id)->Get_Selection_Priority(name_buffer);
        int this_idx=idx;
        idx +=names+3;
        if(this_priority<0) continue;

        // Ties and roundoff are likely, so be careful about it.
        float depth_difference=min_depth-current_min_depth,tolerance=1e-3f;
        if(depth_difference>tolerance) continue;
        if(this_priority<current_priority && depth_difference>-tolerance) continue;
        if(this_priority==current_priority && depth_difference>=0) continue;

        current_priority=this_priority;
        current_min_depth=min_depth;
        current_idx=this_idx;}

    if(selected_object) selected_object->Clear_Selection();
    selected_object=0;

    if(current_idx>=0){
        GLint names=buffer[current_idx];
        int object_id=buffer[current_idx+3];
        PHYSBAM_ASSERT(0<=object_id && object_id<opengl_world.object_list.m && opengl_world.object_list(object_id)->selectable);
        selected_object=opengl_world.object_list(object_id);
        ARRAY_VIEW<GLuint> name_buffer(names-1,&buffer[current_idx+4]);
        if(!selected_object->Set_Selection(name_buffer,modifiers))
            selected_object=0;}

    Selection_Callback();
}
//#####################################################################
// Function Set_Current_Selection
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Set_Current_Selection(OPENGL_OBJECT<T>* object)
{
    if(selected_object){
        selected_object->Clear_Selection();
        selected_object=0;}

    if(object) selected_object=object;
}
//#####################################################################
// Function Selection_Callback
//#####################################################################
template<class T> void BASIC_VISUALIZATION<T>::
Selection_Callback()
{
    Update_OpenGL_Strings();
    glutPostRedisplay();
}
namespace PhysBAM{
template class BASIC_VISUALIZATION<double>;
template class BASIC_VISUALIZATION<float>;
}
