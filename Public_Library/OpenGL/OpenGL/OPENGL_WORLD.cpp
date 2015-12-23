//#####################################################################
// Copyright 2002-2008, Eran Guendelman, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_WORLD
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/constants.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Utilities/TIMER.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <OpenGL/OpenGL/OPENGL_ARCBALL.h>
#include <OpenGL/OpenGL/OPENGL_LIGHT.h>
#include <OpenGL/OpenGL/OPENGL_MOUSE_HANDLER.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_WINDOW_GLUT.h>
#include <OpenGL/OpenGL/OPENGL_WINDOW_PBUFFER.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
namespace PhysBAM{

using ::std::log;
using ::std::pow;

template<class T> inline OPENGL_WORLD<T>*& Opengl_World()
{
    static OPENGL_WORLD<T>* opengl_world=0;
    return opengl_world;
}

//#####################################################################
// Constructor OPENGL_WORLD
//#####################################################################
template<class T> OPENGL_WORLD<T>::
OPENGL_WORLD(STREAM_TYPE stream_type)
    :initialized(false),stream_type(stream_type),smooth_shading(false),ambient_light(OPENGL_COLOR::White()),fovy(50),mode_2d(false),load_names_for_selection(false),
    window(0),fill_mode(DRAW_FILLED),enable_lighting_for_wireframe(false),white_background(false),
    display_strings(true),show_object_names(false),display_object_names_in_corner(false),view_auto_help(false),
    timer_id(0),idle_delay(0),idle_timer(0),view_target_timer(0),frame_counter_timer(0),frames_rendered(0),frames_per_second(0),show_frames_per_second(true),
    left_handed_coordinate_system(false),nearclip_factor(.0625),farclip_factor(4),nearclip(nearclip_factor),farclip(farclip_factor),
    arcball(new OPENGL_ARCBALL<T>(stream_type,*this)),camera_distance(1),arcball_matrix(arcball->Value()),rotation_matrix(arcball->Value()),
    zoom_direction(1),translation_direction(1),oldmousex(0),oldmousey(0),do_mouse_rotation(false),do_mouse_zoom(false),do_mouse_target_xy(false),
    do_mouse_target_z(false),external_mouse_handler(0),shift_was_pressed(false),ctrl_was_pressed(false),
    current_key_binding_category("User-Defined Keys"),current_key_binding_category_priority(1),prompt_mode(false),
    process_hits_cb(0),selection_mode(false),current_selection(0)
{
    if(Opengl_World<T>()!=0) PHYSBAM_FATAL_ERROR(); 
    Opengl_World<T>()=this;

    key_bindings.Resize(0,OPENGL_KEY::MAX_KEY_INDEX,0,OPENGL_KEY::MAX_MODIFIER_INDEX);

    Set_Key_Binding_Category("Default Keys (OPENGL_WORLD)");
    Set_Key_Binding_Category_Priority(1000);
    Bind_Key("^q",{[this](){exit(0);},"Quit"});
    Bind_Key("^w",{[this](){fill_mode=(fill_mode+1)%3;},"Toggle wireframe mode"});

    // TODO: fix full screen
    Bind_Key("^m",{[this](){camera_distance*=(T).75;nearclip*=(T).75;farclip*=(T).75;},"Zoom by 3/4"});
    Bind_Key("^n",{[this](){camera_distance/=(T).75;nearclip/=(T).75;farclip/=(T).75;},"Zoom by 4/3"});

    Bind_Key("^d",{[this]()
            {
                LOG::cout<<"Filename for saving image in .ppm or .jpg (if supported) formats (return for default, - for cancel)\n";
                std::string filename;
                getline(std::cin,filename);
                size_t len=filename.length();
                if(filename[len-1]=='\n'){--len;filename=filename.substr(0,len);}
                if(len<=0) return;
                LOG::cout<<"saving to %s..."<<filename;
                Save_Screen(filename,true);
                LOG::cout<<"\n";
            },"Save screen"});
    Bind_Key("^|",{[this]()
            {
                show_frames_per_second=!show_frames_per_second;
                if(show_frames_per_second){frame_counter_timer=1;Prepare_For_Idle();}
                else frame_counter_timer=0;
            },"Toggle Show Frames/Second"});
    Bind_Key("^g",{[this]()
            {
                if(smooth_shading){
                    smooth_shading=false;
                    for(int i=0;i<object_list.m;i++)
                        if(can_toggle_smooth_shading(i))
                            object_list(i)->Turn_Smooth_Shading_Off();}
                else{
                    smooth_shading=true;
                    for(int i=0;i<object_list.m;i++)
                        if(can_toggle_smooth_shading(i))
                            object_list(i)->Turn_Smooth_Shading_On();}
            },"Toggle smooth shading"});
    Bind_Key("^a",{[this](){show_object_names=!show_object_names;},"Toggle object names"});
    Bind_Key("^t",{[this](){white_background=!white_background;},"Toggle Background"});
    Bind_Key("|",{[this](){window->Request_Resize(640,480);},"Resize Window to 640x480"});
    Bind_Key("?",{[this](){view_auto_help=!view_auto_help;},"Display Help"});
    Bind_Key("<F1>",{[this](){display_strings=!display_strings;},"Toggle display strings"});

    Set_Key_Binding_Category("User-Defined Keys");
    Set_Key_Binding_Category_Priority(1);
    Set_View_Target_Timer(1);
    arcball->world=this;
}
//#####################################################################
// ~OPENGL_WORLD
//#####################################################################
template<class T> OPENGL_WORLD<T>::
~OPENGL_WORLD()
{
    Clear_All_Lights();
    delete window;
    delete arcball;
    Opengl_World<T>()=0;
}
//#####################################################################
// Singleton
//#####################################################################
template<class T> OPENGL_WORLD<T>* OPENGL_WORLD<T>::
Singleton()
{
    PHYSBAM_ASSERT(Opengl_World<T>());
    return Opengl_World<T>();
}
//#####################################################################
// Function Run_Visualization
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Run_Visualization(const std::string& window_title)
{
    Initialize(window_title);
    window->Main_Loop();
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Initialize(const std::string& window_title,const int width,const int height,const bool offscreen)
{
    if(window) delete window;
    if(offscreen)
        window=new OPENGL_WINDOW_PBUFFER<T>(*this,window_title,width,height);
    else
        window=new OPENGL_WINDOW_GLUT<T>(*this,window_title,width,height);
    initialized=true;
    Prepare_For_Idle();
    Set_View_Target_Timer(0);
    Initialize_Glut_Independent();
}
//#####################################################################
// Function Initialize_Glut_Independent
//#####################################################################
// Should call this after initializing a GL context (using glut or
// otherwise)
template<class T> void OPENGL_WORLD<T>::
Initialize_Glut_Independent()
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);

    // These are important to get screen captures to work correctly
    glPixelStorei(GL_PACK_ALIGNMENT,1);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);

    Center_Camera_On_Scene();
}
//#####################################################################
// Function Clear_All_Objects
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Clear_All_Objects()
{
    object_list.Remove_All();
}
//#####################################################################
// Function Add_Object
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Add_Object(OPENGL_OBJECT<T> *object,bool include_bounding_box,bool toggle_smooth_shading)
{
    object_list.Append(object);
    use_bounding_box.Append(include_bounding_box);
    can_toggle_smooth_shading.Append(toggle_smooth_shading);
}
//#####################################################################
// Function Clear_All_Lights
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Clear_All_Lights()
{
    lights.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Add_Light
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Add_Light(OPENGL_LIGHT *light)
{
    lights.Append(light);
}
//#####################################################################
// Function Set_Ambient_Light
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Ambient_Light(T value)
{
    ambient_light=OPENGL_COLOR::Gray(value);
}
//#####################################################################
// Function Set_Ambient_Light
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Ambient_Light(const OPENGL_COLOR& color)
{
    ambient_light=color;
}
//#####################################################################
// Function Set_Key_Binding_Category
//####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Key_Binding_Category(const std::string &category)
{
    current_key_binding_category=category;
}
//#####################################################################
// Function Set_Key_Binding_Category_Priority
//####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Key_Binding_Category_Priority(int priority)
{
    current_key_binding_category_priority=priority;
}
//#####################################################################
// Function Bind_Key
//####################################################################
template<class T> void OPENGL_WORLD<T>::
Bind_Key(const OPENGL_KEY& key,OPENGL_CALLBACK callback)
{
    Unbind_Key(key);
    Append_Bind_Key(key,callback);
}
//#####################################################################
// Function Bind_Key
//####################################################################
template<class T> void OPENGL_WORLD<T>::
Bind_Key(const std::string& key,OPENGL_CALLBACK callback)
{
    Bind_Key(OPENGL_KEY::From_String(key),callback);
}
//#####################################################################
// Function Append_Bind_Key
//####################################################################
template<class T> void OPENGL_WORLD<T>::
Append_Bind_Key(const OPENGL_KEY& key,OPENGL_CALLBACK callback)
{
    key_bindings(key.Index()).Append(callback);
    int index;
    for(index=0;index<key_bindings_by_category.m;index++) if(key_bindings_by_category(index).name==current_key_binding_category) break;
    if(index>=key_bindings_by_category.m){ // Insert in location dependent on priority
        for(index=0;index<key_bindings_by_category.m && key_bindings_by_category(index).priority<=current_key_binding_category_priority;index++){}
        key_bindings_by_category.Insert(OPENGL_KEY_BINDING_CATEGORY(current_key_binding_category,current_key_binding_category_priority),index);}
    key_bindings_by_category(index).key_bindings.Append(PAIR<OPENGL_KEY,OPENGL_CALLBACK>(key,callback));
}
//#####################################################################
// Function Append_Bind_Key
//####################################################################
template<class T> void OPENGL_WORLD<T>::
Append_Bind_Key(const std::string& key,OPENGL_CALLBACK callback)
{
    Append_Bind_Key(OPENGL_KEY::From_String(key),callback);
}
//#####################################################################
// Function Unbind_Key
//####################################################################
template<class T> void OPENGL_WORLD<T>::
Unbind_Key(const OPENGL_KEY& key)
{
    key_bindings(key.Index()).Remove_All();
    for(int i=0;i<key_bindings_by_category.m;i++)
        for(int j=0;j<key_bindings_by_category(i).key_bindings.m;j++)
            if(key_bindings_by_category(i).key_bindings(j).x==key){
                key_bindings_by_category(i).key_bindings(j).y={0,0};
                key_bindings_by_category(i).key_bindings.Remove_Index(j);}
}
//#####################################################################
// Function Unbind_Keys
//####################################################################
template<class T> void OPENGL_WORLD<T>::
Unbind_Keys(const std::string& keys)
{
    ARRAY<OPENGL_KEY> key_list;OPENGL_KEY::Parse_Key_Sequence(keys,key_list);
    for(int i=0;i<key_list.m;i++)Unbind_Key(key_list(i));
}
//#####################################################################
// Function Set_Idle_Callback
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Idle_Callback(OPENGL_CALLBACK callback,const T delay)
{
    bool need_prepare=!idle_timer || idle_timer>delay;
    idle_callback=callback;
    idle_delay=delay;
    idle_timer=idle_delay;
    if(need_prepare) Prepare_For_Idle();
}
//#####################################################################
// Function Set_View_Target_Timer
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_View_Target_Timer(const T view_target_timer_input)
{
    bool need_prepare=!view_target_timer || view_target_timer>view_target_timer_input;
    view_target_timer=view_target_timer_input;
    if(need_prepare) Prepare_For_Idle();
}
//#####################################################################
// Function Prepare_For_Idle
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Prepare_For_Idle()
{
    if(!timer_id) timer_id=TIMER::Singleton()->Register_Timer();
    if(!initialized) return;
    bool use_idle=idle_callback.func && !idle_delay;
    window->Setup_Idle(use_idle);

    // If not use idle then wait awhile without consuming CPU
    if(!use_idle && (idle_timer || view_target_timer || frame_counter_timer)){
        T wait=min(idle_timer?idle_timer:FLT_MAX,view_target_timer?view_target_timer:FLT_MAX,frame_counter_timer?frame_counter_timer:FLT_MAX);
        window->Setup_Timer(wait);}
}
//#####################################################################
// Function Handle_Idle
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Handle_Idle()
{
    double delta_time=TIMER::Singleton()->Peek_And_Reset_Time(timer_id)/1000;
    if(delta_time<0) delta_time=0; // this can happen due to roll-over
    bool need_redisplay=false;
    if(view_target_timer > 0){
        view_target_timer-=delta_time;
        if(view_target_timer<=0){
            need_redisplay=true;
            view_target_timer=0;}}
    if(show_frames_per_second){
        frame_counter_timer-=delta_time;
        if(frame_counter_timer<0){frame_counter_timer=1;frames_per_second=frames_rendered;frames_rendered=0;}}
    if(idle_callback.func){
        idle_timer-=delta_time;
        if(idle_timer<=0){
            need_redisplay=true;
            idle_callback.func();
            idle_timer=idle_delay;}}
    if(need_redisplay) window->Redisplay();
    Prepare_For_Idle();
}
//#####################################################################
// Function Handle_Timer
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Handle_Timer()
{
    if(idle_callback.func || idle_delay) Handle_Idle();
}
//#####################################################################
// Function Handle_Display_Prompt_Only
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Handle_Display_Prompt_Only()
{
    glPushAttrib(GL_ENABLE_BIT);
    glEnable(GL_SCISSOR_TEST);
    glScissor(0,window->Height()-30,window->Width(),30);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    Display_Prompt_Strings();
    glFlush();
    glDisable(GL_SCISSOR_TEST); // workaround for bug in 2.1.2 NVIDIA 177.13
    glPopAttrib();
}
//#####################################################################
// Function Render_World()
//#####################################################################
template<class T> void OPENGL_WORLD<T>::Render_World(bool selecting,bool swap_buffers)
{
    if(prompt_mode){Handle_Display_Prompt_Only();return;}

    glMatrixMode(GL_PROJECTION); // probably doesn't need to be reset each time!
    if(!selecting) glLoadIdentity();
    gluPerspective(fovy,window->Width()/(T)window->Height(),nearclip,farclip);
    if(left_handed_coordinate_system){ // convert to a left-handed coordinate system
        glScalef(-1,1,1);
        glFrontFace(GL_CW);}
    else
        glFrontFace(GL_CCW);

    // Set the background
    OPENGL_COLOR wireframe_color;
    if(white_background){glClearColor(1,1,1,0);wireframe_color=OPENGL_COLOR::Black();}
    else{glClearColor(0,0,0,0);wireframe_color=OPENGL_COLOR::White();}

    // Clear z-buffer and color buffer
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode (GL_MODELVIEW); // establish camera coordinate frame in scene
    glLoadIdentity();
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambient_light.rgba);
    glEnable(GL_LIGHTING);
    for(int i=0;i<lights.m;i++) lights(i)->Send_To_GL_Pipeline(i);

    glTranslatef(0,0,-camera_distance);
    glMultMatrix(arcball_matrix.x); // assumes MATRIX_4X4::x is a column-major array
    glMultMatrix(rotation_matrix.x); // (as OpenGL expects)
    OpenGL_Translate(-target_position);

    Update_Clipping_Planes();

    if(fill_mode==DRAW_FILLED)
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    else if(fill_mode==DRAW_WIREFRAME){
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        if(!enable_lighting_for_wireframe){
            glDisable(GL_LIGHTING);wireframe_color.Send_To_GL_Pipeline();}}
    else if(!selecting){
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0,1.0);
        for(int i=0;i<object_list.m;i++) if((!selecting && object_list(i)->visible) || object_list(i)->selectable) object_list(i)->Display();
        glDisable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        if(!enable_lighting_for_wireframe){
            glDisable(GL_LIGHTING);wireframe_color.Send_To_GL_Pipeline();}}

    for(int i=0;i<object_list.m;i++)
        if((!selecting && object_list(i)->visible) || object_list(i)->selectable){
            if(load_names_for_selection) glLoadName(i);
            object_list(i)->Display();}
    if(view_target_timer>0) Opengl_World<T>()->Display_Target();
    glEnable(GL_LIGHTING);

    if(!selecting){
        if(view_auto_help){
            Display_Auto_Help();}
        else{
            if(prompt_mode) Display_Prompt_Strings();
            else Display_Strings();

            if(show_object_names){
                if(display_object_names_in_corner) Display_Object_Names_In_Corner();
                else Display_Object_Names();}}}
    ++frames_rendered;

    glDisable(GL_LIGHTING);
    GLenum gl_error=glGetError();
    if(gl_error !=GL_NO_ERROR) LOG::cout<<"OpenGL Error: "<<gluErrorString(gl_error)<<std::endl;
}
//#####################################################################
// Function Update_Clipping_Planes
//#####################################################################
template<class T> void OPENGL_WORLD<T>::Set_External_Mouse_Handler(OPENGL_MOUSE_HANDLER* mouse_handler)
{
    external_mouse_handler=mouse_handler;
}
//#####################################################################
// Function Update_Clipping_Planes
//#####################################################################
template<class T> void OPENGL_WORLD<T>::Update_Clipping_Planes()
{
    for(int i=0;i<clipping_planes.m;i++)
        if(clipping_planes(i)) OpenGL_Clip_Plane(GL_CLIP_PLANE0+i,*clipping_planes(i));
}
//#####################################################################
// Function Add_Clipping_Plane
//#####################################################################
template<class T> GLenum OPENGL_WORLD<T>::Add_Clipping_Plane(const PLANE<T> &plane)
{
    int index=-1;
    for(int i=0;i<clipping_planes.m;i++)
        if(!clipping_planes(i)){index=i;break;}
    if(index==-1)
        index=clipping_planes.Append(new PLANE<T>(plane));
    else clipping_planes(index)=new PLANE<T>(plane);
    return GL_CLIP_PLANE0+index;
}
//#####################################################################
// Function Set_Clipping_Plane
//#####################################################################
template<class T> void OPENGL_WORLD<T>::Set_Clipping_Plane(GLenum id,const PLANE<T> &plane)
{
    int index=id-GL_CLIP_PLANE0;
    PHYSBAM_ASSERT(clipping_planes(index));
    *clipping_planes(index)=plane;
}
//#####################################################################
// Function Remove_Clipping_Plane
//#####################################################################
template<class T> void OPENGL_WORLD<T>::Remove_Clipping_Plane(GLenum id)
{
    int index=id-GL_CLIP_PLANE0;
    PHYSBAM_ASSERT(clipping_planes(index));
    delete clipping_planes(index);clipping_planes(index)=0;
}
//#####################################################################
// Function Remove_All_Clipping_Planes
//#####################################################################
template<class T> void OPENGL_WORLD<T>::Remove_All_Clipping_Planes()
{
    for(int i=0;i<clipping_planes.m;i++) delete clipping_planes(i);
    clipping_planes.Remove_All();
}
//#####################################################################
// Function Print_Hits
//#####################################################################
// This version of the function is useful for debugging purposes
// to show you the contents of the hits buffer
template<class T> void OPENGL_WORLD<T>::Print_Hits(GLint hits,GLuint buffer[])
{
    int idx=0;
    for(int i=0;i<(int)hits;i++){
        GLint names=buffer[idx];
        unsigned int denom=0xffffffff;
        double mindepth=(double)buffer[idx+1]/denom;
        double maxdepth=(double)buffer[idx+2]/denom;
        LOG::cout<<"Hit "<<i<<" [depth range="<<mindepth<<","<<maxdepth<<"]: ";
        idx +=3;
        for(int j=0;j<(int)names;j++){
            LOG::cout<<buffer[idx]<<" ";
            idx++;
        }
        LOG::cout<<std::endl;
    }
}
//#####################################################################
// Function Get_Number_Of_Valid_Hits
//#####################################################################
template<class T> int OPENGL_WORLD<T>::Get_Number_Of_Valid_Hits(GLint hits,GLuint buffer[],int buff_size)
{
    if(hits>=0) return hits; // nothing to do
    int i=0,idx=0;
    while(idx<buff_size){
        GLint names=buffer[idx];idx+=3+names;
        if(idx>=buff_size) break;
        i++;}
    return i;
}
//#####################################################################
// Function Get_Selections
//#####################################################################
template<class T> void OPENGL_WORLD<T>::Get_Selections(ARRAY<OPENGL_SELECTION<T>* > &selections,GLint hits,GLuint buffer[])
{
    selections.Remove_All();
    int idx=0;
    for(int i=0;i<(int)hits;i++){
        GLint names=buffer[idx];
        int object_id=buffer[idx+3];
        PHYSBAM_ASSERT(0<=object_id && object_id<object_list.m && object_list(object_id)->selectable);
        OPENGL_SELECTION<T>* selection=object_list(object_id)->Get_Selection(&buffer[idx+4],names-1);
        if(selection){
            unsigned int denom=0xffffffff;
            selection->min_depth=(T)buffer[idx+1]/denom;
            selection->max_depth=(T)buffer[idx+2]/denom;
            selection->hide=shift_was_pressed && ctrl_was_pressed;
            selections.Append(selection);
        }
        idx +=names+3;
    }
}
//#####################################################################
// Function Handle_Reshape_Main
//#####################################################################
// glut independent
template<class T> void OPENGL_WORLD<T>::
Handle_Reshape_Main()
{
    if(window){
        glViewport(0,0,(GLsizei)window->Width(),(GLsizei)window->Height());
        window->Redisplay();}
}

//#####################################################################
// Function Handle_Keypress_Main
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Handle_Keypress_Main(const OPENGL_KEY& key,int x,int y)
{
    VECTOR<int,2> index=key.Index();

    if(key_bindings(index).m >=1)
    {
        for(int i=0;i<key_bindings(index).m;i++)
            key_bindings(index)(i).func();
        window->Redisplay();
    }
}
//#####################################################################
// Handle_Keypress_Prompt
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Handle_Keypress_Prompt(unsigned char key)
{
    if(isprint(key)) prompt_response.push_back(key);
    else if(key==8 && !prompt_response.empty()) prompt_response.resize(prompt_response.size()-1); // BACKSPACE
    else if(key==13 || key==27) { // ENTER or ESC
        // If ESC pressed we set prompt_response to null to indicate aborted prompt
        if(key==27){prompt_response.clear();prompt_response_success=false;}
        glPopAttrib();
        prompt_mode=false;
        prompt_response_cb.func();}
    window->Redisplay();
}
//#####################################################################
// Function Handle_Click
//#####################################################################
// glut independent, although button and state should be GLUT defines
template<class T> void OPENGL_WORLD<T>::
Handle_Click_Main(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed)
{
    if(external_mouse_handler){
        external_mouse_handler->Handle_Click(button,state,x,y,ctrl_pressed,shift_pressed);
        window->Redisplay();}
    
    switch(button){
        case GLUT_LEFT_BUTTON:
            // Add selection stuff
            if((shift_pressed && state==GLUT_DOWN) || selection_mode){
                const int buff_size=512;
                GLuint selectBuf[buff_size];
                GLint viewport[4];
                shift_was_pressed=shift_pressed!=0;
                ctrl_was_pressed=ctrl_pressed!=0;
                glGetIntegerv(GL_VIEWPORT,viewport);
                glSelectBuffer(buff_size,selectBuf);
                (void) glRenderMode (GL_SELECT);

                glInitNames();
                glPushName(0);
                glMatrixMode(GL_PROJECTION);
                glPushMatrix();
                glLoadIdentity();
                gluPickMatrix((GLfloat) x,(GLfloat) (viewport[3]-y),
                              5.0,5.0,viewport);
                // glOrtho(0.0,10.0,0.0,10.0,0.0,10.0);
                Render_World(true,false);
                glMatrixMode(GL_PROJECTION);
                glPopMatrix();
                glFlush();

                GLint hits=glRenderMode (GL_RENDER);
                if(hits<0) hits=Get_Number_Of_Valid_Hits(hits,selectBuf,buff_size); // if had buffer overflow, get number of valid hit records
                if(process_hits_cb) process_hits_cb(hits,selectBuf);
            } else
            {
                VECTOR<T,2> mouseVector=Convert_Mouse_Coordinates(x,y);
                if(state==GLUT_UP){
                    shift_was_pressed=false;
                    if(do_mouse_rotation){
                        arcball->End_Drag(mouseVector); // indictes that dragging should end
                        arcball_matrix=arcball->Value(); // extracts the current matrix transform
                        // rotation stored in rotation_matrix.
                        rotation_matrix=arcball_matrix*rotation_matrix;
                        arcball_matrix=MATRIX<T,4>::Identity_Matrix();
                    }
                    do_mouse_rotation=do_mouse_target_xy=false;
                } else if(state==GLUT_DOWN){
                    if(mode_2d || ctrl_pressed) Prepare_For_Target_XY_Drag();
                    else{
                        // Update
                        arcball->Begin_Drag(mouseVector); // indicates that dragging should begin
                        do_mouse_rotation=true;
                        Set_View_Target_Timer(1);
                    }
                }
            }
            break;
        case GLUT_RIGHT_BUTTON:
            if(state==GLUT_UP){do_mouse_zoom=do_mouse_target_z=false;}
            else if(state==GLUT_DOWN){
                if(ctrl_pressed){do_mouse_target_z=true;Set_View_Target_Timer(1);}
                else do_mouse_zoom=true;}
            break;
        case GLUT_MIDDLE_BUTTON:
            if(state==GLUT_UP){do_mouse_target_xy=false;}
            else if(state==GLUT_DOWN)  Prepare_For_Target_XY_Drag();
            break;}
    oldmousex=x;
    oldmousey=y;
}
//#####################################################################
// Function Handle_Drag_Main
//#####################################################################
// glut independent
template<class T> void OPENGL_WORLD<T>::
Handle_Drag_Main(int x,int y)
{
    if(external_mouse_handler){external_mouse_handler->Handle_Drag(x,y);return;}
    
    GLint viewport[4];
    GLdouble mvmatrix[16]={   1,0,0,0,
                                0,1,0,0,
                                0,0,1,0,
                                0,0,0,1 },projmatrix[16];
    GLint realy; //  OpenGL y coordinate position
    GLdouble wx,wy,wz; //  first point returned world x,y,z coords
    GLdouble win_x,win_y,win_z;
    if(selection_mode && current_selection!=0){
        // User Selects Guidance Points that are part of the system.
        // Find the unprojected points
        glGetIntegerv (GL_VIEWPORT,viewport);
        glGetDoublev (GL_MODELVIEW_MATRIX,mvmatrix);
        glGetDoublev (GL_PROJECTION_MATRIX,projmatrix);

        //  note viewport[3] is height of window in pixels
        realy=viewport[3]-(GLint) y-1;
        gluProject((GLdouble) current_selection->x,(GLdouble) current_selection->y,(GLdouble) current_selection->z,
            mvmatrix,projmatrix,viewport,&win_x,&win_y,&win_z);
        gluUnProject((GLdouble) x,(GLdouble) realy,win_z,
            mvmatrix,projmatrix,viewport,&wx,&wy,&wz);
        
        current_selection->x=wx;
        current_selection->y=wy;
        current_selection->z=wz;}
    else{
        if(do_mouse_rotation){
            T length=(T (x-oldmousex)*(x-oldmousex)+(y-oldmousey)*(y-oldmousey))/abs(window->Width()*window->Height());
            if(length > 0.0000001){
                VECTOR<T,2> v=Convert_Mouse_Coordinates(x,y);
                arcball->Update(v); //Alters the internal state of the arcball
                arcball_matrix=arcball->Value(); //reads the matrix from the arcball
            }
            Set_View_Target_Timer(.5f);
        }
        if(do_mouse_zoom){
            T factor=pow(1.01,zoom_direction*(y-oldmousey));
            camera_distance*=factor;
            nearclip=nearclip_factor*camera_distance;
            farclip=farclip_factor*camera_distance;}
        if(do_mouse_target_xy){
            T dx=oldmousex-x;
            T dy=y-oldmousey;
            target_position+=dx*target_x_drag_vector+dy*target_y_drag_vector;
            Set_View_Target_Timer(.5f);
        }
        if(do_mouse_target_z){
            TV view_forward,view_up,view_right;
            Get_View_Frame(view_forward,view_up,view_right);

            T dy=y-oldmousey;
            target_position +=dy*view_forward*0.001*camera_distance;
            Set_View_Target_Timer(.5f);
        }
    }
    oldmousex=x;
    oldmousey=y;

    window->Redisplay();
}
//#####################################################################
// Function Prepare_For_Target_XY_Drag
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Prepare_For_Target_XY_Drag()
{
    do_mouse_target_xy=true;Set_View_Target_Timer(1);

    // Find the unprojected points
    GLint viewport[4];
    GLdouble mvmatrix[16]={ 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1 },projmatrix[16];
    GLdouble wx,wy,wz;
    GLdouble win_x,win_y,win_z;

    glGetIntegerv (GL_VIEWPORT,viewport);
    glGetDoublev (GL_MODELVIEW_MATRIX,mvmatrix);
    glGetDoublev (GL_PROJECTION_MATRIX,projmatrix);
    gluProject((GLdouble) target_position.x,(GLdouble) target_position.y,(GLdouble) target_position.z,
        mvmatrix,projmatrix,viewport,&win_x,&win_y,&win_z);
    gluUnProject((GLdouble) win_x+1,(GLdouble) win_y,win_z,
        mvmatrix,projmatrix,viewport,&wx,&wy,&wz);
    target_x_drag_vector=TV(wx,wy,wz)-target_position;
    gluUnProject((GLdouble) win_x,(GLdouble) win_y+1,win_z,
        mvmatrix,projmatrix,viewport,&wx,&wy,&wz);
    target_y_drag_vector=TV(wx,wy,wz)-target_position;
}
//#####################################################################
// Function Display_Target
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Display_Target()
{
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glTranslatef(target_position.x,target_position.y,target_position.z);
    GLfloat clear_color[4];
    glGetFloatv(GL_COLOR_CLEAR_VALUE,clear_color);
    glColor3f(1-clear_color[0],1-clear_color[1],1-clear_color[2]);
    T logcd=log(camera_distance)/log(20.);
    T smallsize=pow(20.0,floor(logcd-.5));
    glutWireCube(smallsize);
    glutWireCube(20*smallsize);
    glPopMatrix();
    glPopAttrib();
}
//#####################################################################
// Function Display_Auto_Help
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Display_Auto_Help()
{
    ARRAY<std::string> strings1,strings2;
    for(int i=0;i<key_bindings_by_category.m;i++){
        if(key_bindings_by_category(i).key_bindings.m>0){
            strings1.Append(key_bindings_by_category(i).name+":");
            strings2.Append("");
            for(int j=0;j<key_bindings_by_category(i).key_bindings.m;j++){
                std::ostringstream string_stream1;
                string_stream1.flags(std::ios::right);
                string_stream1.width(5);
                string_stream1<<key_bindings_by_category(i).key_bindings(j).x.Name();
                string_stream1<<":";
                strings1.Append(string_stream1.str());
                strings2.Append(std::string("      ")+key_bindings_by_category(i).key_bindings(j).y.help);}}}

    Display_Strings(strings2,OPENGL_COLOR::Yellow(),true,0,13,GLUT_BITMAP_8_BY_13);
    Display_Strings(strings1,OPENGL_COLOR::White(),false,0,13,GLUT_BITMAP_8_BY_13);
}
//#####################################################################
// Function Display_Strings
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Draw_Transparent_Text_Box(const ARRAY<std::string> &strings,const VECTOR<int,2> &top_left_corner,int vspace,void *font,const OPENGL_COLOR &color)
{
    int max_string_length=0;for(int i=0;i<strings.m;i++) max_string_length=max(max_string_length,glutBitmapLength(font,(const unsigned char *)strings(i).c_str()));
    int num_lines=strings.m;
    OPENGL_SHAPES::Draw_Translucent_Stripe(top_left_corner.x,top_left_corner.y,max_string_length+vspace,-(num_lines+1)*vspace,color);
}
//#####################################################################
// Function Display_Strings
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Display_Strings(const ARRAY<std::string> &strings,const OPENGL_COLOR &color,bool draw_transparent_box,int horizontal_offset,int vspace,void *font)
{
    if(!strings.m || !display_strings) return;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    int g_Height=window->Height();
    gluOrtho2D(0,window->Width(),0,g_Height);

    if(draw_transparent_box) Draw_Transparent_Text_Box(strings,VECTOR<int,2>(horizontal_offset,window->Height()),vspace,font,OPENGL_COLOR::Gray(0,0.5));

    g_Height-=vspace;
    color.Send_To_GL_Pipeline();
    for(int i=0;i<strings.m;i++){
        OpenGL_String(VECTOR<T,2>(horizontal_offset,g_Height),strings(i),font);
        g_Height-=vspace;
    }

    // set openGL states back to the way they were
    glPopAttrib();

    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}
//#####################################################################
// Function Display_Strings
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Display_Strings(bool draw_transparent_box)
{
    if(!display_strings) return;
    static OPENGL_COLOR text_color=OPENGL_COLOR::White();
    if(show_frames_per_second){
        ARRAY<std::string> strings_with_fps;
        strings_with_fps.Append(LOG::sprintf("%dfps",frames_per_second));
        strings_with_fps.Append_Elements(strings_to_print);
        Display_Strings(strings_with_fps,text_color,draw_transparent_box);}
    else
        Display_Strings(strings_to_print,text_color,draw_transparent_box);
}
//#####################################################################
// Function Display_Object_Names
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Display_Object_Names()
{
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glColor3f(1,1,1);

    for(int i=0;i<object_list.m;i++)
        if(object_list(i)->name.length() && object_list(i)->show_name)
            OpenGL_String(object_list(i)->frame->t,object_list(i)->name,GLUT_BITMAP_HELVETICA_12);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
}
//#####################################################################
// Function Display_Object_Names
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Display_Object_Names_In_Corner()
{
    Display_Object_Names();

    ARRAY<std::string> strings;

    if(show_frames_per_second) strings.Append(LOG::sprintf("%dfps",frames_per_second));

    Display_Strings(strings,OPENGL_COLOR::White(),true,0,12,GLUT_BITMAP_HELVETICA_12);
}
//#####################################################################
// Function Add_String
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Add_String(const std::string& s)
{
    // Separate s into separate lines (in case s contains carriage returns)
    for(std::string::size_type start=0;start<s.length();){
        std::string::size_type end=s.find('\n',start);
        strings_to_print.Append(s.substr(start,end-start));
        if(end==std::string::npos) break;start=end+1;}
}
//#####################################################################
// Function Add_String
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Clear_Strings()
{
    strings_to_print.Remove_All();
}
//#####################################################################
// Function Set_Zoom_Direction
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Zoom_Direction(bool up_is_zoom_in)
{
    zoom_direction=(up_is_zoom_in) ? 1 : -1;
}
//#####################################################################
// Function Set_Translation_Direction
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Translation_Direction(bool up_is_move_in)
{
    translation_direction=(up_is_move_in) ? 1 : -1;
}
//#####################################################################
// Function Set_Lighting_For_Wireframe
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Lighting_For_Wireframe(bool enable_flag)
{
    enable_lighting_for_wireframe=enable_flag;
}
//#####################################################################
// Function Set_2D_Mode
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_2D_Mode(bool mode)
{
    mode_2d=mode;
}
//#####################################################################
// Function Set_Process_Hits_Callback
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Process_Hits_Callback(PROCESS_HITS_CB process_hits_cb_input)
{
    process_hits_cb=process_hits_cb_input;
}
//#####################################################################
// Function Get_View_Frame
//   Returned vectors are already normalized (or at least should be!).
//   view_forward is a vector in world space which in the current view
//   appears to point into screen.
//   view_up appears to point up, and view_right points right.
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Get_View_Frame(TV &view_forward,TV &view_up,TV &view_right)
{
    MATRIX<T,4> matrix=arcball_matrix*rotation_matrix;

    // We want to find the vector v such that matrix * v = (0,0,-1)
    // (because OpenGL has negative z pointing into screen)
    // Since matrix^-1 = matrix^T, this is simply the negative of
    // the third row of matrix.
    view_forward=TV(-matrix(2,0),-matrix(2,1),-matrix(2,2));

    // We want to find the vector v such that matrix * v = (0,1,0)
    view_up=TV(matrix(1,0),matrix(1,1),matrix(1,2));

    // We want to find the vector v such that matrix * v = (-1,0,0)
    // (negative to compensate for the glScalef(-1,1,1) to get the LHS
    view_right=TV(-matrix(0,0),-matrix(0,1),-matrix(0,2));
}

//#####################################################################
// Function Set_View_Frame
//   Assumes vectors are normalized!
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_View_Frame(const TV &view_forward,const TV &view_up,const TV &view_right)
{
    arcball_matrix=MATRIX<T,4>::Identity_Matrix();
    rotation_matrix=MATRIX<T,4>(-view_right.x,view_up.x,-view_forward.x,0,
                                    -view_right.y,view_up.y,-view_forward.y,0,
                                    -view_right.z,view_up.z,-view_forward.z,0,0,0,0,1);
}
//#####################################################################
// Function Get_Look_At
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Get_Look_At(TV &camera,TV &target,TV &up)
{
    TV view_forward,view_right;
    Get_View_Frame(view_forward,up,view_right);

    target=target_position;
    camera=target_position-camera_distance*view_forward;
}
//#####################################################################
// Function Get_Look_At
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Look_At(const TV &camera,const TV &target,const TV &up)
{
    TV view_forward=target-camera;
    TV view_up=(up-up.Projected(view_forward)).Normalized();
    TV view_right=TV::Cross_Product(view_up,view_forward);
    view_right.Normalize();

    target_position=target;
    camera_distance=view_forward.Normalize();

    Set_View_Frame(view_forward,view_up,view_right);

    nearclip=nearclip_factor*camera_distance;
    farclip=farclip_factor*camera_distance;
}
//#####################################################################
// Function Save_View
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Save_View(const std::string& filename,bool verbose)
{
    TV camera,target,up;
    Get_Look_At(camera,target,up);

    if(verbose)
    {
        LOG::cout<<"Writing view to "<<filename<<std::endl;
        LOG::cout<<"Camera: "<<camera<<std::endl;
        LOG::cout<<"Target: "<<target<<std::endl;
        LOG::cout<<"Up: "<<up<<std::endl;
    }

    std::ofstream output(filename.c_str());
    output<<camera<<std::endl<<target<<std::endl<<up<<std::endl;
    output.close();

    output.open((filename+"_render").c_str());
    output<<"\tLocation=\t"<<camera<<std::endl;
    output<<"\tLook_At=\t"<<target<<std::endl;
    output<<"\tPseudo_Up=\t"<<up<<std::endl;
    output<<"\tField_Of_View=\t"<<2*180/pi*atan(tan(0.5*fovy*pi/180.0)*window->Width()/window->Height())<<std::endl; // convert fovy to fovx
    output<<"\tFocal_Distance=\t.1"<<std::endl;
    output<<"\tAspect_Ratio=\t"<<(T)window->Width()/window->Height()<<std::endl;
    output.close();
}
//#####################################################################
// Function Load_View
//#####################################################################
template<class T> bool OPENGL_WORLD<T>::
Load_View(const std::string& filename,bool verbose)
{
    std::ifstream input(filename.c_str());
    if(input)
    {
        TV camera,target,up(0,1,0);
        // Older camera scripts may be missing up vector
        input>>camera>>target>>up;
        if(verbose)
        {
            LOG::cout<<"Loading view from "<<filename<<std::endl;
            LOG::cout<<"Camera: "<<camera<<std::endl;
            LOG::cout<<"Target: "<<target<<std::endl;
            LOG::cout<<"Up: "<<up<<std::endl;
        }
        Set_Look_At(camera,target,up);
        return true;
    }
    else
    {
        if(verbose) LOG::cout<<"Saved view not found "<<filename<<std::endl;
        return false;
    }
}
//#####################################################################
// Function Ray_Through_Normalized_Image_Coordinate
//#####################################################################
template<class T> RAY<VECTOR<T,3> > OPENGL_WORLD<T>::
Ray_Through_Normalized_Image_Coordinate(VECTOR<T,2> coordinates)
{
    MATRIX<T,4> matrix=arcball_matrix*rotation_matrix;
    TV view_forward(-matrix(3,1),-matrix(3,2),-matrix(3,3));
    TV X(coordinates.x*(T)window->Width()/(T)window->Height(),coordinates.y,-1/tan(.5f*fovy*(T)pi/180));
    TV position=target_position-camera_distance*view_forward;
    return RAY<VECTOR<T,3> >(position,matrix.Transposed().Homogeneous_Times(X));
}
//#####################################################################
// Function Set_Left_Handed
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Set_Left_Handed(const bool left_handed)
{
    // To get left handed coordinates, we apply a 180 degree rotation here and do glScalef(-1,1,1) in Render_World
    if(left_handed_coordinate_system!=left_handed){
        rotation_matrix=MATRIX<T,4>::Rotation_Matrix_Y_Axis(pi)*rotation_matrix;
        left_handed_coordinate_system=left_handed;}
}
//#####################################################################
// Function Get_Camera_Position
//#####################################################################
template<class T> VECTOR<T,3> OPENGL_WORLD<T>::
Get_Camera_Position()
{
    TV view_forward,view_up,view_right;
    Get_View_Frame(view_forward,view_up,view_right);

    return target_position-camera_distance*view_forward;
}
//#####################################################################
// Function Get_Target_Position
//#####################################################################
template<class T> VECTOR<T,3> OPENGL_WORLD<T>::
Get_Target_Position()
{
    return target_position;
}
//#####################################################################
// Function Scene_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_WORLD<T>::
Scene_Bounding_Box()
{
    RANGE<TV> bounding_box;
    bool first=true;
    for(int i=0;i<object_list.m;i++) if(use_bounding_box(i) && object_list(i)->Use_Bounding_Box()){
        if(first){bounding_box=object_list(i)->Bounding_Box();first=false;}
        else bounding_box=RANGE<TV>::Combine(bounding_box,object_list(i)->Bounding_Box());}
    return bounding_box;
}
//#####################################################################
// Function Center_Camera_On_Box
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Center_Camera_On_Bounding_Box(const RANGE<TV>& bounding_box,const bool adjust_distance)
{
    target_position=bounding_box.Center();
    if(!adjust_distance) return;

    // Ensure viewing frustum includes bounding box
    if(bounding_box.Edge_Lengths().x*window->Height()/window->Width() > bounding_box.Edge_Lengths().y)
        camera_distance=.5*bounding_box.Edge_Lengths().x*window->Height()/(window->Width()*tan(.5f*fovy*(T)pi/180)); // wide
    else
        camera_distance=.5*bounding_box.Edge_Lengths().y/tan(.5f*fovy*(T)pi/180); // tall
    camera_distance+=.5*bounding_box.Edge_Lengths().z;
    camera_distance=max(camera_distance,(T).1);

    nearclip=nearclip_factor*camera_distance;
    farclip=farclip_factor*camera_distance;
}
//#####################################################################
// Function Reset_Camera_Orientation
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Reset_Camera_Orientation(const bool reset_up_vector_only)
{
    arcball_matrix=MATRIX<T,4>::Identity_Matrix();
    if(reset_up_vector_only){
        TV camera,target,up;Get_Look_At(camera,target,up);Set_Look_At(camera,target,TV(0,1,0));}
    else rotation_matrix=MATRIX<T,4>::Rotation_Matrix_Y_Axis(left_handed_coordinate_system?pi:0);
}
//#####################################################################
// Function Center_Camera_On_Scene
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Center_Camera_On_Scene()
{
    RANGE<TV> bounding_box=Scene_Bounding_Box();
    if(bounding_box.Empty()) bounding_box=RANGE<TV>::Centered_Box();
    Center_Camera_On_Bounding_Box(bounding_box,true);
}
//#####################################################################
// Function Convert_Mouse_Coordinates
//#####################################################################
// This function converts mouse space pixel coordinates to the normalize coordinates the arcball expects
template<class T> VECTOR<T,2> OPENGL_WORLD<T>::
Convert_Mouse_Coordinates(int x,int y)
{
    if(left_handed_coordinate_system) x=window->Width()-x-1;
    VECTOR<T,2> coord;
    if(window->Width()>=window->Height()) {
        coord.x=((T)x/window->Width()-0.5)*2*(window->Width()/window->Height());
        coord.y=-((T)y/window->Height()-0.5)*2;
    } else {
        coord.x=((T)x/window->Width()-0.5)*2;
        coord.y=-((T)y/window->Height()-0.5)*2*(window->Height()/window->Width());
    }
    return coord;
}

//#####################################################################
// Function Save_Screen
//#####################################################################
template<class T> void OPENGL_WORLD<T>::
Save_Screen(const std::string& filename,const bool use_back_buffer,int jpeg_quality)
{
    ARRAY<VECTOR<T,4>,VECTOR<int,2> > image;
    Get_Image(image,use_back_buffer);
    IMAGE<T>::Write(filename,image);
}
//#####################################################################
// Function Get_Image
//#####################################################################
template<class T> template<int d> void OPENGL_WORLD<T>::
Get_Image(ARRAY<VECTOR<T,d>,VECTOR<int,2> > &image,const bool use_back_buffer)
{
    // Assuming GLubyte is same type as unsigned char
    STATIC_ASSERT(sizeof(GLfloat)==sizeof(float));

    // Sanity check: width and height should match current viewport settings
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT,vp);
    PHYSBAM_ASSERT(window->Width()==vp[2] && window->Height()==vp[3]);

    ARRAY<VECTOR<float,d>,VECTOR<int,2> > temporary_image(0,window->Height(),0,window->Width());
    image.Resize(0,window->Width(),0,window->Height()); // temporary is row major
    glReadBuffer(use_back_buffer?GL_BACK:GL_FRONT);
    glReadPixels(0,0,window->Width(),window->Height(),d==3?GL_RGB:GL_RGBA,GL_FLOAT,temporary_image.array.Get_Array_Pointer());
    for(int i=0;i<window->Width();i++)
        for(int j=0;j<window->Height();j++)
            image(i,j)=VECTOR<T,d>(temporary_image(j,i)); // swap to column major
}
//#####################################################################
// Function Display_Prompt_Strings
//#####################################################################
template<class T> void OPENGL_WORLD<T>::Display_Prompt_Strings()
{
    static OPENGL_COLOR prompt_color=OPENGL_COLOR::Red();

    ARRAY<std::string> strings;
    strings.Append(prompt+" "+prompt_response);
    Display_Strings(strings,prompt_color,true);
}
//#####################################################################
// Function Prompt_User
//#####################################################################
template<class T> void OPENGL_WORLD<T>::Prompt_User(const std::string& prompt_input,OPENGL_CALLBACK prompt_response_cb_input,const std::string& default_response)
{
    prompt=prompt_input;
    prompt_response=default_response;
    prompt_response_success=true;
    prompt_response_cb=prompt_response_cb_input;

    glPushAttrib(GL_COLOR_BUFFER_BIT);
    glDrawBuffer(GL_FRONT);
    prompt_mode=true;
    view_auto_help=false; // force it off
}
//#####################################################################
template void OPENGL_WORLD<float>::Get_Image(ARRAY<VECTOR<float,3>,VECTOR<int,2> > &image,const bool use_back_buffer);
template void OPENGL_WORLD<float>::Get_Image(ARRAY<VECTOR<float,4>,VECTOR<int,2> > &image,const bool use_back_buffer);
template void OPENGL_WORLD<double>::Get_Image(ARRAY<VECTOR<double,3>,VECTOR<int,2> > &image,const bool use_back_buffer);
template void OPENGL_WORLD<double>::Get_Image(ARRAY<VECTOR<double,4>,VECTOR<int,2> > &image,const bool use_back_buffer);
}
namespace PhysBAM{
template class OPENGL_WORLD<double>;
template class OPENGL_WORLD<float>;
}
