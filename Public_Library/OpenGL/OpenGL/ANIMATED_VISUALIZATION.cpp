//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Read_Write/STRING_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Images/MOV_FILE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <climits>
#include <cstdlib>
#include <fstream>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> ANIMATED_VISUALIZATION<T>::
ANIMATED_VISUALIZATION(VIEWER_DIR& viewer_dir)
    :viewer_dir(viewer_dir),substeps_level(0),play(false),loop(false),fixed_frame_rate(false),start_frame(0),stop_frame(INT_MAX),
    frame_rate(24),jpeg_quality(95)
{
    next_frame_cb={[this](){Next_Frame();},"Next frame"};
    prev_frame_cb={[this](){Prev_Frame();},"Prev frame"};
    goto_frame_cb={[this](){Goto_Frame();},"Goto frame"};
    reset_cb={[this](){Reset();},"Reset"};
    toggle_play_cb={[this](){Toggle_Play();},"Play/Pause"};
    toggle_loop_cb={[this](){Toggle_Loop();},"Loop"};
    toggle_fixed_frame_rate_cb={[this](){Toggle_Fixed_Frame_Rate();},"Toggle fixed frame rate"};
    goto_last_frame_cb={[this](){Goto_Last_Frame();},"Goto last frame"};
    goto_frame_prompt_cb={[this](){Goto_Frame_Prompt();},"goto_frame_prompt"};
    capture_frames_cb={[this](){Capture_Frames();},"Capture frames (to images)"};
    capture_frames_prompt_cb={[this](){Capture_Frames_Prompt();},"capture_frames_prompt"};
    null_cb={0,0};

    if(MOV_WRITER<float>::Enabled()) saved_frame_filename="capture.mov";
    else if(IMAGE<float>::Is_Supported(".png")) saved_frame_filename="capture.%05d.png";
    else if(IMAGE<float>::Is_Supported(".jpg")) saved_frame_filename="capture.%05d.jpg";
    else saved_frame_filename="capture.%05d.ppm";
}
//#####################################################################
// Function Add_Arguments
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS& parse_args)
{
    BASIC_VISUALIZATION<T>::Add_Arguments(parse_args);

    parse_args.Add("-jpeg_quality",&jpeg_quality,"quality","jpeg quality settings");
    parse_args.Add("-so",&saved_frame_filename,"file","save frames output");
    parse_args.Add("-start_frame",&start_frame,"frame","start frame");
    parse_args.Add("-stop_frame",&stop_frame,"frame","stop frame");
    parse_args.Add("-fps",&frame_rate,"rate","frames per second");
}
//#####################################################################
// Function Parse_Arguments
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS& parse_args)
{
    BASIC_VISUALIZATION<T>::Parse_Arguments(parse_args);
}
//#####################################################################
// Function Add_OpenGL_Initialization
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Add_OpenGL_Initialization()
{
    BASIC_VISUALIZATION<T>::Add_OpenGL_Initialization();

    Update_OpenGL_Strings();
}
//#####################################################################
// Function Initialize_Components_And_Key_Bindings
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Initialize_Components_And_Key_Bindings()
{
    viewer_dir.Advance_Directory(0);
    BASIC_VISUALIZATION<T>::Initialize_Components_And_Key_Bindings();

    opengl_world.Set_Key_Binding_Category("Default Keys (ANIMATED_VISUALIZATION)");

    opengl_world.Bind_Key('p',toggle_play_cb);
    opengl_world.Bind_Key('P',toggle_loop_cb);
    opengl_world.Bind_Key('r',reset_cb);
    opengl_world.Bind_Key('s',next_frame_cb);
    opengl_world.Bind_Key('S',next_frame_cb);
    opengl_world.Bind_Key("^s",prev_frame_cb);
    opengl_world.Bind_Key('g',goto_frame_cb);
    opengl_world.Bind_Key('f',toggle_fixed_frame_rate_cb);
    opengl_world.Bind_Key('z',goto_last_frame_cb);
    opengl_world.Bind_Key("^d",capture_frames_cb);

    opengl_world.Set_Key_Binding_Category("User-Defined Keys");
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool ANIMATED_VISUALIZATION<T>::
Valid_Frame()
{
    ARRAY<int> last_frame_stack;
    viewer_dir.Read_Last_Frame(last_frame_stack);
    return !LEXICOGRAPHIC_COMPARE()(last_frame_stack,viewer_dir.frame_stack);
}
//#####################################################################
// Function Goto_Start_Frame
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Goto_Start_Frame()
{
    viewer_dir.Set(0);
    Set_Frame();
}
//#####################################################################
// Function Goto_Last_Frame
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Goto_Last_Frame()
{
    viewer_dir.Read_Last_Frame();
    if(Directory_Exists(viewer_dir.current_directory))
        Set_Frame();
}
//#####################################################################
// Function Render_Offscreen
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Render_Offscreen()
{
    Capture_Frames(saved_frame_filename,start_frame,stop_frame,jpeg_quality,false);
}
OPENGL_EPS_OUTPUT<float>* PhysBAM::opengl_eps_output=0;
//#####################################################################
// Function Capture_Frames
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Capture_Frames(const std::string& filename_pattern,int capture_start_frame,int capture_end_frame,int jpeg_quality,bool swap_buffers)
{
    bool use_eps=IEnds_With(filename_pattern,".eps");
    bool movie=IEnds_With(filename_pattern,".mov");
    MOV_WRITER<T>* mov=0;
    if(movie && MOV_WRITER<T>::Enabled()){
        mov=new MOV_WRITER<T>(filename_pattern,24);
        LOG::cout<<"Capturing to quicktime '"<<filename_pattern<<"'"<<std::endl;}

    LOG::cout<<"Capturing frames "<<capture_start_frame<<" to ";
    if(capture_end_frame==INT_MAX) LOG::cout<<"last valid frame";else LOG::cout<<capture_end_frame;
    LOG::cout<<" into '"<<filename_pattern<<"'"<<std::endl;
    for(int f=capture_start_frame;f<=capture_end_frame;f++)
    {
        viewer_dir.Set(f);
        if(!Valid_Frame()) break;
        Set_Frame();
        if(use_eps){
            std::string filename=LOG::sprintf(filename_pattern.c_str(),f);
            opengl_eps_output=new OPENGL_EPS_OUTPUT<float>(filename);}
        opengl_world.Render_World(false,swap_buffers);
        glFinish();
        if(mov){
            LOG::cout<<"  Frame "<<f<<std::endl;
            ARRAY<VECTOR<T,3>,VECTOR<int,2> > image;
            opengl_world.Get_Image(image,swap_buffers);
            if(mov) mov->Add_Frame(image);}
        else if(!use_eps){
            std::string filename=LOG::sprintf(filename_pattern.c_str(),f);
            LOG::cout<<"Capturing frame "<<f<<" to "<<filename<<std::endl;
            opengl_world.Save_Screen(filename,swap_buffers,jpeg_quality);}
        else{
            delete opengl_eps_output;
            opengl_eps_output=0;}
    }
    delete mov;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Set_Frame()
{
    Pre_Frame_Extra();

    for(int i=0;i<component_list.m;i++){
#ifdef NDEBUG
        int attempts=0;
        bool done=false;
        while(!done && attempts++ < 10){
            try{
                component_list(i)->Set_Frame();
                done=true;}
            catch(std::exception& error){
                LOG::cerr<<"Read error: "<<error.what()<<std::endl;
                LOG::cerr<<"Retrying in 1/10 s..."<<std::endl;
                PROCESS_UTILITIES::Sleep(.1);}}
#else
        // don't do try-catch in debug to allow gdb to catch it
        component_list(i)->Set_Frame();
#endif
    }

    // Update/Invalidate selections
    if(selected_object && selected_object->Destroy_Selection_After_Frame_Change())
        selected_object=0;

    Set_Frame_Extra();
    Update_OpenGL_Strings();
}
//#####################################################################
// Function Update_OpenGL_Strings
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Update_OpenGL_Strings()
{
    opengl_world.Clear_Strings();
    std::ostringstream stream;
    stream<<"frame ";
    viewer_dir.frame_stack.Write_Raw(stream);
    if(!frame_title.empty()) stream<<": "<<frame_title;
    opengl_world.Add_String(stream.str());
    BASIC_VISUALIZATION<T>::Update_OpenGL_Strings();
}
//#####################################################################
// Function Next_Frame
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Next_Frame()
{
    bool valid=viewer_dir.Find_Next_Directory(substeps_level);
    if(!valid && loop) viewer_dir.Set(0);
    if(valid || loop){
        Set_Frame();
        opengl_world.Set_Idle_Callback(play?next_frame_cb:null_cb,fixed_frame_rate?(float)1/frame_rate:0);}
    else if(play) opengl_world.Set_Idle_Callback(next_frame_cb,.2);
}
//#####################################################################
// Function Prev_Frame
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Prev_Frame()
{
    if(viewer_dir.Find_Prev_Directory(substeps_level))
        Set_Frame();
}
//#####################################################################
// Function Goto_Frame_Prompt
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Goto_Frame_Prompt()
{
    if(opengl_world.prompt_response.empty()) return;
    viewer_dir.Set(opengl_world.prompt_response);
    if(!Directory_Exists(viewer_dir.current_directory)) return;
    Set_Frame();
}
//#####################################################################
// Function Goto_Frame
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Goto_Frame()
{
    opengl_world.Prompt_User("Goto frame: ",goto_frame_prompt_cb,"");
}
//#####################################################################
// Function Reset
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Reset()
{
    viewer_dir.Set(0);
    Set_Frame();
    if(play) Toggle_Play(); // Stop playing
}
//#####################################################################
// Function Toggle_Play
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Toggle_Play()
{
    play=!play;
    opengl_world.Set_Idle_Callback(play?next_frame_cb:null_cb,fixed_frame_rate?(float)1/frame_rate:0);
}
//#####################################################################
// Function Toggle_Loop
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Toggle_Loop()
{
    loop=!loop;
}
//#####################################################################
// Function Toggle_Fixed_Frame_Rate
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Toggle_Fixed_Frame_Rate()
{
    fixed_frame_rate=!fixed_frame_rate;
    if(fixed_frame_rate) last_frame_time=0;
    opengl_world.Set_Idle_Callback(play?next_frame_cb:null_cb,fixed_frame_rate?(float)1/frame_rate:0);
}
//#####################################################################
// Function Capture_Frames_Prompt
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Capture_Frames_Prompt()
{
    if(!opengl_world.prompt_response_success) return;

    bool done=false;
    if(capture_frames_prompt_state.step==0){
        if(!opengl_world.prompt_response.empty()) capture_frames_prompt_state.filename_pattern=opengl_world.prompt_response;
        capture_frames_prompt_state.step=1;
        opengl_world.Prompt_User(LOG::sprintf("Start frame [%d]: ",start_frame),capture_frames_prompt_cb,"");}
    else if(capture_frames_prompt_state.step==1){
        if(!opengl_world.prompt_response.empty()) String_To_Value(opengl_world.prompt_response,capture_frames_prompt_state.start_frame);
        capture_frames_prompt_state.step=2;
        opengl_world.Prompt_User("End frame [last valid frame]: ",capture_frames_prompt_cb,"");}
    else if(capture_frames_prompt_state.step==2){
        if(!opengl_world.prompt_response.empty()) String_To_Value(opengl_world.prompt_response,capture_frames_prompt_state.end_frame);
        capture_frames_prompt_state.step=3;
        if(toupper(Get_File_Extension(capture_frames_prompt_state.filename_pattern))=="JPG")
            opengl_world.Prompt_User(LOG::sprintf("JPEG quality [%d]: ",jpeg_quality),capture_frames_prompt_cb,"");
        else done=true;}
    else if(capture_frames_prompt_state.step==3){
        if(opengl_world.prompt_response.empty()) String_To_Value(opengl_world.prompt_response,capture_frames_prompt_state.jpeg_quality);
        done=true;}

    if(done){
        Capture_Frames(capture_frames_prompt_state.filename_pattern,capture_frames_prompt_state.start_frame,capture_frames_prompt_state.end_frame,
            capture_frames_prompt_state.jpeg_quality);
        if(system(LOG::sprintf("pbp %s &",capture_frames_prompt_state.filename_pattern.c_str()).c_str())){};}
}
//#####################################################################
// Function Capture_Frames
//#####################################################################
template<class T> void ANIMATED_VISUALIZATION<T>::
Capture_Frames()
{
    // Fill with defaults
    capture_frames_prompt_state.filename_pattern=saved_frame_filename;
    capture_frames_prompt_state.start_frame=start_frame;
    capture_frames_prompt_state.end_frame=INT_MAX;
    capture_frames_prompt_state.jpeg_quality=jpeg_quality;
    capture_frames_prompt_state.step=0;
    opengl_world.Prompt_User("Capture filename [" + saved_frame_filename + "]: ",capture_frames_prompt_cb,"");
}
namespace PhysBAM{
template class ANIMATED_VISUALIZATION<double>;
template class ANIMATED_VISUALIZATION<float>;
}
