//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VIDEO_WINDOW__
#define __VIDEO_WINDOW__

#ifdef WIN32
#include <windows.h>
#endif
#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glu.h>
#else
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/clamp.h>
#include <Tools/Parsing/STRING_UTILITIES.h>
#include <Tools/Utilities/TIMER.h>
#include <string>
#include "VIDEO.h"
#include "VIDEO_VIEW.h"
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>

namespace PhysBAM{

class VIDEO_WINDOW:public Fl_Gl_Window
{
private:
    int frame;
    int next_frame;
    int in_frame;
    int out_frame;
    int frame_rate;
    bool playing;
    bool looping;
    ARRAY<VIDEO_VIEW<float>*> views;
    ARRAY<int> view_order;
    ARRAY<std::string> view_names;
    int timer_id;
    int selected;
    ARRAY<std::string> dummy_names;
    ARRAY<VECTOR<int,2> > dummy_sequences;

    void (*update_ui_callback)(void*);
    void* update_ui_callback_data;

public:
    float frames_per_second;
    float frames_per_second_accumulator;
    int update_count;

    VIDEO_WINDOW(int x_input,int y_input,int w_input,int h_input,const char *l_input=0)
        :Fl_Gl_Window(x_input,y_input,w_input,h_input,l_input),next_frame(-1),frame_rate(24),playing(true),looping(true),selected(0),frames_per_second_accumulator(0),update_count(0)
    {
        timer_id=TIMER::Singleton()->Register_Timer();
    }

    ~VIDEO_WINDOW()
    {
        TIMER::Singleton()->Release_Timer(timer_id);
        views.Delete_Pointers_And_Clean_Memory();
    }

    bool Valid()
    {return views.m>0;}

    void Set_Callback(void (*f)(void*),void* data)
    {update_ui_callback=f;update_ui_callback_data=data;}
    
    bool Is_Movie(const std::string& filename) const
    {const char* extensions[]={".avi",".mov",".qt",".asf",".wmv",".mpg",".mpeg4",".mp4",0};
    for(int i=0;extensions[i]!=0;i++) if(boost::iends_with(filename,extensions[i])) return true;
    return false;}

    bool Is_Sequence(const std::string& filename,std::string& pattern,int& in,int& out) const
    {static boost::regex re("^(.*\\.)(-?[0-9]+)+-([0-9]+)(..+)$");
    boost::cmatch matches;
    if(boost::regex_match(filename.c_str(),matches,re)){
        in=atoi(matches[2].str().c_str());out=atoi(matches[3].str().c_str());pattern=matches[1]+"%d"+matches[4];
        return true;}
    return false;}

    VIDEO_READER<float>* Open_Video(std::string filename){
        VIDEO_READER<float>* video=0;std::string pattern;int in=0,out=0;
        if(Is_Movie(filename)){
            LOG::cout<<"    Reading Video "<<filename<<std::endl;
            video=VIDEO_READER<float>::Read(filename);}
        else if(Is_Sequence(filename,pattern,in,out)){
            LOG::cout<<"    Reading Image Sequence "<<filename<<std::endl;
            video=VIDEO_READER<float>::Read(pattern,in,out);}
        else{
            LOG::cout<<"    Reading Image "<<filename<<std::endl;
            video=VIDEO_READER<float>::Read(filename,1,1);}
        return video;
    }

    std::string Load_Files(ARRAY<ARRAY<std::string> > filenames_of_sequences)
    {
        for(int sequence_index=0;sequence_index<filenames_of_sequences.m;sequence_index++){
            LOG::cout<<"Loading Sequence "<<sequence_index<<std::endl;
            const ARRAY<std::string>& filenames=filenames_of_sequences(sequence_index);
            ARRAY<VIDEO_READER<float>*> videos;
            ARRAY<std::string> names;
            for(int i=0;i<filenames.m;i++){
                VIDEO_READER<float>* video=Open_Video(filenames(i));
                if(video){
                    names.Append(filenames(i));
                    videos.Append(video);}}
            int view_num=views.Append(new VIDEO_VIEW<float>(videos,names));
            view_order.Append(view_num);
            view_names.Append(str(boost::format("Seq %d")%sequence_index));}
        if(views.m>=1){
            frame_rate=views(1)->Frame_Rate();
            if(filenames_of_sequences.m) selected=1;
            frame=in_frame=1;out_frame=views(1)->Frame_Count();
            for(int i=0;i<views.m;i++) Center_Video(i);
            Frame(In());
            damage(1);}
        return Valid()?"Loaded":"[No Video Loaded]";}

    std::string Add_File(std::string filename,const bool new_sequence=false){
        if(!views.m || new_sequence){
            VIDEO_READER<float>* video=Open_Video(filename);
            if(video){
                ARRAY<VIDEO_READER<float>*> videos;videos.Append(video);
                ARRAY<std::string> names;names.Append(filename);
                int num=views.Append(new VIDEO_VIEW<float>(videos,names));
                view_names.Append(str(boost::format("Seq %d")%num));
                view_order.Append(num);
                frame_rate=views(1)->Frame_Rate();
                frame=in_frame=1;out_frame=views(1)->Frame_Count();
                for(int i=0;i<views.m;i++) Center_Video(i);
                Frame(In());
                damage(1);}
            else return "Failed to open video";}
        else if(selected){
            VIDEO_READER<float>* video=Open_Video(filename);
            if(video) views(selected)->Add_Video(video,filename);
            else return "Failed to open video";}
        else return "No selected view";
        return "";
    }

    void Resize_Parent(Fl_Window* window)
    {int extra_w=window->w()-w(),extra_h=window->h()-h();
    window->resize(window->x(),window->y(),views(1)->Width()+extra_w,views(1)->Height()+extra_h);
    for(int i=0;i<views.m;i++){Zoom_Video(i);Center_Video(i);}}

    void Center_Video(const int view_index)
    {views(view_index)->transform.Column(4)=(.5f*VECTOR<float,2>(w(),h())).Append(0).Append(1);damage(1);}
    
    void Zoom_Video(const int view_index)
    {views(view_index)->transform(1,1)=views(view_index)->transform(2,2)=views(view_index)->transform(3,3)=(float)w()/views(view_index)->Width();damage(1);}

    void Zoom_Actual(const int view_index)
    {views(view_index)->transform(1,1)=views(view_index)->transform(2,2)=views(view_index)->transform(3,3)=1;damage(1);}

    void Center_Video_Current()
    {if(selected) Center_Video(selected);}
    
    void Zoom_Video_Current()
    {if(selected) Zoom_Video(selected);}

    void Zoom_Actual_Current()
    {if(selected) Zoom_Actual(selected);}

    //#####################################################################
    // Function draw
    //#####################################################################
    void draw()
    {if(!valid()){glDisable(GL_DEPTH_TEST);glViewport(0,0,w()-1,h()-1);glEnable(GL_TEXTURE_2D);}
    glMatrixMode(GL_PROJECTION);glLoadIdentity();glMatrixMode(GL_MODELVIEW);glLoadIdentity();
    gluOrtho2D(0,w()-1,0,h()-1);
    glClearColor(.1,.1,.1,1);glClear(GL_COLOR_BUFFER_BIT);
    if(Valid()) for(int i=view_order.m;i>=1;i--){int v=view_order(i);views(v)->Draw(v==selected);}}

    //#####################################################################
    // Function handle
    //#####################################################################
    int handle(int event)
    {static bool translating=false;static bool scaling=false;static VECTOR<int,2> origin;static VIDEO_VIEW<float>* view=0;
    if(event==FL_MOUSEWHEEL){
        //if(Fl::event_dy()>0) video_view->transform=video_view->transform*MATRIX<float,4>::Scale_Matrix(1.25f*VECTOR<float,3>(1,1,1));
        //else video_view->transform=video_view->transform*MATRIX<float,4>::Scale_Matrix(.8f*VECTOR<float,3>(1,1,1));
        //damage(1);
        return 1;}
    else if(event==FL_PUSH){
        origin=VECTOR<int,2>(Fl::event_x(),h()-Fl::event_y()-1);
        view=0;
        for(int i=0;i<view_order.m;i++){int v=view_order(i);LOG::cout<<"v is "<<v<<std::endl;
            if(views(v)->Inside(VECTOR<float,2>(origin))){
                LOG::cout<<"i am inside"<<std::endl;
                view=views(v);Select_View(v);update_ui_callback(update_ui_callback_data);break;}}
        if(view==0){selected=0;return 1;}
        if(Fl::event_button()==2)translating=true;
        else if(Fl::event_button()==3)scaling=true;
        return 1;}
    else if(event==FL_DRAG){
        if(Fl::event_button()==2&&translating){
            VECTOR<float,2> translation=VECTOR<float,2>(1,1)*(VECTOR<float,2>(Fl::event_x(),h()-Fl::event_y()-1)-VECTOR<float,2>(origin));
            view->transform=MATRIX<float,4>::Translation_Matrix(translation.Append(0))*view->transform;}
        else if(Fl::event_button()==3&&scaling){
            float scale=pow(2.f,.01f*(origin.y-(h()-Fl::event_y()-1)));
            MATRIX<float,4> translate=MATRIX<float,4>::Translation_Matrix(view->transform.Column(4).Remove_Index(1)-.5f*VECTOR<float,3>(w(),h(),0));
            MATRIX<float,4> translate_inverse=MATRIX<float,4>::Translation_Matrix(-view->transform.Column(4).Remove_Index(1)+.5f*VECTOR<float,3>(w(),h(),0));
            MATRIX<float,4> center=MATRIX<float,4>::Translation_Matrix(.5f*VECTOR<float,3>(w(),h(),0));
            MATRIX<float,4> center_inverse=MATRIX<float,4>::Translation_Matrix(-.5f*VECTOR<float,3>(w(),h(),0));
            //view->transform=view->transform*MATRIX<float,4>::Scale_Matrix(scale*VECTOR<float,3>(1,1,1))*translate*MATRIX<float,4>::Scale_Matrix(scale*VECTOR<float,3>(1,1,1))*translate_inverse;
            //view->transform=center*MATRIX<float,4>::Scale_Matrix(scale*VECTOR<float,3>(1,1,1))*center_inverse*view->transform;
            view->transform=view->transform*MATRIX<float,4>::Scale_Matrix(scale*VECTOR<float,3>(1,1,1));
        }

        origin=VECTOR<int,2>(Fl::event_x(),h()-Fl::event_y()-1);
        damage(1);
        return 1;}
    else if(event==FL_RELEASE){
        translating=false;scaling=false;
        return 1;}
        return Fl_Gl_Window::handle(event);}

    void Update_Frame()
    {if(Valid()){
        if(next_frame!=-1){
            Frame(next_frame);}
        else if(playing){
            if(Frame()==Out() && looping) next_frame=In();
            else next_frame=Frame()+1;
            frames_per_second_accumulator+=(float)1000/TIMER::Singleton()->Peek_And_Reset_Time(timer_id);
            if(++update_count>30){frames_per_second=frames_per_second_accumulator/update_count;update_count=0;frames_per_second_accumulator=0;}
            Frame(next_frame);}}}

    std::string Title() const
    {std::string title="";for(int i=0;i<views.m;i++) title+=" Seq: "+views(i)->Name(); return title;}


    //#####################################################################
    // Query Frames and Ranges
    //#####################################################################
    int Frame() const
    {return frame;}
    
    int In() const
    {return in_frame;}

    int Out()const
    {return out_frame;}
    
    int Min() const
    {return 1;}
    
    int Max() const
    {return views.m?views(1)->Frame_Count():1;}

    int Frame_Rate()
    {return frame_rate;}

    bool Playing()
    {return playing;}

    const ARRAY<std::string>& View_Names()
    {return view_names;}

    int Selected_View() const
    {return selected;}

    const void Select_View(const int i)
    {int j=view_order.Find(i);
    LOG::cout<<"view order "<<view_order<<std::endl;;
    for(int k=j-1;k>=1;k--) view_order(k+1)=view_order(k);
    view_order(1)=i;
    LOG::cout<<"view order "<<view_order<<std::endl;;
    selected=clamp(i,1,views.m);damage(1);}

    const ARRAY<VECTOR<int,2> >& Video_Sequences()
    {if(views.m>0){
        int active=selected?selected:1;
        return views(active)->video_to_frames;}
    else return dummy_sequences;}

    const ARRAY<std::string>& Video_Names()
    {if(views.m>0){
        int active=selected?selected:1;
        return views(active)->names;}
    else return dummy_names;}

    //#####################################################################
    // Choose Ranges
    //#####################################################################
    void Loop_Next(const int offset)
    {if(selected){
        VECTOR<int,2> range(views(selected)->Range(offset));
        LOG::cout<<"range is "<<range<<std::endl;
        In(range[1]);Out(range[2]);}}
    
    void Loop_All()
    {In(Min());Out(Max());}

    //#####################################################################
    // Set Frame Ranges
    //#####################################################################
    void Frame_Rate(int frame_rate_input)
    {frame_rate=frame_rate_input;}
    
    void Lazy_Frame(int frame_number_input)
    {next_frame=frame_number_input;}

    void Frame(int frame_number_input)
    {if(Valid()){
         frame=clamp(frame_number_input,In(),Out());for(int i=0;i<views.m;i++) views(i)->Set_Frame(frame);damage(1);next_frame=-1;}
    }
    
    void In(const int in_frame_input)
    {if(Valid()){in_frame=clamp(in_frame_input,Min(),Max());Frame(frame);}}
    
    void Out(const int out_frame_input)
    {if(Valid()){out_frame=clamp(out_frame_input,Min(),Max());Frame(frame);}}

    void Play()
    {playing=!playing;}
    
    void Stop()
    {playing=false;Frame(In());};

//#####################################################################
};
}
#endif
