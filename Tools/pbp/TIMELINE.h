//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TIMELINE__
#define __TIMELINE__

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
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <cstdio>
#include <string>
#include "VIDEO.h"
#include "VIDEO_VIEW.h"
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Gl_Window.H>

using namespace PhysBAM;

template<class T>
class TIMELINE:public Fl_Widget
{
private:
    PhysBAM::INTERVAL<int> frame_range;
    PhysBAM::GRID<VECTOR<T,1> > frame_grid;
    PhysBAM::INTERVAL<int> widget_range;
    PhysBAM::GRID<VECTOR<T,1> > widget_grid;
    PhysBAM::ARRAY<PhysBAM::VECTOR<int,2> > ranges;
    PhysBAM::ARRAY<std::string> names;
    int frame;
    int in,out;
    enum {BIND_NONE,BIND_LEFT,BIND_CENTER,BIND_RIGHT,BIND_FRAME} bind_location;

public:
    TIMELINE(int x_input,int y_input,int w_input,int h_input,const char *l_input=0)
        :Fl_Widget(x_input,y_input,w_input,h_input,l_input),frame(0),in(1),out(1),bind_location(BIND_NONE)
    {box(FL_DOWN_BOX);Range(0,1,0);}

    ~TIMELINE()
    {}

    void Init_Sequence(const PhysBAM::ARRAY<PhysBAM::VECTOR<int,2> >& ranges_input,const ARRAY<std::string>& names_input)
    {ranges=ranges_input;names=names_input;
        LOG::cout<<ranges<<std::endl;
        LOG::cout<<names<<std::endl;
    }

    void Range(const int min,const int max,const int frame_input)
    {frame_range=PhysBAM::INTERVAL<int>(min,max);
    frame_grid.Initialize(max-min+1,0,1,true);frame=frame_input;
    damage(FL_DAMAGE_ALL);}

    void Frame(int frame_raw)
    {frame=PhysBAM::clamp(frame_raw,in,out);}

    int Frame() const
    {return frame;}

    void Subrange(const int min,const int max)
    {in=PhysBAM::clamp(min,frame_range.min_corner,frame_range.max_corner);out=PhysBAM::clamp(max,frame_range.min_corner,frame_range.max_corner);damage(FL_DAMAGE_ALL);}

    VECTOR<int,2> Subrange() const
    {return VECTOR<int,2>(in,out);}

    int Global(const int local)
    {return frame_grid.Index(widget_grid.X(local-widget_range.min_corner)).x;}

    int Local(const int global)
    {return widget_grid.Index(frame_grid.X(global)).x+widget_range.min_corner;}

    int Local_DX(const int global,const T scale=0)
    {return widget_grid.Index(frame_grid.X(global)+frame_grid.dX.x*scale).x+widget_range.min_corner;}

    int handle(int event)
    {int event_frame=PhysBAM::clamp(Global(Fl::event_x()),frame_range.min_corner,frame_range.max_corner);int real_frame=Local(frame);static int saved_frame;
    switch(event){
        case FL_PUSH:
            if(abs(Fl::event_x()-real_frame)<10 || Fl::event_y()-y()<15){
                bind_location=BIND_FRAME;
                Frame(event_frame);
                Update();}
            else if(abs(Fl::event_x()-Local(in))<5){
                bind_location=BIND_LEFT;
                Subrange(min(event_frame,out),out);
                Update();}
            else if(abs(Fl::event_x()-Local(out))<5){
                bind_location=BIND_RIGHT;
                Subrange(in,max(event_frame,in));
                Update();}
            else if(in<event_frame && out>event_frame){
                bind_location=BIND_CENTER;
                saved_frame=event_frame;
                Update();}
            return 1;
            break;
        case FL_DRAG:
            if(bind_location!=BIND_NONE){
                if(bind_location == BIND_FRAME) Frame(event_frame);
                else if(bind_location == BIND_LEFT) Subrange(min(event_frame,out),out);
                else if(bind_location == BIND_RIGHT) Subrange(in,max(event_frame,in));
                else if(bind_location == BIND_CENTER){
                    int offset=event_frame-saved_frame;saved_frame=event_frame;
                    if(offset>0) offset=min(offset,min(out+offset,frame_range.max_corner)-out);
                    else if(offset<0) offset=max(offset,max(in+offset,frame_range.min_corner)-in);
                    Subrange(in+offset,out+offset);}
                Update();
            }
            return 1;
            break;
        case FL_MOVE:
        case FL_ENTER: 
            if(abs(Fl::event_x()-real_frame)<10 || Fl::event_y()-y()<15) fl_cursor(FL_CURSOR_CROSS);
            else if(abs(Fl::event_x()-Local(in))<5) fl_cursor(FL_CURSOR_WE);
            else if(abs(Fl::event_x()-Local(out))<5) fl_cursor(FL_CURSOR_WE);
            else if(in<event_frame && out>event_frame) fl_cursor(FL_CURSOR_HAND);
            else fl_cursor(FL_CURSOR_DEFAULT);
            return 1;
            break;
        case FL_LEAVE:
            fl_cursor(FL_CURSOR_DEFAULT);
            return 1;
        break;

        case FL_RELEASE:
            bind_location=BIND_NONE;
            return 1;
            break;
    }
            
    
    return Fl_Widget::handle(event);}

protected:
    void draw()
    {
        static char buf[10];
        draw_box();
        
        int padding=10;
        widget_range=PhysBAM::INTERVAL<int>(this->x()+Fl::box_dx(box())+padding,this->w()+this->x()-Fl::box_dw(box())-padding);
        widget_grid.Initialize(widget_range.Size()+1,0,1,true);
        int y_offset=y()+Fl::box_dy(box());

        fl_color(130,255,130);
        for(int i=1;i<=ranges.m;i+=2){
            int xmin=Local_DX(ranges(i)[1],(T)-.5),xmax=Local_DX(ranges(i)[2],(T).5);
            //LOG::cout<<ranges(i)<<" xmin = "<<xmin<<" xmax = "<<xmax<<std::endl;
            fl_rectf(xmin,y_offset,xmax-xmin,10);}
        fl_color(100,200,100);
        for(int i=2;i<=ranges.m;i+=2){
            int xmin=Local_DX(ranges(i)[1],(T)-.5),xmax=Local_DX(ranges(i)[2],(T).5);
            fl_rectf(xmin,y_offset,xmax-xmin,10);}

        if(bind_location==BIND_CENTER) fl_color(150,150,0);else fl_color(80,140,200);
        fl_rectf(Local(in),y_offset+30,Local(out)-Local(in),5);
        if(bind_location==BIND_LEFT) fl_color(150,150,0);else fl_color(120,180,255);
        fl_rectf(Local(in)-2,y_offset+25,4,15);
        if(bind_location==BIND_RIGHT) fl_color(150,150,0);else fl_color(120,180,255);
        fl_rectf(Local(out)-2,y_offset+25,4,15);

        if(bind_location==BIND_FRAME) fl_color(150,150,0);else fl_color(0,100,255);
        fl_rectf(Local(frame)-2,y_offset,4,25);
        fl_rectf(Local(frame)-15,y_offset+25,30,15);
        fl_color(255,255,255);
        sprintf(buf,"%d",frame);
        fl_font(FL_COURIER,9);
        fl_draw(buf,Local(frame)-15,y_offset+25,30,15,FL_ALIGN_CENTER);
        fl_color(0,100,0);
        int budget=max(1,frame_grid.counts.x/(w()/10));int count=0;
        for(int i=1;i<=frame_grid.counts.x;i+=budget,count++){
            int offset=5;
            if(count%5==0) offset=0;
            fl_line(Local(i),y_offset+offset,Local(i),y_offset+10);}
        fl_line(Local_DX(1,(T)-.5),y_offset+10,Local_DX(frame_grid.counts.x,(T).5),y_offset+10);
        fl_color(0,0,0);
        for(int i=1;i<=frame_grid.counts.x;i+=5*budget){
            sprintf(buf,"%d",i);
            fl_draw(buf,Local(i),y_offset+20);}
    }

    void Update()
    {redraw();do_callback();}
};
#endif
