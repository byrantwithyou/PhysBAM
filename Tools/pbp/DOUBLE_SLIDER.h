//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DOUBLE_SLIDER__
#define __DOUBLE_SLIDER__

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
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <string>
#include "VIDEO.h"
#include "VIDEO_VIEW.h"
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Gl_Window.H>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>

template<class T>
class DOUBLE_SLIDER:public Fl_Widget
{
private:
    PhysBAM::GRID<PhysBAM::VECTOR<T,1> > pixel_grid;
    PhysBAM::BOX<PhysBAM::VECTOR<T,1> > range;
    PhysBAM::BOX<PhysBAM::VECTOR<T,1> > subrange;
    int pixel_tolerance;
    T tolerance;

public:
    DOUBLE_SLIDER(int x_input,int y_input,int w_input,int h_input,const char *l_input=0)
        :Fl_Widget(x_input,y_input,w_input,h_input,l_input),pixel_tolerance(3)
    {box(FL_DOWN_BOX);Range(0,1);Subrange(.25,.75);}

    ~DOUBLE_SLIDER()
    {}

    void Range(const T min,const T max)
    {range=PhysBAM::BOX<PhysBAM::VECTOR<T,1> >(min,max);pixel_grid.Initialize(w()-Fl::box_dw(box()),range.xmin,range.xmax);tolerance=pixel_tolerance*pixel_grid.dx;Subrange(subrange.xmin,subrange.xmax);damage(1);}

    void Subrange(const T min,const T max)
    {subrange=PhysBAM::BOX<PhysBAM::VECTOR<T,1> >(PhysBAM::clamp(min,range.xmin,range.xmax),PhysBAM::clamp(max,range.xmin,range.xmax));damage(1);}

    PhysBAM::BOX<PhysBAM::VECTOR<T,1> > Subrange() const
    {return subrange;}

    int handle(int event)
    {T value=pixel_grid.x(Fl::event_x()-x()-Fl::box_dx(box()));static enum {BIND_NONE,BIND_LEFT,BIND_CENTER,BIND_RIGHT} bind_location=BIND_NONE;static T old_value;
    switch(event){
      case FL_PUSH:
        if(range.Lazy_Inside(value)){
            if(Fl::event_clicks()){Subrange(range.xmin,range.xmax);Update();return 1;}
            else if(std::abs(value-subrange.xmin)<tolerance) bind_location=BIND_LEFT;
            else if(std::abs(value-subrange.xmax)<tolerance) bind_location=BIND_RIGHT;
            else if(subrange.Lazy_Inside(value)) bind_location=BIND_CENTER;
            old_value=value;}
        return 1;
        break;
      case FL_DRAG:
        if(bind_location!=BIND_NONE){
            T delta=value-old_value;old_value=value;
            T new_xmin=subrange.xmin,new_xmax=subrange.xmax;
            if(bind_location == BIND_CENTER){
                if(delta>0){T newsub=subrange.xmax+delta;delta=PhysBAM::clamp_max(newsub,range.xmax)-subrange.xmax;}
                if(delta<0){T newsub=subrange.xmin+delta;delta=PhysBAM::clamp_min(newsub,range.xmin)-subrange.xmin;}}
            if(bind_location != BIND_LEFT) new_xmax+=delta;
            if(bind_location != BIND_RIGHT) new_xmin+=delta;
            Subrange(new_xmin,new_xmax);
            Update();}
        return 1;
        break;
      case FL_MOVE:
      case FL_ENTER:
        if(std::abs(value-subrange.xmin)<tolerance || std::abs(value-subrange.xmax)<tolerance) fl_cursor(FL_CURSOR_WE);
        else fl_cursor(FL_CURSOR_DEFAULT);
        return 1;
        break;
      case FL_LEAVE:
        fl_cursor(FL_CURSOR_DEFAULT);
        return 1;
        break;
      case FL_RELEASE:
        if(bind_location != BIND_NONE){
            T delta=value-old_value;
            T new_xmin=subrange.xmin,new_xmax=subrange.xmax;
            if(bind_location != BIND_LEFT) new_xmax+=delta;
            if(bind_location != BIND_RIGHT) new_xmin+=delta;
            Subrange(new_xmin,new_xmax);
            Update();
            bind_location=BIND_NONE;}
        return 1;
      default:
        return Fl_Widget::handle(event);}}

protected:
    void draw()
    {
        int x = this->x()+Fl::box_dx(box());
        int y = this->y()+Fl::box_dy(box());
        int w = this->w()-Fl::box_dw(box());
        int h = this->h()-Fl::box_dh(box());

        draw_box();
        int left_pixel=pixel_grid.Cell(subrange.xmin,0).x;
        int pixel_width=pixel_grid.Cell(subrange.xmax,0).x-left_pixel+1;
        draw_box(FL_UP_BOX,left_pixel+x,y,pixel_width,h,selection_color());
        draw_label(x,y,w,h);
    }

    void Update()
    {redraw();do_callback();}
    
//#####################################################################
};
#endif
