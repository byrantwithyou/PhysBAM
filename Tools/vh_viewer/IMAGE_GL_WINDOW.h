//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __IMAGE_GL_WINDOW__
#define __IMAGE_GL_WINDOW__

#include <Fl/Fl.h>
#include <FL/Fl_Counter.h>
#include <FL/Fl_Gl_Window.h>
#include <FL/Fl_Window.h>
#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glu.h>
#else
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TEXTURED_RECT.h>
#include "../vh_segment/VH_LEVELSET_BUILDER.h"
#include "VIEWER_WINDOW.h"

namespace PhysBAM{

template<class T>
void OpenGL_Vertex(VECTOR<T,2>& X)
{glVertex2f(X.x,X.y);}

template<class T> class VIEWER_WINDOW;
template<class T>
class IMAGE_GL_WINDOW:public Fl_Gl_Window
{
public:
    GRID<VECTOR<T,2> > grid;
    VECTOR<int,2> old_event;
    OPENGL_TEXTURED_RECT rectangle;
    OPENGL_TEXTURE texture;
    T scale,scale_factor;VECTOR<T,2> translate;
    VIEWER_WINDOW<T>& parent;
    bool update_image;
    OPENGL_COLOR* opengl_image;
    
    IMAGE_GL_WINDOW(int x, int y, int w, int h,VIEWER_WINDOW<T>& parent_input) 
        :Fl_Gl_Window(x,y,w,h,"My GL Window"),scale(1),scale_factor(.01),translate(0,0),parent(parent_input),update_image(true),opengl_image(0)
    {
        grid.Initialize(w,h,-640,640,-640,640);
    }

    void Update_Image()
    {//parent.visible_human.Convert_To_Image(image);
    if(!opengl_image) opengl_image=new OPENGL_COLOR[parent.visible_human.image_resolution.x*parent.visible_human.image_resolution.z];
    int t=0;for(int ij=0;ij<parent.visible_human.image_resolution.z;ij++) for(int i=0;i<parent.visible_human.image_resolution.x;i++){
        int tissue_id=parent.visible_human.image(i,ij);
        opengl_image[t++]=Color_Map(tissue_id);
        int entry;
        if(parent.tissues_hash.Get(tissue_id,entry)) opengl_image[t-1]=OPENGL_COLOR(1,1,1,1);}
    texture.Update_Texture(opengl_image);rectangle.width=parent.visible_human.image_resolution.x;rectangle.height=parent.visible_human.image_resolution.z;
    //delete opengl_image;
    rectangle.Set_Texture(&texture);}

private:
    void Setup_Viewing_Transform(){
    glMatrixMode(GL_PROJECTION);glLoadIdentity();
    gluOrtho2D(-w()/2,w()/2,-h()/2,h()/2);
    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);glLoadIdentity();
    glScalef(scale,scale,1);
    glTranslatef(translate.x,translate.y,0);}

    void Setup_Texture()
    {if(update_image){
        update_image=false;
        texture.Initialize(parent.visible_human.image_resolution.x,parent.visible_human.image_resolution.z,0);
        Update_Image();}}
    
    void draw() 
    {Setup_Texture();
    Setup_Viewing_Transform();
    glClearColor(.5,.5,.5,0);glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    rectangle.Display();}

    int handle(int event)
    {VECTOR<int,2> index(Fl::event_x()+1,h()-Fl::event_y());
    if(Fl::event_button()==1){
        if(event==FL_PUSH || event==FL_DRAG){
            VECTOR<double,3> value;
            Setup_Viewing_Transform();
            glTranslatef(-rectangle.width/2,-rectangle.height/2,0);
            GLdouble modelview[16],projection[16];GLint viewport[4];
            glGetIntegerv(GL_VIEWPORT,viewport);glGetDoublev(GL_PROJECTION_MATRIX,projection);glGetDoublev(GL_MODELVIEW_MATRIX,modelview);
            gluUnProject(index.x,index.y,1,modelview,projection,viewport,&value.x,&value.y,&value.z);value+=VECTOR<double,3>(1,1,1);
            parent.Update_Hover_Point(VECTOR<int,2>((int)value.x,(int)value.y));
        }
        damage(1);
        return 1;}
    if(Fl::event_button()==2){
        if(event==FL_PUSH){old_event=index;}
        else if(event==FL_DRAG){translate+=(T)1/scale*VECTOR<T,2>(index-old_event);old_event=index;damage(1);}
        return 1;}
    if(Fl::event_button()==3){
        if(event==FL_PUSH){old_event=index;}
        else if(event==FL_DRAG){scale-=(old_event.y-index.y)*scale_factor;old_event=index;damage(1);}
        damage(1);
        return 1;}
    return 0;}

    OPENGL_COLOR Color_Map(unsigned short int image_index)
    {static VECTOR<unsigned short,3> masks(0x001f,0x07e0,0xf800);
    static VECTOR<T,3> maxes=VECTOR<T,3>(1<<5,1<<6,1<<5)-VECTOR<T,3>(1,1,1);
    static VECTOR<T,3> color_scale=VECTOR<T,3>(255,255,255)/maxes;
    static T color_map_factor=(T)1/(T)256;
    VECTOR<T,3> color=color_map_factor*(color_scale*VECTOR<T,3>((image_index&masks.x),(image_index&masks.y)>>5,(image_index&masks.z)>>11));
    color.x-=floor(color.x);color.y-=floor(color.y);color.z-=floor(color.z);
    return OPENGL_COLOR(color.x,color.y,color.z,(T)1);}

//#####################################################################
};
}
#endif
