//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VIDEO_VIEW__
#define __VIDEO_VIEW__
// Function 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <cstring>
#include "VIDEO.h"
#ifdef WIN32
#include <windows.h>
#endif
#ifndef __APPLE__
#include <GL/gl.h>
#else
#include <OpenGL/gl.h>
#endif

namespace PhysBAM{

template<class T>
class VIDEO_VIEW
{
public:
    MATRIX<T,4> transform;
private:
    ARRAY<VIDEO_READER<T>*> videos;
    ARRAY<int> frame_to_video;
    ARRAY<int> frame_to_local_frame;
public:
    ARRAY<std::string> names;
    ARRAY<VECTOR<int,2> > video_to_frames;
private:

    GLuint texture;
    unsigned char* frame_buffer;
    bool texture_initialized;
    T width_fraction,height_fraction;
    bool need_copy;
    const int width,height;
    const int frame_rate;
    int frame;

public:
    VIDEO_VIEW(const ARRAY<VIDEO_READER<T>*>& videos,const ARRAY<std::string>& names)
        :transform(MATRIX<T,4>::Identity_Matrix()),videos(videos),names(names),texture(-1),frame_buffer(0),texture_initialized(false),need_copy(false),
        width(videos(1)->Width()),height(videos(1)->Height()),frame_rate(videos(1)->Frame_Rate()),frame(1)
    {Update_Video_Pointers();
    frame_buffer=new unsigned char[width*height*3];}

    ~VIDEO_VIEW()
    {delete[] frame_buffer;videos.Delete_Pointers_And_Clean_Memory();}

    void Add_Video(VIDEO_READER<T>* video,std::string name)
    {videos.Append(video);names.Append(name);
    Update_Video_Pointers();}


    void Update_Video_Pointers()
    {frame_to_video.Remove_All();frame_to_local_frame.Remove_All();video_to_frames.Remove_All();
    for(int i=0;i<videos.m;i++){
        VECTOR<int,2> frames(frame_to_video.m+1,0);
        for(int frame=1;frame<=videos(i)->Frame_Count();frame++){frame_to_video.Append(i);frame_to_local_frame.Append(frame);}
        frames[2]=frame_to_video.m;video_to_frames.Append(frames);}}
        
    void Initialize_Texture()
    {glGenTextures(1,&texture);
    glBindTexture(GL_TEXTURE_2D,texture);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    int power_of_two_width=integer_log(width),power_of_two_height=integer_log(height);
    if(1<<power_of_two_width!=width)power_of_two_width++;
    if(1<<power_of_two_height!=height)power_of_two_height++;
    power_of_two_width=1<<power_of_two_width;power_of_two_height=1<<power_of_two_height;
    int texture_dimension=max(power_of_two_width,power_of_two_height);
    int array_size=texture_dimension*texture_dimension*3;
    float* dummy_image=new float[array_size];memset(dummy_image,0,array_size*sizeof(float));
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,texture_dimension,texture_dimension,0,GL_RGB,GL_FLOAT,dummy_image);
    delete [] dummy_image;
    width_fraction=(T)width/texture_dimension;
    height_fraction=(T)height/texture_dimension;
    texture_initialized=true;}

    bool Inside(const VECTOR<T,2>& pixel) const
    {VECTOR<T,2> object_space_pixel=(transform.Inverse().Homogeneous_Times(pixel.Append(0))).Remove_Index(3);
        std::cout<<"pixel "<<pixel<<" object space pixel "<<object_space_pixel<<std::endl;
    return RANGE<VECTOR<T,2> >(-width/2,width/2,-height/2,height/2).Lazy_Inside(object_space_pixel);}

    void Draw(const bool selected)
    {if(!texture_initialized) Initialize_Texture();
    if(need_copy){need_copy=false;
        glPixelStorei(GL_UNPACK_ALIGNMENT,1);glPixelStorei(GL_UNPACK_ALIGNMENT,1);
        glTexSubImage2D(GL_TEXTURE_2D,0,0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,frame_buffer);}
    glPushMatrix();
    glMultMatrixf(transform.x);
    //glTranslatef(origin.x,origin.y,0);
    //glScalef(scale,scale,1);
    glTranslatef(-(T)width/2,-(T)height/2,0); 
    glDisable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    if(selected) glColor3f(1,1,0);else glColor3f(.4,.6,.8);
    glVertex3f(-2,-2,0);
    glVertex3f(2+width,-2,0);
    glVertex3f((T)2+width,(T)2+height,0);
    glVertex3f(-2,(T)2+height,0);
    glEnd();
    glColor3f(1,1,1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,texture);
    glBegin(GL_QUADS);
    glTexCoord2f(0,height_fraction);glVertex3f(0,0,0);
    glTexCoord2f(width_fraction,height_fraction);glVertex3f((T)width,0,0);
    glTexCoord2f(width_fraction,0);glVertex3f((T)width,(T)height,0);
    glTexCoord2f(0,0);glVertex3f(0,height,0);
    glEnd();

    glPopMatrix();}

    int Frame_Count() const
    {return frame_to_video.m;}

    int Frame_Rate() const
    {return frame_rate;}

    void Set_Frame(const int frame)
    {if(frame>frame_to_video.m) return;
    this->frame=frame;
    videos(frame_to_video(frame))->Frame(frame_to_local_frame(frame),frame_buffer);need_copy=true;}
    
    int Width() const
    {return width;}
 
    int Height() const
    {return height;}
    
    std::string Name() const
    {return names(frame_to_video(frame));}

    VECTOR<int,2> Range(const int offset) const
    {return video_to_frames((frame_to_video(frame)+offset-1)%video_to_frames.m+1);}

//#####################################################################
};
}
#endif
