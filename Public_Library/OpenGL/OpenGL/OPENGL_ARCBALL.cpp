//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL/OPENGL_ARCBALL.h>
#include <OpenGL/OpenGL/OPENGL_WINDOW.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
namespace PhysBAM{
#define LG_NSEGS 4
#define NSEGS (1<<LG_NSEGS)
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_ARCBALL<T>::
OPENGL_ARCBALL(STREAM_TYPE stream_type,OPENGL_WORLD<T>& opengl_world_input)
    :opengl_world(opengl_world_input),dragging(false),rotation_axis(-1)
{
    x_axis=OPENGL_COLOR::Red();
    y_axis=OPENGL_COLOR::Green();
    z_axis=OPENGL_COLOR::Blue();
    outer_rim=OPENGL_COLOR::Ground_Tan();
    highlight=OPENGL_COLOR::Yellow();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_ARCBALL<T>::
Reinitialize()
{
    sphere=SPHERE<TV>();
    qNow=qDown=qDrag=ROTATION<TV>();
    center=vDown=VECTOR<T,2>();
    vFrom=vTo=TV();
}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_ARCBALL<T>::
Update(const VECTOR<T,2> &vNow)
{
    Update_World(vNow);
}
//#####################################################################
// Function Update_World
//#####################################################################
template<class T> void OPENGL_ARCBALL<T>::
Update_World(const VECTOR<T,2> &vNow)
{
    vFrom=MouseOnSphere(vDown,center,sphere.radius);
    vTo=MouseOnSphere(vNow,center,sphere.radius);
    if(dragging){
        qDrag=Qt_FromBallPoints(vFrom,vTo);
        qNow=qDrag*qDown;}
}
//#####################################################################
// Function DrawColor
//#####################################################################
template<class T> void OPENGL_ARCBALL<T>::
DrawColor(const OPENGL_COLOR &color,int roation_axis,int my_axis) const
{
    if(rotation_axis==my_axis) highlight.Send_To_GL_Pipeline();
    else color.Send_To_GL_Pipeline();
}
//#####################################################################
// Function DrawAnyArcWorld
//#####################################################################
template<class T> void OPENGL_ARCBALL<T>::
DrawAnyArcWorld(const TV &vFrom,const TV &vTo) const
{
    assert(vFrom!=vTo);
    int i;
    TV pts[NSEGS+1];
    double dot;
    pts[0]=vFrom;
    pts[1]=pts[NSEGS]=vTo;
    for(i=0;i<LG_NSEGS;i++) pts[1]=Bisect_Vectors(pts[0],pts[1]);
    dot=2*TV::Dot_Product(pts[0],pts[1]);
    for(i=2;i<NSEGS;i++) pts[i]=(pts[i-1]*dot)-pts[i-2];
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glDepthMask(0);
    glDisable(GL_DEPTH_TEST);
    OpenGL_Begin(GL_LINE_STRIP);
    for(i=0;i<NSEGS;i++){
        pts[i].y*=(T)world->window->Width()/(T)world->window->Height();
        OpenGL_Vertex((T)4.8*sphere.radius*pts[i]);}
    OpenGL_End();

    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glDepthMask(1);
    glEnable(GL_DEPTH_TEST);
}
//#####################################################################
// Function DrawAnyArcObj
//#####################################################################
template<class T> void OPENGL_ARCBALL<T>::
DrawAnyArcObj(const TV &vFrom,const TV &vTo) const
{
    assert(vFrom!=vTo);
    int i;
    T radius=sphere.radius*world->Get_Camera_Position().Magnitude();
    TV pts[NSEGS+1];
    double dot;
    pts[0]=vFrom;
    pts[1]=pts[NSEGS]=vTo;
    for(i=0;i<LG_NSEGS;i++) pts[1]=Bisect_Vectors(pts[0],pts[1]);
    dot=2*TV::Dot_Product(pts[0],pts[1]);
    for(i=2;i<NSEGS;i++) pts[i]=(pts[i-1]*dot)-pts[i-2];
    glDepthMask(0);
    glDisable(GL_DEPTH_TEST);
    OpenGL_Begin(GL_LINE_STRIP);
    TV camera=TV(world->Get_Camera_Position()-world->Get_Target_Position()).Normalized();
    for(i=0;i<NSEGS;i++){
        if(qDown.Rotate(pts[i]).Dot(camera)>-1e-8 || qNow.Rotate(pts[i]).Dot(camera)>1e-8)
            OpenGL_Vertex(qNow.Rotate(radius*pts[i]));
        else{
            OpenGL_End();
            OpenGL_Begin(GL_LINE_STRIP);}}
    OpenGL_End();

    glDepthMask(1);
    glEnable(GL_DEPTH_TEST);
}
//#####################################################################
// Function DrawDragArc
//#####################################################################
template<class T> void OPENGL_ARCBALL<T>::
DrawDragArc() const
{
    if(dragging&&vFrom!=vTo){
        OPENGL_COLOR::Gray().Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINE_STRIP);
        OpenGL_Vertex(vFrom);
        OpenGL_Vertex(TV());
        OpenGL_Vertex(vTo);
        OpenGL_End();}
}
//#####################################################################
// Function Circ
//#####################################################################
template<class T> void OPENGL_ARCBALL<T>::
Circ() const
{
    TV p(0,1,0),m(1,0,0),n(0,0,1);
    DrawColor(outer_rim,rotation_axis,4);
    DrawAnyArcWorld(p,m);
    DrawAnyArcWorld(m,-p);
    DrawAnyArcWorld(-p,-m);
    DrawAnyArcWorld(-m,p);
    DrawColor(z_axis,rotation_axis,2);
    DrawAnyArc(p,m);
    DrawAnyArc(m,-p);
    DrawAnyArc(-p,-m);
    DrawAnyArc(-m,p);
    DrawColor(x_axis,rotation_axis,0);
    DrawAnyArc(p,n);
    DrawAnyArc(n,-p);
    DrawAnyArc(-p,-n);
    DrawAnyArc(-n,p);
    DrawColor(y_axis,rotation_axis,1);
    DrawAnyArc(m,n);
    DrawAnyArc(n,-m);
    DrawAnyArc(-m,-n);
    DrawAnyArc(-n,m);
}
//#####################################################################
// Function MouseOnSphere
//#####################################################################
template<class T> auto OPENGL_ARCBALL<T>::
MouseOnSphere(const VECTOR<T,2> &mouse,const VECTOR<T,2> &ballCenter,double ballRadius) -> TV
{
    TV ballMouse;
    T mag;
    ballMouse.x=T((mouse.x-ballCenter.x)/ballRadius);
    ballMouse.y=T((mouse.y-ballCenter.y)/ballRadius);
    mag=ballMouse.Magnitude_Squared();
    if(mag>1){
        T scale=1/sqrt(mag);
        ballMouse.x*=scale;
        ballMouse.y*=scale;
        ballMouse.z=0;}
    else ballMouse.z=sqrt(1-mag);
    return ballMouse;
}
//#####################################################################
// Function Qt_FromBallPoints
//#####################################################################
template<class T> auto OPENGL_ARCBALL<T>::
Qt_FromBallPoints(const TV &from,const TV &to) -> ROTATION<TV>
{
    return ROTATION<TV>::From_Rotated_Vector(from,to);
}
//#####################################################################
// Function Bisect_Vectors
//#####################################################################
template<class T> auto OPENGL_ARCBALL<T>::
Bisect_Vectors(const TV &v1,const TV &v2) const -> TV
{
    TV v=v1+v2;
    T normal=v.Magnitude_Squared();
    if(normal<1e-5) v=TV(0,0,1);
    return v.Normalized();
}
template class OPENGL_ARCBALL<double>;
template class OPENGL_ARCBALL<float>;
}
