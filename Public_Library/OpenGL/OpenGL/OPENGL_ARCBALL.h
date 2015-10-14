//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __H_OPENGL_ARCBALL__
#define __H_OPENGL_ARCBALL__

#include <Tools/Matrices/MATRIX_4X4.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM{

#define LG_NSEGS 4
#define NSEGS (1<<LG_NSEGS)

template<class T> class OPENGL_WORLD;

template<class T>
class OPENGL_ARCBALL:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_WORLD<T>& opengl_world;

    SPHERE<TV> sphere;
    ROTATION<TV> qNow,qDown,qDrag;
    VECTOR<T,2> center,vDown;
    TV vFrom,vTo,vrFrom,vrTo;
    MATRIX<T,4> mNow,mDown,mDeltaNow; //note that these only work in world space
    bool dragging,use_sphere_center,use_object_space;
    OPENGL_WORLD<T>* world;
    int rotation_axis;
    OPENGL_COLOR x_axis,y_axis,z_axis,outer_rim,highlight;

    OPENGL_ARCBALL(STREAM_TYPE stream_type,OPENGL_WORLD<T>& opengl_world_input);
    void Reinitialize();
    void Update(const VECTOR<T,2> &vNow);
    void Update_World(const VECTOR<T,2> &vNow);
    void Update_Obj(const VECTOR<T,2> &vNow);
    MATRIX<T,4> Value() {return mDeltaNow;}

    void Begin_Drag(const VECTOR<T,2> &vNow)
    {dragging=true;vDown=vNow;Update(vNow);}

    void End_Drag(const VECTOR<T,2> &vNow)
    {dragging=false;Update(vNow);qDown=qNow;mDown=mNow;rotation_axis=-1;}
    
    void Display() const
    {glDisable(GL_LIGHTING);DrawOuterRing();DrawResultArc();DrawDragArc();glEnable(GL_LIGHTING);}
    
private:
    void DrawColor(const OPENGL_COLOR &color,int roation_axis,int my_axis) const;
    void DrawAnyArc(const TV &vFrom,const TV &vTo) const
    {assert(vFrom!=vTo);DrawAnyArcObj(vFrom,vTo);return;}
    void DrawAnyArcWorld(const TV &vFrom,const TV &vTo) const;
    void DrawAnyArcObj(const TV &vFrom,const TV &vTo) const;
    void DrawHalfArc(const TV &n) const;
    void DrawOuterRing() const
    {Circ();}
    void DrawDragArc() const;
    void DrawResultArc() const
    {/*RESCOLOR();if(vrFrom!=vrTo) DrawAnyArc(vrFrom,vrTo);*/}
    void Circ() const;
    TV MouseOnSphere(const VECTOR<T,2> &mouse,const VECTOR<T,2> &ballCenter,double ballRadius);
    TV MouseOnSphere(const VECTOR<T,2> &mouse,const SPHERE<TV> &sphere,const TV &previous_result);
    ROTATION<TV> Qt_FromBallPoints(const TV &from,const TV &to);
    void Qt_ToBallPoints(const ROTATION<TV> &q,TV &arcFrom,TV &arcTo);
    TV Bisect_Vectors(const TV &v1,const TV &v2) const;
};
}
#endif
