//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_EPS_OUTPUT
//#####################################################################
#ifndef __OPENGL_EPS_OUTPUT__
#define __OPENGL_EPS_OUTPUT__
#include <Core/Arrays/ARRAY.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{
template<class T>
class OPENGL_EPS_OUTPUT
{
protected:
    typedef VECTOR<T,3> TV;
    std::ostream& stream;
public:
    int glmode;
    ARRAY<TV> buffer;

    OPENGL_EPS_OUTPUT(const std::string& filename);
    virtual ~OPENGL_EPS_OUTPUT();

    void Begin(int mode);
    void End();
    template<class T2> void Vertex(const VECTOR<T2,3>& p);
    template<class T2> void Vertex(const VECTOR<T2,2>& p);
    template<class T2> void Vertex(const VECTOR<T2,1>& p);
    void Head();
    void Tail();
    template<class T2> void Emit(const VECTOR<T2,3>& p);
    void Emit(const char* p);
    void Set_Color(const TV& color);
    void Draw_Arrays(int mode,int dimension,int length,const void* vertices);
    template<class T2,int d> void Draw_Bezier(const VECTOR<VECTOR<T2,d>,4>& pts);
protected:
    void Draw_Point(const TV& p);
    void Draw_Line(const TV& a,const TV& b);
    void Draw_Polygon(int i,int n);
    void Transform_Buffer();
};
}
#endif
