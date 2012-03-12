//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EPS_FILE_GEOMETRY__
#define __EPS_FILE_GEOMETRY__

#include <PhysBAM_Tools/Images/EPS_FILE.h>
#include <string>
namespace PhysBAM{

template<class T> class TRIANGLE_2D;
template<class T> class SEGMENT_2D;
template<class TV> class SPHERE;

template<class T>
class EPS_FILE_GEOMETRY:public EPS_FILE<T>
{
    typedef VECTOR<T,2> TV;
    typedef EPS_FILE<T> BASE;
public:
    using BASE::Line_Color;using BASE::stream;using BASE::Draw_Object;

    EPS_FILE_GEOMETRY(const std::string& filename,const RANGE<TV>& box=RANGE<TV>(TV(),TV(500,500)))
        :EPS_FILE<T>(filename,box)
    {}

    virtual ~EPS_FILE_GEOMETRY()
    {}

    template<class T_OBJECT>
    void Draw_Object_Colored(const T_OBJECT& object,const VECTOR<T,3>& color)
    {(*stream)<<"gsave ";Line_Color(color);Draw_Object(object);(*stream)<<"grestore"<<std::endl;}

    void Fill_Object(const TRIANGLE_2D<T>& tri, std::string color)
    {(*stream)<<color<<" setrgbcolor"<<std::endl<<
            "newpath"<<std::endl<<
            tri.X(0).x<<" "<<tri.X(0).y<<" moveto"<<std::endl<<
            tri.X(1).x<<" "<<tri.X(1).y<<" lineto"<<std::endl<<
            tri.X(2).x<<" "<<tri.X(2).y<<" lineto"<<std::endl<<
            "closepath"<<std::endl<<
            "fill"<<std::endl<<
            "0 0 0 setrgbcolor"<<std::endl<<
            std::endl;}
    
    void Draw_Object(const TRIANGLE_2D<T>& tri)
    {Draw_Line(tri.X[0],tri.X[1]);Draw_Line(tri.X[1],tri.X[2]);Draw_Line(tri.X[2],tri.X[0]);}

    void Draw_Object(const SEGMENT_2D<T>& seg)
    {Draw_Line(seg.x1,seg.x2);}

    void Draw_Object(const SPHERE<TV>& circle)
    {Emit_Point(circle.center);(*stream)<<circle.radius<<" 0 360 arc stroke"<<std::endl;Bound(circle.center-circle.radius);Bound(circle.center+circle.radius);}

//#####################################################################
};
}
#endif
