//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_IMAGE
//#####################################################################
#ifndef __VECTOR_IMAGE__
#define __VECTOR_IMAGE__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <fstream>
#include <string>
namespace PhysBAM{

template<class T>
class VECTOR_IMAGE
{
public:
    typedef VECTOR<T,2> TV;
    std::ofstream stream;
    RANGE<TV> bounding_box;
    RANGE<TV> output_box;
    bool fixed_bounding_box;

    struct FORMATTING
    {
        VECTOR<T,3> line_color,fill_color;
        T line_width,point_radius,line_opacity,fill_opacity;
        int line_style; // 0=none, 1=solid, 2=dotted, 3=dashed
        int fill_style; // 0=none, 1=solid
        const char* arrow_style; // LaTeX
        std::string misc; // LaTeX

        FORMATTING();
    } def_format, cur_format;

    VECTOR_IMAGE(const std::string& filename,const RANGE<TV>& box=RANGE<TV>(TV(),TV(500,500)));
    virtual ~VECTOR_IMAGE();

    void Use_Fixed_Bounding_Box(const RANGE<TV>& box)
    {bounding_box=box;fixed_bounding_box=true;}

    void Bound(const TV& pt);

    void Draw_Object(const TV &a) // point
    {Draw_Object(a,cur_format.point_radius);}
    void Draw_Object(const TV &a,const TV &b); // line
    void Draw_Object(const TV &a,const TV &b,const TV &c); // triangle
    void Draw_Object(ARRAY_VIEW<TV> pts); // polygon
    void Draw_Object(ARRAY_VIEW<TV> outside,ARRAY_VIEW<ARRAY_VIEW<TV> > holes); // polygon with holes
    template <int d> void Draw_Object(const VECTOR<TV,d>& pts) // polygon
    {Draw_Object(ARRAY_VIEW<TV>(d,const_cast<TV*>(&pts(0))));}
    void Draw_Object(const TV &pt,T radius); // circle
    void Draw_Object(const RANGE<TV>& box);

protected:
    virtual void Emit_Object(const TV &a,const TV &b)=0; // line
    virtual void Emit_Object(const TV &a,const TV &b,const TV &c)=0; // triangle
    virtual void Emit_Object(ARRAY_VIEW<TV> pts)=0; // polygon
    virtual void Emit_Object(ARRAY_VIEW<TV> pts,ARRAY_VIEW<ARRAY_VIEW<TV> > holes)=0; // polygon with holes
    virtual void Emit_Object(const TV &pt,T radius)=0; // circle
    virtual void Emit_Object(const RANGE<TV>& box)=0;
//#####################################################################
};
}
#endif
