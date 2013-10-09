//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EPS_FILE
//#####################################################################
#ifndef __EPS_FILE__
#define __EPS_FILE__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Images/VECTOR_IMAGE.h>
#include <iostream>
#include <string>
namespace PhysBAM{

template<class T>
class EPS_FILE:public VECTOR_IMAGE<T>
{
public:
    using VECTOR_IMAGE<T>::stream;using VECTOR_IMAGE<T>::output_box;using VECTOR_IMAGE<T>::bounding_box;using VECTOR_IMAGE<T>::cur_format;
    using VECTOR_IMAGE<T>::Bound;
    typedef VECTOR<T,2> TV;typedef VECTOR<T,3> TV3;
    int head_offset;

    TV3 effective_color;
    T effective_line_width,effective_point_radius,effective_line_opacity,effective_fill_opacity;

    EPS_FILE(const std::string& filename,const RANGE<TV>& box=RANGE<TV>(TV(),TV(500,500)));
    ~EPS_FILE();

    void Emit(const std::string& str);
    void Emit(const TV &pt);
    void Emit(const TV3 &col);

    void Update_Effective_Formatting();
    void Mt(const TV &pt);
    void Lt(const TV &pt);
    void Stroke();
    void Fill();

    void Emit_Gouraud_Triangle(const TV &a,const TV &b,const TV &c,const TV3 &col_a,const TV3 &col_b,const TV3 &col_c);
    void Use_Normal_Cap();
    void Use_Round_Cap();
    void Use_Extended_Cap();
    void Use_Miter_Join();
    void Use_Round_Join();
    void Use_Bevel_Join();

protected:
    virtual void Emit_Object(const TV &a,const TV &b); // line
    virtual void Emit_Object(const TV &a,const TV &b,const TV &c); // triangle
    virtual void Emit_Object(ARRAY_VIEW<TV> pts); // polygon
    virtual void Emit_Object(ARRAY_VIEW<TV> pts,ARRAY_VIEW<ARRAY_VIEW<TV> > holes); // polygon with holes
    virtual void Emit_Object(const TV &pt,T radius); // circle
    virtual void Emit_Object(const RANGE<TV>& box);
    void Compute_Transform(T& scale,TV& shift);
//#####################################################################
};
}
#endif
