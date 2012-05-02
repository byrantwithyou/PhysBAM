//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEX_FILE
//#####################################################################
#ifndef __TEX_FILE__
#define __TEX_FILE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Images/VECTOR_IMAGE.h>
#include <iostream>
#include <string>
namespace PhysBAM{

template<class T>
class TEX_FILE:public VECTOR_IMAGE<T>
{
public:
    using VECTOR_IMAGE<T>::stream;using VECTOR_IMAGE<T>::output_box;using VECTOR_IMAGE<T>::bounding_box;using VECTOR_IMAGE<T>::cur_format;
    typedef VECTOR<T,2> TV;
    int unit_offset,frame_offset;
    VECTOR<T,3> effective_line_color,effective_fill_color;

    TEX_FILE(const std::string& filename,const RANGE<TV>& box=RANGE<TV>(TV(),TV(500,500)));
    ~TEX_FILE();

    void Emit(const TV &pt);
    bool Emit_Options(bool line,bool fill);

    void Update_Effective_Formatting();

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
