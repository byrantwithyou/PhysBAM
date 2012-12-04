//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Images/VECTOR_IMAGE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> VECTOR_IMAGE<T>::FORMATTING::
FORMATTING()
    :line_width(1),point_radius(3),line_opacity(1),fill_opacity(1),line_style(1),fill_style(0),arrow_style(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> VECTOR_IMAGE<T>::
VECTOR_IMAGE(const std::string& filename,const RANGE<TV>& box)
    :stream(filename.c_str()),bounding_box(RANGE<TV>::Empty_Box()),output_box(box),fixed_bounding_box(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> VECTOR_IMAGE<T>::
~VECTOR_IMAGE()
{
}
//#####################################################################
// Function Bound
//#####################################################################
template<class T> void VECTOR_IMAGE<T>::
Bound(const TV& pt)
{
    if(!fixed_bounding_box) bounding_box.Enlarge_To_Include_Point(pt);
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void VECTOR_IMAGE<T>::
Draw_Object(const TV &a,const TV &b)
{
    Bound(a);
    Bound(b);
    Emit_Object(a,b);
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void VECTOR_IMAGE<T>::
Draw_Object(const TV &a,const TV &b,const TV &c)
{
    Bound(a);
    Bound(b);
    Bound(c);
    Emit_Object(a,b,c);
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void VECTOR_IMAGE<T>::
Draw_Object(ARRAY_VIEW<TV> pts)
{
    for(int i=0;i<pts.Size();i++) Bound(pts(i));
    Emit_Object(pts);
}
template<class T> void VECTOR_IMAGE<T>::
Draw_Object(ARRAY_VIEW<TV> outside,ARRAY_VIEW<ARRAY_VIEW<TV> > holes)
{
    for(int i=0;i<outside.Size();i++) Bound(outside(i));
    Emit_Object(outside,holes);
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void VECTOR_IMAGE<T>::
Draw_Object(const TV &pt,T radius)
{
    Bound(pt-radius);
    Bound(pt+radius);
    Emit_Object(pt,radius);
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void VECTOR_IMAGE<T>::
Draw_Object(const RANGE<TV>& box)
{
    Bound(box.min_corner);
    Bound(box.max_corner);
    Emit_Object(box);
}
namespace PhysBAM{
template class VECTOR_IMAGE<float>;
template class VECTOR_IMAGE<double>;
}
