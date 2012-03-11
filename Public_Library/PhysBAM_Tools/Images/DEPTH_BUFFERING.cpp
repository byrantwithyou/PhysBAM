//#####################################################################
// Copyright 2012, Alexey, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Images/DEPTH_BUFFERING.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> DISPLAY_PRIMITIVE<T>::
DISPLAY_PRIMITIVE(const TV &a,const int style):style(style)
{ 
    type=POINT;
    vertices(0)=a;
    vertices(1)=a;
    vertices(2)=a;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> DISPLAY_PRIMITIVE<T>::
DISPLAY_PRIMITIVE(const TV &a,const TV &b,const int style):style(style)
{ 
    type=SEGMENT;
    vertices(0)=a;
    vertices(1)=b;
    vertices(2)=b;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> DISPLAY_PRIMITIVE<T>::
DISPLAY_PRIMITIVE(const TV &a,const TV &b,const TV &c,const int style):style(style)
{ 
    type=TRIANGLE;
    vertices(0)=a;
    vertices(1)=b;
    vertices(2)=c;
}
//#####################################################################
// Function Initialize_Elements
//#####################################################################
template<class T> void DISPLAY_PRIMITIVE<T>::
Initialize_Elements()
{
    elements.Clean_Memory();
    elements.Append(vertices);
}
//#####################################################################
// Function Initialize_Bounding_Box
//#####################################################################
template<class T> void DISPLAY_PRIMITIVE<T>::
Initialize_Bounding_Box()
{
    bounding_box=RANGE<TV2>::Bounding_Box(
        DB::Project(vertices(0)),
        DB::Project(vertices(1)),
        DB::Project(vertices(2)));
}
//#####################################################################
// Function Cut_By_Primitive
//#####################################################################
template<class T> void DISPLAY_PRIMITIVE<T>::
Cut_By_Primitive(const DISPLAY_PRIMITIVE<T> &p)
{
    if(!(type==POINT||p.type==POINT||RANGE<TV2>::Intersect(bounding_box,p.bounding_box).Empty())){
        switch(p.type){
            case SEGMENT:{
                Cut_By_Plane(DB::Cutting_Plane(p.vertices(0),p.vertices(1)));
            }break;
            case TRIANGLE:{
                Cut_By_Plane(DB::Cutting_Plane(p.vertices(0),p.vertices(1)));
                Cut_By_Plane(DB::Cutting_Plane(p.vertices(1),p.vertices(2)));
                Cut_By_Plane(DB::Cutting_Plane(p.vertices(2),p.vertices(0)));
                Cut_By_Plane(PLANE<T>(p.vertices(0),p.vertices(1),p.vertices(2)));
            }break;
            default:;}}
}
//#####################################################################
// Function Cut_By_Plane
//#####################################################################
template<class T> void DISPLAY_PRIMITIVE<T>::
Cut_By_Plane(const PLANE<T> &p)
{
    
}
//#####################################################################
// Function Cutting_Plane
//#####################################################################
template<class T> PLANE<T> DEPTH_BUFFERING<T>::
Cutting_Plane(const TV &a,const TV &b)
{ 
    return PLANE<T>(Embed(Project(a-b).Rotate_Clockwise_90().Normalized()),a);
}
//#####################################################################
// Function Add_Element
//#####################################################################
template<class T> int DEPTH_BUFFERING<T>::
Add_Element(const TV &a,const int style)
{ 
    primitives.Append(DISPLAY_PRIMITIVE<T>(a,style));
    return primitives.Size()-1;
}
//#####################################################################
// Function Add_Element
//#####################################################################
template<class T> int DEPTH_BUFFERING<T>::
Add_Element(const TV &a,const TV &b,const int style)
{
    primitives.Append(DISPLAY_PRIMITIVE<T>(a,b,style));
    return primitives.Size()-1;
}
//#####################################################################
// Function Add_Element
//#####################################################################
template<class T> int DEPTH_BUFFERING<T>::
Add_Element(const TV &a,const TV &b,const TV &c,const int style)
{
    primitives.Append(DISPLAY_PRIMITIVE<T>(a,b,c,style));
    return primitives.Size()-1;
}
//#####################################################################
// Function Process_Primitives
//#####################################################################
template<class T> void DEPTH_BUFFERING<T>::
Process_Primitives()
{
    for(int i=0;i<primitives.m;i++) {
        primitives(i).Initialize_Elements();
        primitives(i).Initialize_Bounding_Box();
    };
    for(int i=0;i<primitives.m;i++) for(int j=0;j<primitives.m;j++) if(i!=j) primitives(i).Cut_By_Primitive(primitives(j));
}
//#####################################################################
template class DISPLAY_PRIMITIVE<float>;
template class DEPTH_BUFFERING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DISPLAY_PRIMITIVE<double>;
template class DEPTH_BUFFERING<double>;
#endif
