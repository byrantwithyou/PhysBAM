//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Images/DEPTH_BUFFERING.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> PRIMITIVE<T>::
PRIMITIVE()
{ 
    type=EMPTY;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> PRIMITIVE<T>::
PRIMITIVE(const TV &a,const int style):style(style)
{ 
    type=POINT;
    vertices(0)=a;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> PRIMITIVE<T>::
PRIMITIVE(const TV &a,const TV &b,const int style):style(style)
{ 
    type=SEGMENT;
    vertices(0)=a;
    vertices(1)=b;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> PRIMITIVE<T>::
PRIMITIVE(const TV &a,const TV &b,const TV &c,const int style):style(style)
{ 
    type=TRIANGLE;
    vertices(0)=a;
    vertices(1)=b;
    vertices(2)=c;
}
//#####################################################################
// Function Add_Element
//#####################################################################
template<class T> int DEPTH_BUFFERING<T>::
Add_Element(const TV &a,const int style)
{ 
    primitives.Append(PRIMITIVE<T>(a,style));
    return primitives.Size()-1;
}
template<class T> int DEPTH_BUFFERING<T>::
//#####################################################################
// Function Add_Element
//#####################################################################
Add_Element(const TV &a,const TV &b,const int style)
{
    primitives.Append(PRIMITIVE<T>(a,b,style));
    return primitives.Size()-1;
}
//#####################################################################
// Function Add_Element
//#####################################################################
template<class T> int DEPTH_BUFFERING<T>::
Add_Element(const TV &a,const TV &b,const TV &c,const int style)
{
    primitives.Append(PRIMITIVE<T>(a,b,c,style));
    return primitives.Size()-1;
}
//#####################################################################
// Function Process_Primitives
//#####################################################################
template<class T> void DEPTH_BUFFERING<T>::
Process_Primitives(const TV &a,const TV &b,const TV &c,const int style)
{
    for(int i=0;i<primitives.m;i++){
        primitives(i).Remove_All();
        primitives(i).Clean_Memory();}
    
    
}
//#####################################################################
template class PRIMITIVE<float>;
template class DEPTH_BUFFERING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PRIMITIVE<double>;
template class DEPTH_BUFFERING<double>;
#endif
