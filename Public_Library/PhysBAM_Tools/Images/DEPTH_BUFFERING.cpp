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
    vertices(0)=a;
    type=POINT;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> PRIMITIVE<T>::
PRIMITIVE(const TV &a,const TV &b,const int style):style(style)
{ 
    vertices(0)=a;
    vertices(1)=b;
    type=SEGMENT;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> PRIMITIVE<T>::
PRIMITIVE(const TV &a,const TV &b,const TV &c,const int style):style(style)
{ 
    vertices(0)=a;
    vertices(1)=b;
    vertices(2)=c;
    type=TRIANGLE;
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
template class PRIMITIVE<float>;
template class DEPTH_BUFFERING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PRIMITIVE<double>;
template class DEPTH_BUFFERING<double>;
#endif
