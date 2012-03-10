//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEPTH_BUFFERING
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __DEPTH_BUFFERING__
#define __DEPTH_BUFFERING__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{

template<class T>
class PRIMITIVE
{
public:

    enum TYPE{EMPTY,POINT,SEGMENT,TRIANGLE};
    typedef VECTOR<T,3> TV;
    
    VECTOR<TV,3> vertices;
    TYPE type;
    int style;

    ARRAY<VECTOR<TV, 3> > elements;
    
    PRIMITIVE();
    PRIMITIVE(const TV &a,const int style);
    PRIMITIVE(const TV &a,const TV &b,const int style);
    PRIMITIVE(const TV &a,const TV &b,const TV &c,const int style);
};

template<class T>
class DEPTH_BUFFERING
{
    typedef VECTOR<T,3> TV;

    ARRAY<PRIMITIVE<T> > primitives;

public:

    DEPTH_BUFFERING(){}
    ~DEPTH_BUFFERING(){}
    
    int Add_Element(const TV &a,const int style);
    int Add_Element(const TV &a,const TV &b,const int style);
    int Add_Element(const TV &a,const TV &b,const TV &c,const int style);

//#####################################################################
};
}
#endif
#endif
