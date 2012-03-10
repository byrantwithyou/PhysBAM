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
struct PRIMITIVE
{
    enum TYPE{EMPTY=0,POINT=1,SEGMENT,TRIANGLE};
    typedef VECTOR<T,3> TV;
    
    VECTOR<TV,3> vertices;
    TYPE type;
    
    PRIMITIVE() {type=EMPTY;}
    PRIMITIVE(const TV &a) {vertices(0)=a; type=POINT;}
    PRIMITIVE(const TV &a,const TV &b) {vertices(0)=a; vertices(1)=b; type=SEGMENT;}
    PRIMITIVE(const TV &a,const TV &b,const TV &c) {vertices(0)=a; vertices(1)=b; vertices(3)=c; type=TRIANGLE;}
};

template<class T>
class DEPTH_BUFFERING
{
    typedef VECTOR<T,3> TV;

    ARRAY<PRIMITIVE<T> > primitives;

public:

    DEPTH_BUFFERING(){}
    ~DEPTH_BUFFERING(){}
    
    int Add_Element(const TV &a) {primitives.Append(PRIMITIVE<T>(a)); return primitives.Size()-1;}
    int Add_Element(const TV &a,const TV &b) {primitives.Append(PRIMITIVE<T>(a,b)); return primitives.Size()-1;}
    int Add_Element(const TV &a,const TV &b,const TV &c) {primitives.Append(PRIMITIVE<T>(a,b,c)); return primitives.Size()-1;}

//#####################################################################
};
}
#endif
#endif
