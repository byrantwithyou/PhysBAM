//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEPTH_BUFFERING
//#####################################################################
#ifndef __DEPTH_BUFFERING__
#define __DEPTH_BUFFERING__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>

namespace PhysBAM{

template<class T> class DEPTH_BUFFERING;

template<class T>
class DISPLAY_PRIMITIVE
{
public:

    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;
    typedef DEPTH_BUFFERING<T> DB;

    enum TYPE{EMPTY=0,POINT,SEGMENT,TRIANGLE};
    
    VECTOR<TV,3> vertices;
    TYPE type;
    int style;

    RANGE<TV2> bounding_box;
    ARRAY<VECTOR<TV, 3> > elements;

    DISPLAY_PRIMITIVE() {type=EMPTY;}
    DISPLAY_PRIMITIVE(const TV &a,const int style);
    DISPLAY_PRIMITIVE(const TV &a,const TV &b,const int style);
    DISPLAY_PRIMITIVE(const TV &a,const TV &b,const TV &c,const int style);

    void Initialize_Bounding_Box();
    void Initialize_Elements();
    
    void Cut_By_Primitive(const DISPLAY_PRIMITIVE<T> &p);
    void Cut_By_Plane(const PLANE<T> &p);
};

template<class T>
class DEPTH_BUFFERING
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;

    ARRAY<DISPLAY_PRIMITIVE<T> > primitives;

public:

    DEPTH_BUFFERING(){}
    ~DEPTH_BUFFERING(){}

    static VECTOR<T,2> Project(const VECTOR<T,3> &a) {return VECTOR<T,2>(a(0),a(1));}
    static VECTOR<T,3> Embed  (const VECTOR<T,2> &a) {return VECTOR<T,3>(a(0),a(1),0);}
    static PLANE<T> Get_Cutting_Plane(const TV &a,const TV &b);

    int Add_Element(const TV &a,const int style);
    int Add_Element(const TV &a,const TV &b,const int style);
    int Add_Element(const TV &a,const TV &b,const TV &c,const int style);

    void Process_Primitives();
//#####################################################################
};
}
#endif
