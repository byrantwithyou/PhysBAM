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

#define DB_TOL_PLANE    1e-7
#define DB_TOL_POINT    2e-2
#define DB_TOL_SEGMENT  1e-2
#define DB_SIZE_POINT   2e-2
#define DB_SIZE_SEGMENT 1e-2

namespace PhysBAM{

template<class T> class DEPTH_BUFFERING;

template<class T>
class DISPLAY_PRIMITIVE
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;
    typedef DEPTH_BUFFERING<T> DB;

public:

    enum TYPE{EMPTY=0,POINT,SEGMENT,TRIANGLE};
    
    VECTOR<TV,3> vertices;
    TYPE type;
    int style;
    RANGE<TV2> bounding_box;

    DISPLAY_PRIMITIVE() {type=EMPTY;}
    DISPLAY_PRIMITIVE(const TV &a,const int style);
    DISPLAY_PRIMITIVE(const TV &a,const TV &b,const int style);
    DISPLAY_PRIMITIVE(const TV &a,const TV &b,const TV &c,const int style);
    void Initialize_Bounding_Box();
    VECTOR<T,3> Centroid();
};

template<class T>
class DISPLAY_PRIMITIVE_CUTTING:public DISPLAY_PRIMITIVE<T>
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;
    typedef DEPTH_BUFFERING<T> DB;
    typedef DISPLAY_PRIMITIVE<T> BASE;

    using BASE::EMPTY;using BASE::POINT;using BASE::SEGMENT;using BASE::TRIANGLE;
    using BASE::vertices;using BASE::type;using BASE::style;using BASE::bounding_box;

public:

    DISPLAY_PRIMITIVE_CUTTING():DISPLAY_PRIMITIVE<T>(){}
    DISPLAY_PRIMITIVE_CUTTING(const TV &a,const int style):DISPLAY_PRIMITIVE<T>(a,style){}
    DISPLAY_PRIMITIVE_CUTTING(const TV &a,const TV &b,const int style):DISPLAY_PRIMITIVE<T>(a,b,style){}
    DISPLAY_PRIMITIVE_CUTTING(const TV &a,const TV &b,const TV &c,const int style):DISPLAY_PRIMITIVE<T>(a,b,c,style){}
        
    ARRAY<VECTOR<TV, 3> > elements;
    void Initialize_Elements();
    void Cut_By_Primitive(const DISPLAY_PRIMITIVE_CUTTING<T> &p);
    void Cut_By_Plane(const PLANE<T> &p,T tol=DB_TOL_PLANE);
};

template<class T>
class DISPLAY_PRIMITIVE_ORDERING:public DISPLAY_PRIMITIVE<T>
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;
    typedef DEPTH_BUFFERING<T> DB;
    typedef DISPLAY_PRIMITIVE<T> BASE;

    using BASE::EMPTY;using BASE::POINT;using BASE::SEGMENT;using BASE::TRIANGLE;
    using BASE::vertices;using BASE::type;using BASE::style;using BASE::bounding_box;

public:
    
    int parent;

    DISPLAY_PRIMITIVE_ORDERING():DISPLAY_PRIMITIVE<T>(){}
    DISPLAY_PRIMITIVE_ORDERING(const TV &a,const int style,const int parent):DISPLAY_PRIMITIVE<T>(a,style),parent(parent){}
    DISPLAY_PRIMITIVE_ORDERING(const TV &a,const TV &b,const int style,const int parent):DISPLAY_PRIMITIVE<T>(a,b,style),parent(parent){}
    DISPLAY_PRIMITIVE_ORDERING(const TV &a,const TV &b,const TV &c,const int style,const int parent):DISPLAY_PRIMITIVE<T>(a,b,c,style),parent(parent){}
};

template<class T>
class DEPTH_BUFFERING
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;

public:

    ARRAY<DISPLAY_PRIMITIVE_CUTTING<T> > primitives_cutting;
    ARRAY<DISPLAY_PRIMITIVE_ORDERING<T> > primitives_ordering;

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
