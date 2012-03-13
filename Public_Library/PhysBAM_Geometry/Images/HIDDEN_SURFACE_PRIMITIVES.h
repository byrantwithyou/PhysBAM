//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIDDEN_SURFACE_PRIMITIVES
//#####################################################################
#ifndef __HIDDEN_SURFACE_PRIMITIVES__
#define __HIDDEN_SURFACE_PRIMITIVES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>

namespace PhysBAM{

template<class T>
class SURFACE_PRIMITIVE
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;

public:

    int num_vertices;
    VECTOR<TV,3> vertices;
    int parent;
    RANGE<TV2> bounding_box;

    SURFACE_PRIMITIVE() {num_vertices=0;}
    SURFACE_PRIMITIVE(const TV &a,const int pa);
    SURFACE_PRIMITIVE(const TV &a,const TV &b,const int pa);
    SURFACE_PRIMITIVE(const TV &a,const TV &b,const TV &c,const int pa);
};

template<class T>
class HIDDEN_SURFACE_PRIMITIVES
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;

public:

    ARRAY<SURFACE_PRIMITIVE<T> > primitives;
    ARRAY<int> order;

    void Initialize(DIRECTED_GRAPH<>& dg);
    void Divide_Primitive(int divide,int cutter,ARRAY<int>& inside,ARRAY<int>& outside){}
    bool Test_Edge(int a,int b){return true;} // 0=no edge, 1=a->b
    void Emit_Node(int a) {order.Append(a);}
    bool Projections_Intersect(int a,int b);
};
}
#endif
