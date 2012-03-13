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

    int num_veritces;
    VECTOR<TV,3> vertices;
    int parent;
    RANGE<TV2> bounding_box;

    SURFACE_PRIMITIVE() {type=EMPTY;}
    SURFACE_PRIMITIVE(const TV &a,const int style);
    SURFACE_PRIMITIVE(const TV &a,const TV &b,const int style);
    SURFACE_PRIMITIVE(const TV &a,const TV &b,const TV &c,const int style);

    void Initialize_Bounding_Box();
};

template<class T>
class HIDDEN_SURFACE_PRIMITIVES
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,2> TV2;

public:
    void Initialize(DIRECTED_GRAPH<>& dg);
    void Divide_Primitive(int divide,int cutter,ARRAY<int>& inside,ARRAY<int>& outside);
    bool Test_Edge(int a,int b); // 0=no edge, 1=a->b
    int Emit_Node(int a);
};
}
#endif
