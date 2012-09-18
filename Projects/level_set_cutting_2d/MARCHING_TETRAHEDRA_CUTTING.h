//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHING_TETRAHEDRA_CUTTING
//#####################################################################
#ifndef __MARCHING_TETRAHEDRA_CUTTING__
#define __MARCHING_TETRAHEDRA_CUTTING__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<int d>
struct MARCHING_TETRAHEDRA_CUTTING_CASE
{
    enum WORKAROUND {max_side_elements=d};
    unsigned char comp[4];
    unsigned short elements[2][max_side_elements];
};

template<class TV>
class MARCHING_TETRAHEDRA_CUTTING
{
public:
    enum WORKAROUND {max_side_elements=TV::m};
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<int,TV::m+1> E;typedef VECTOR<int,2> S;typedef VECTOR<T,2> P;

    MARCHING_TETRAHEDRA_CUTTING() {}
    ~MARCHING_TETRAHEDRA_CUTTING() {}

    static const ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<TV::m> >& Case_Table();
    static void Initialize_Case_Table(ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<2> >& table);
    static void Initialize_Case_Table(ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<3> >& table);
    static void Initialize_Neighbor_Cases(ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<TV::m> >& table, int c);
    static void Query_Case(ARRAY<E>& parents,ARRAY<E>& children,ARRAY<E>& split_parents,
        const ARRAY<T>& phi,ARRAY<PAIR<S,T> >& weights);
//#####################################################################
};
}
#endif
