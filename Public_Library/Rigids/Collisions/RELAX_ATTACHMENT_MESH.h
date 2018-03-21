//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RELAX_ATTACHMENT_MESH
//#####################################################################
#ifndef __RELAX_ATTACHMENT_MESH__
#define __RELAX_ATTACHMENT_MESH__

#include <Core/Matrices/MATRIX.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV>
class RELAX_ATTACHMENT_MESH
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
public:
    struct DIFF_ENTRY
    {
        int e; // element of relaxed attachment point
        TV Y; // relaxed attachment point
        VECTOR<MATRIX<T,TV::m>,5> dYdI; // Dependence on X, Z, A, B, C
    };

    ARRAY<DIFF_ENTRY> diff_entry;
    TV w; // barycentric coords of relaxed attachment point
    VECTOR<MATRIX<T,TV::m>,5> dwdI; // Dependence on X, Z, A, B, C
    bool active;
    int e; // Final attachment point element
    TV Y; // Final attachment point
    
    void Relax(int e0,const TV& w0,const TV& Z,const T_SURFACE& surface,int ex_pt,T friction);
};
}
#endif
