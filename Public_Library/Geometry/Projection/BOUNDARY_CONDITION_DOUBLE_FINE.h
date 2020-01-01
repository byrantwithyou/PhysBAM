//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CONDITION_DOUBLE_FINE
//#####################################################################
#ifndef __BOUNDARY_CONDITION_DOUBLE_FINE__
#define __BOUNDARY_CONDITION_DOUBLE_FINE__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV> class LEVELSET;

template<class TV>
class BOUNDARY_CONDITION_DOUBLE_FINE
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;

public:
    enum bc_enum {bc_slip=-3,bc_noslip=-2,bc_free=-1};
    GRID<TV> mac_grid; // original grid
    GRID<TV> grid; // double fine, node grid.
    int ghost;
    ARRAY<char,TV_INT> bc_type;
    HASHTABLE<TV_INT,PAIR<T,int> > bc_p;
    HASHTABLE<FACE_INDEX<TV::m>,TRIPLE<T,int,bool> > bc_u;

    // Sufficient data to solve BC data extending out ghost cells.
    BOUNDARY_CONDITION_DOUBLE_FINE(const GRID<TV>& mac_grid,int ghost);
    BOUNDARY_CONDITION_DOUBLE_FINE(const BOUNDARY_CONDITION_DOUBLE_FINE&) = delete;
    void operator=(const BOUNDARY_CONDITION_DOUBLE_FINE&) = delete;
    ~BOUNDARY_CONDITION_DOUBLE_FINE()=default;

    void Reset(char type);

    void Set(const LEVELSET<TV>& ls,char type,std::function<T(const TV& X)> f=0,bool thin=false,bool invert=false,T contour=0);
    void Set(const T_SURFACE& surface,char type,std::function<T(const TV& X, int e)> f=0,bool thin=false);

    void Set_Pressure_Boundary_Conditions(ARRAY<bool,TV_INT>& psi_D,
        ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& p) const;
    void Set_Viscosity_Boundary_Conditions(ARRAY<bool,TV_INT>& psi_D,
        ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& u,int axis) const;
//#####################################################################
};
}
#endif
