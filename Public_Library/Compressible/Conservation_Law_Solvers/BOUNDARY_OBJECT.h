//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_OBJECT  
//#####################################################################
#ifndef __BOUNDARY_OBJECT__
#define __BOUNDARY_OBJECT__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_CALLBACKS.h>
namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;

template<class TV,class TV_DIMENSION>
class BOUNDARY_OBJECT
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename ARRAY<T,TV_INT>::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_DIMENSION_SCALAR::ELEMENT T_ARRAYS_ELEMENT;
public:
    BOUNDARY_OBJECT();
    virtual ~BOUNDARY_OBJECT();

//#####################################################################
    void Get_State_At_Location(const GRID<VECTOR<T,1> >& grid_1d,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_1d,const T location,const VECTOR<int,2>& region_boundaries,TV_DIMENSION& u_1d);
    void Fill_Ghost_Cells_Neumann(const GRID<VECTOR<T,1> >& grid_1d,ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_1d,const T_FACE_ARRAYS_SCALAR& face_velocities,const VECTOR<int,TV::m-1>& node_lower_dimension,const int axis,
        const int ghost_cells,const bool use_exact_neumann_face_location,const VECTOR<int,2>& domain,const VECTOR<int,2>& region_boundaries,const VECTOR<bool,2>& psi_N,
        CONSERVATION_CALLBACKS<T>* callbacks);
    virtual void Apply_Neumann_Boundary_Condition(TV_DIMENSION& u_1d,const T neumann_face_velocity,const int axis)=0;
    virtual void Apply_Neumann_Boundary_Condition(T_ARRAYS_ELEMENT& u_1d,const TV& normal,const T object_velocity_normal_component)=0;
//#####################################################################
};
}
#endif
