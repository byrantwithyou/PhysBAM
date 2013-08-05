//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM__
#define __PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM__

#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/PROJECTION_REFINEMENT_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM:public PROJECTION_REFINEMENT_UNIFORM<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;typedef FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<TV> T_FACE_LOOKUP_FIRE_MULTIPHASE;
public:
    typedef PROJECTION_REFINEMENT_UNIFORM<TV> BASE;
    using BASE::fine_mpi_grid;using BASE::collidable_solver;using BASE::elliptic_solver;using BASE::p;using BASE::coarse_mpi_grid;using BASE::solid_wall;using BASE::fine_psi_N;using BASE::poisson;
    using BASE::coarse_grid;using BASE::local_grid;using BASE::fast_local_projection;using BASE::fine_grid;using BASE::coarse_scale;using BASE::domain_boundary;using BASE::face_velocities_save;
    using BASE::Set_Beta_Face_For_Boundary_Conditions;using BASE::Map_Fine_To_Local_Boundary_For_Cell;using BASE::Map_Fine_To_Local_Interior_For_Cell;
    using BASE::Map_Fine_To_Local_Boundaries_For_Cell;using BASE::Map_Local_To_Fine_Interior_For_Cell;
    
    BOUNDARY<TV,T> *boundary,*phi_boundary;
public:
    LINEAR_INTERPOLATION_UNIFORM<TV,T> phi_interpolation;
    PROJECTION_DYNAMICS_UNIFORM<TV> levelset_projection;
    LEVELSET<TV> &levelset,coarse_levelset;
    ARRAY<T,TV_INT> phi_ghost,coarse_phi,local_phi;
    bool surface_solve;
    int buffer;
public:

    PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM(const GRID<TV>& mac_grid,LEVELSET<TV>& levelset_input,const int scale,const T alpha=1,const bool use_surface_solve=true,const bool flame_input=false,const bool multiphase=false,const bool use_variable_beta=false,const bool use_poisson=false);
    virtual ~PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM();

//#####################################################################
    virtual void Initialize_Grid(const GRID<TV>& mac_grid);
    virtual void Set_Coarse_Boundary_Conditions(ARRAY<T,FACE_INDEX<TV::m> >& coarse_face_velocities);
    bool Set_Local_Boundary_Conditions(GRID<TV>& local_grid,PROJECTION_UNIFORM<TV>& local_projection,TV_INT coarse_index);
    void Set_Local_Phi_From_Fine_Phi(GRID<TV>& local_mac_grid,ARRAY<T,TV_INT>& local_phi,const ARRAY<T,TV_INT>& fine_phi,TV_INT cell_index);
    virtual void Local_Projection_PCG(ARRAY<T,FACE_INDEX<TV::m> >& fine_face_velocities,GRID<TV>& local_grid,ARRAY<T,FACE_INDEX<TV::m> >& local_face_velocities,FAST_PROJECTION_DYNAMICS_UNIFORM<TV>& local_projection,const T dt,const T time,TV_INT cell_index);
    bool Contains_Inside(TV_INT cell_index,const ARRAY<T,TV_INT>& levelset_phi,int buffer);
    bool Contains_Outside(TV_INT cell_index,const ARRAY<T,TV_INT>& levelset_phi,int buffer);
    void Set_Coarse_Phi_From_Fine_Phi(ARRAY<T,TV_INT>& coarse_phi,const ARRAY<T,TV_INT>& fine_phi);
    void Set_Levelset_Boundary_Conditions(const GRID<TV>& levelset_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& levelset_velocities,const ARRAY<T,TV_INT>& levelset_phi,const T time);
    void Map_Fine_To_Levelset_For_Constraints(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    virtual void Map_Fine_To_Coarse(ARRAY<T,FACE_INDEX<TV::m> >& coarse_face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    virtual void Map_Coarse_To_Fine(const ARRAY<T,FACE_INDEX<TV::m> >& coarse_face_velocities,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
//#####################################################################
};
}
#endif
