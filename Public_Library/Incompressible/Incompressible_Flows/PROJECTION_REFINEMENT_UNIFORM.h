//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_REFINEMENT_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_REFINEMENT_UNIFORM__
#define __PROJECTION_REFINEMENT_UNIFORM__

#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <Dynamics/Incompressible_Flows/FAST_PROJECTION_DYNAMICS_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class PROJECTION_REFINEMENT_UNIFORM:public PROJECTION_DYNAMICS_UNIFORM<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
public:
    typedef PROJECTION_DYNAMICS_UNIFORM<TV> BASE;
    using BASE::poisson;using BASE::elliptic_solver;using BASE::p;
    
    THREAD_QUEUE* thread_queue;
public:
    MPI_UNIFORM_GRID<TV> *coarse_mpi_grid,*fine_mpi_grid;
    FAST_PROJECTION_DYNAMICS_UNIFORM<TV> fast_local_projection;
    T_FACE_ARRAYS_SCALAR coarse_face_velocities,coarse_face_velocities_save,face_velocities_save,local_face_velocities;
    T_FACE_ARRAYS_SCALAR &beta_face;
    T_FACE_ARRAYS_BOOL fine_psi_N;
    GRID<TV> coarse_grid,fine_grid,local_grid;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary,solid_wall;
    T scale_face_inverse,alpha;
public:
    int coarse_scale;

    PROJECTION_REFINEMENT_UNIFORM(const GRID<TV>& mac_grid,const int scale,const T alpha=1,const bool flame_input=false,const bool multiphase=false,const bool use_variable_beta=false,const bool use_poisson=false);
    PROJECTION_REFINEMENT_UNIFORM(const GRID<TV>& mac_grid,LEVELSET<TV>& levelset_input,const int scale,const T alpha=1);
    virtual ~PROJECTION_REFINEMENT_UNIFORM();

//#####################################################################
    virtual void Initialize_Grid(const GRID<TV>& mac_grid);
    void Average_Velocities_From_Fine_To_Coarse(ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& fine_face_velocities);
    void Set_Beta_Face_For_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& coarse_face_velocities);
    virtual void Set_Coarse_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& coarse_face_velocities);
    virtual void Map_Fine_To_Coarse(T_FACE_ARRAYS_SCALAR& coarse_face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities);
    void Map_Coarse_To_Fine(const T_FACE_ARRAYS_SCALAR& coarse_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities);
    void Map_Fine_To_Local_Boundary_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities,TV_INT coarse_cell_index);
    void Map_Fine_To_Local_Interior_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities,TV_INT coarse_cell_index,bool zero_out);
    bool Map_Fine_To_Local_Boundaries_For_Cell(GRID<TV>& local_mac_grid,ARRAY<bool,FACE_INDEX<TV::dimension> >& local_psi_N,TV_INT cell_index);
    void Map_Local_To_Fine_Interior_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities,TV_INT cell_index);
    virtual void Local_Projection_PCG(T_FACE_ARRAYS_SCALAR& fine_face_velocities,GRID<TV>& local_grid,T_FACE_ARRAYS_SCALAR& local_face_velocities,FAST_PROJECTION_DYNAMICS_UNIFORM<TV>& local_projection,const T dt,const T time,TV_INT cell_index);
    void Threaded_Local_Projection_PCG(T_FACE_ARRAYS_SCALAR& fine_face_velocities,const T dt,const T time);
    void Local_Projection_PCG(T_FACE_ARRAYS_SCALAR& fine_face_velocities,const T dt,const T time);
    virtual void Map_Coarse_To_Fine(const T_FACE_ARRAYS_SCALAR& coarse_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    virtual void Make_Divergence_Free(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
//#####################################################################
};
}
#endif
