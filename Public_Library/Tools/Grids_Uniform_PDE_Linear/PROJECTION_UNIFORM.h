//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Michael Lentine, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_UNIFORM__
#define __PROJECTION_UNIFORM__

#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Grids_PDE_Linear/PROJECTION.h>
#include <Tools/Grids_Uniform/SIDED_FACE_INDEX.h>
#include <Tools/Grids_Uniform_PDE_Linear/POISSON_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class PROJECTION_UNIFORM:public PROJECTION<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
public:
    typedef PROJECTION<T> BASE;
    using PROJECTION<T>::use_non_zero_divergence;
    
    GRID<TV> p_grid; // p_grid is a cell centered MAC grid
    ARRAY<T,TV_INT> p;
    ARRAY<T,TV_INT> p_save_for_projection;
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_save_for_projection;
    LAPLACE_UNIFORM<TV>* elliptic_solver;
    LAPLACE_UNIFORM<TV>* laplace; 
    POISSON_UNIFORM<TV>* poisson;     
    ARRAY<T,TV_INT> divergence; // use this to set up a non-zero divergence
    bool use_divergence_multiplier;
    ARRAY<T,TV_INT> divergence_multiplier;
    THREAD_QUEUE* thread_queue;

    PROJECTION_UNIFORM(const GRID<TV>& mac_grid,const bool use_variable_beta=false,const bool use_poisson=false,THREAD_QUEUE* thread_queue=0);
protected:
    PROJECTION_UNIFORM(THREAD_QUEUE* thread_queue_input=0);
public:
    virtual ~PROJECTION_UNIFORM();

    void Use_Non_Zero_Divergence(const bool use_non_zero_divergence_input=true)
    {use_non_zero_divergence=use_non_zero_divergence_input;
    if(use_non_zero_divergence) divergence.Resize(p_grid.Domain_Indices());else divergence.Clean_Memory();}

    void Use_Divergence_Multiplier(const bool use_divergence_multiplier_input=true)
    {use_divergence_multiplier=use_divergence_multiplier_input;
    if(use_divergence_multiplier) divergence_multiplier.Resize(p_grid.Domain_Indices(3));else divergence_multiplier.Clean_Memory();}

//#####################################################################
    virtual void Initialize_Grid(const GRID<TV>& mac_grid);
    virtual void Make_Divergence_Free(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    virtual void Calculate_Kinetic_Energy_Error(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<TV,TV_INT>& kinetic_energy_error);
    void Zero_Out_Neumann_Pocket_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    virtual void Apply_Pressure(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,bool scale_by_dt=false);
    void Enforce_Velocity_Compatibility(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void Set_Up_For_Projection(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void Restore_After_Projection(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void Exchange_Pressures_For_Projection();
    void Compute_Divergence(const T_FACE_LOOKUP& face_lookup,LAPLACE_UNIFORM<TV>* solver);
    void Compute_Divergence_Threaded(RANGE<TV_INT>& domain,const T_FACE_LOOKUP& face_lookup,LAPLACE_UNIFORM<TV>* solver);
//#####################################################################
};
}
#endif
