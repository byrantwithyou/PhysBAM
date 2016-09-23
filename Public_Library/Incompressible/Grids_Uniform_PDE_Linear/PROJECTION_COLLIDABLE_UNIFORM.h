//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_COLLIDABLE_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_COLLIDABLE_UNIFORM__
#define __PROJECTION_COLLIDABLE_UNIFORM__

#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class PROJECTION_COLLIDABLE_UNIFORM:public PROJECTION_UNIFORM<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
public:
    typedef PROJECTION_UNIFORM<TV> BASE;
    using BASE::p_grid;using BASE::p;using BASE::elliptic_solver;using BASE::poisson;using BASE::laplace;using BASE::Zero_Out_Neumann_Pocket_Velocities;
    
    LAPLACE_COLLIDABLE<TV>* collidable_solver;
    LAPLACE_COLLIDABLE_UNIFORM<TV>* laplace_collidable; 
    POISSON_COLLIDABLE_UNIFORM<TV>* poisson_collidable;

    PROJECTION_COLLIDABLE_UNIFORM(const GRID<TV>& mac_grid,const bool multiphase,const bool use_poisson,const bool use_variable_beta,THREAD_QUEUE* thread_queue=0);
    PROJECTION_COLLIDABLE_UNIFORM(const GRID<TV>& mac_grid,LEVELSET<TV>& levelset_input);
    virtual ~PROJECTION_COLLIDABLE_UNIFORM();

//#####################################################################
    virtual void Initialize_Grid(const GRID<TV>& mac_grid) override;
    virtual void Apply_Pressure(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,bool scale_by_dt=false) override;
//#####################################################################
};
}
#endif
