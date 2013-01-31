//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_PROJECTION_DYNAMICS_UNIFORM  
//#####################################################################
#ifndef __FAST_PROJECTION_DYNAMICS_UNIFORM__
#define __FAST_PROJECTION_DYNAMICS_UNIFORM__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID> class DETONATION_SHOCK_DYNAMICS;

template<class T_GRID>
class FAST_PROJECTION_DYNAMICS_UNIFORM:public PROJECTION_DYNAMICS_UNIFORM<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_ARRAYS_SLIP T_FACE_ARRAYS_SLIP_SCALAR;
    typedef UNIFORM_GRID_ITERATOR_FACE<TV> FACE_ITERATOR;typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename FIRE_INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP_FIRE_MULTIPHASE T_FACE_LOOKUP_FIRE_MULTIPHASE;
public:
    typedef PROJECTION_DYNAMICS_UNIFORM<T_GRID> BASE;using BASE::p_grid;using BASE::elliptic_solver;using BASE::Compute_Divergence;using BASE::Apply_Pressure;
    
    SPARSE_MATRIX_FLAT_NXN<T> A;
    ARRAY<T> b;
    ARRAY<int,TV_INT> cell_index_to_matrix_index;
    ARRAY<TV_INT,int> matrix_index_to_cell_index;

    FAST_PROJECTION_DYNAMICS_UNIFORM(const int scale,const bool flame_input=false,const bool multiphase=false,const bool use_variable_beta=false,const bool use_poisson=false);
    FAST_PROJECTION_DYNAMICS_UNIFORM(const int scale,LEVELSET<TV>& levelset_input);
    virtual ~FAST_PROJECTION_DYNAMICS_UNIFORM();

//#####################################################################
    virtual void Initialize_Grid(const T_GRID& mac_grid);
    void Make_Divergence_Free_Fast(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
//#####################################################################
};
}
#endif
