//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_COLLIDABLE_UNIFORM
//#####################################################################
#ifndef __LAPLACE_COLLIDABLE_UNIFORM__
#define __LAPLACE_COLLIDABLE_UNIFORM__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Poisson/LAPLACE_UNIFORM.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/LAPLACE_COLLIDABLE.h>
namespace PhysBAM{

template<class T> class PCG_SPARSE;
template<class T> class SPARSE_MATRIX_FLAT_MXN;

template<class TV>
class LAPLACE_COLLIDABLE_UNIFORM:public LAPLACE_UNIFORM<TV>,public LAPLACE_COLLIDABLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
    typedef INTERPOLATION_UNIFORM<TV,T> T_INTERPOLATION_SCALAR;
public:
    typedef LAPLACE_UNIFORM<TV> BASE;
    typedef LAPLACE_COLLIDABLE<TV> COLLIDABLE_BASE;
    typedef TV VECTOR_T;
    typedef GRID<TV> GRID_T;

    using BASE::grid;using BASE::psi_N;using BASE::psi_D;using BASE::u;
    using LAPLACE<T>::tolerance;using LAPLACE<T>::number_of_regions;using LAPLACE<T>::solve_neumann_regions;
    using COLLIDABLE_BASE::second_order_cut_cell_method;using COLLIDABLE_BASE::second_order_cut_cell_threshold;
    using COLLIDABLE_BASE::levelset;using COLLIDABLE_BASE::u_interface;using BASE::filled_region_colors;

    //LEVELSET<TV>* levelset; // used in second order accurate cut cell method
    //ARRAY<T,FACE_INDEX<TV::m> > u_interface; // interface boundary condition - 2nd order method
protected:
    ARRAY<T,TV_INT> phi_default;
    LEVELSET<TV>* levelset_default;
public:

    LAPLACE_COLLIDABLE_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input);
    LAPLACE_COLLIDABLE_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,LEVELSET<TV>& cell_centered_levelset,const bool initialize_grid,const bool multiphase_input,
        const bool enforce_compatibility_input);
    virtual ~LAPLACE_COLLIDABLE_UNIFORM();

    void Use_Internal_Level_Set()
    {levelset=levelset_default;phi_default.Resize(grid.Domain_Indices(1));}

    void Use_External_Level_Set(LEVELSET<TV>& cell_centered_levelset) override
    {levelset=&cell_centered_levelset;phi_default.Clean_Memory();}

//#####################################################################
    virtual void Initialize_Grid(const GRID<TV>& mac_grid_input) override;
    virtual void Set_Up_Second_Order_Cut_Cell_Method(const bool use_second_order_cut_cell_method_input=true) override;
    void Find_A(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index) override;
private:
    void Apply_Second_Order_Cut_Cell_Method(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index);
//#####################################################################
};
}
#endif
