//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Math_Tools/sqr.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Utilities/DEBUG_CAST.h>
#include <iomanip>
#include "MPM_POISSON_SYSTEM.h"
#include "MPM_POISSON_VECTOR.h"
#include "MPM_PROJECTION.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_POISSON_SYSTEM<TV>::
MPM_POISSON_SYSTEM(MPM_PROJECTION<TV>& proj)
    :BASE(false,true),proj(proj)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_POISSON_SYSTEM<TV>::
~MPM_POISSON_SYSTEM()
{}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MPM_POISSON_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    ARRAY<T,TV_INT>& rr=debug_cast<MPM_POISSON_VECTOR<TV>&>(result).v;
    const ARRAY<T,TV_INT>& xx=debug_cast<const MPM_POISSON_VECTOR<TV>&>(x).v;
    T one_over_h_square=sqr((T)1/proj.mac_grid.dX.Min());
    T one_over_dt=(T)1.0/proj.sim.dt;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+proj.mac_grid.counts));it.Valid();it.Next()){
        if(proj.cell_dirichlet(it.index) || proj.cell_neumann(it.index))
            rr(it.index)=(T)0;
        else{ // cell is fluid
            rr(it.index)=(T)0;
            for(int d=0;d<TV::m;d++){
                TV_INT left_cell_index=it.index;left_cell_index(d)--;
                TV_INT right_cell_index=it.index;right_cell_index(d)++;
                FACE_INDEX<TV::m> left_face_index(d,proj.mac_grid.First_Face_Index_In_Cell(d,it.index));
                if(proj.face_masses(left_face_index)>proj.sim.min_mass){
                    if(proj.cell_dirichlet(left_cell_index)) rr(it.index)+=xx(it.index)/proj.face_masses(left_face_index);
                    else if(proj.cell_neumann(left_cell_index)) rr(it.index)+=T(0);
                    else rr(it.index)+=(xx(it.index)-xx(left_cell_index))/proj.face_masses(left_face_index);
                }
                FACE_INDEX<TV::m> right_face_index(d,proj.mac_grid.Second_Face_Index_In_Cell(d,it.index));
                if(proj.face_masses(right_face_index)>proj.sim.min_mass){
                    if(proj.cell_dirichlet(right_cell_index)) rr(it.index)+=xx(it.index)/proj.face_masses(right_face_index);
                    else if(proj.cell_neumann(right_cell_index)) rr(it.index)+=T(0);
                    else rr(it.index)+=(xx(it.index)-xx(right_cell_index))/proj.face_masses(right_face_index);
                }}
            rr(it.index)*=proj.sim.dt*one_over_h_square;
            rr(it.index)+=one_over_dt*proj.one_over_lambda_J(it.index)*xx(it.index);}}
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void MPM_POISSON_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    ARRAY<T,TV_INT>& xx=debug_cast<MPM_POISSON_VECTOR<TV>&>(x).v;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),proj.mac_grid.counts));it.Valid();it.Next())
        if(proj.cell_dirichlet(it.index) || proj.cell_neumann(it.index))
            xx(it.index)=T(0);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MPM_POISSON_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const ARRAY<T,TV_INT>& xx=debug_cast<const MPM_POISSON_VECTOR<TV>&>(x).v;
    const ARRAY<T,TV_INT>& yy=debug_cast<const MPM_POISSON_VECTOR<TV>&>(y).v;
    double r=0;
    for(int i=0;i<xx.array.m;i++) r+=xx.array(i)*yy.array(i);
    return r;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MPM_POISSON_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    const ARRAY<T,TV_INT>& xx=debug_cast<const MPM_POISSON_VECTOR<TV>&>(x).v;
//    LOG::cout<<"resudual: ";
//    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),proj.mac_grid.counts));it.Valid();it.Next())
//        LOG::cout<<xx(it.index)<<" ";
//    LOG::cout<<std::endl;
    
    
    return xx.array.Maximum_Magnitude();
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void MPM_POISSON_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void MPM_POISSON_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void MPM_POISSON_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
}
//#####################################################################
template class MPM_POISSON_SYSTEM<VECTOR<float,2> >;
template class MPM_POISSON_SYSTEM<VECTOR<float,3> >;
template class MPM_POISSON_SYSTEM<VECTOR<double,2> >;
template class MPM_POISSON_SYSTEM<VECTOR<double,3> >;
}
