//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_CUTTING.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_FLUID_SYSTEM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_FLUID_SYSTEM<TV>::
INTERFACE_FLUID_SYSTEM(const GRID<TV>& grid_input,GRID<TV>& coarse_grid_input,ARRAY<T,TV_INT>& phi_input)
    :BASE(false,false),grid(grid_input),coarse_grid(coarse_grid_input),phi(coarse_grid_input,phi_input,0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_FLUID_SYSTEM<TV>::
~INTERFACE_FLUID_SYSTEM()
{
    for(int i=0;i<TV::m;i++) delete index_map_u[i];
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class TV> void INTERFACE_FLUID_SYSTEM<TV>::
Set_Matrix(const VECTOR<T,2>& mu)
{
    BASIS_STENCIL_UNIFORM<TV,0> p_stencil(grid.dX);
    BASIS_STENCIL_UNIFORM<TV,1> *u_stencil[TV::m],*udx_stencil[TV::m][TV::m];
    p_stencil.Set_Center();
    p_stencil.Set_Constant_Stencil();
    p_stencil.Dice_Stencil();
    for(int i=0;i<TV::m;i++){
        u_stencil[i]=new BASIS_STENCIL_UNIFORM<TV,1>(grid.dX);
        u_stencil[i]->Set_Face(i);
        u_stencil[i]->Set_Multilinear_Stencil();
        u_stencil[i]->Dice_Stencil();
        for(int j=0;j<TV::m;j++){
            udx_stencil[i][j]=new BASIS_STENCIL_UNIFORM<TV,1>(*u_stencil[i]);
            udx_stencil[i][j]->Differentiate(j);
            udx_stencil[i][j]->Dice_Stencil();}}

    MARCHING_CUBES<TV>::Create_Surface(object,coarse_grid,phi.phi);

    index_map_p=new CELL_MAPPING<TV>(grid);
    index_map_p->periodic.Fill(true);
    for(int i=0;i<TV::m;i++){
        index_map_u[i]=new CELL_MAPPING<TV>(grid);
        index_map_u[i]->periodic.Fill(true);}

    SYSTEM_MATRIX_HELPER<T> helper;
    RANGE<TV_INT> boundary_conditions;
    boundary_conditions.min_corner.Fill(BASIS_INTEGRATION_CUTTING<TV,2>::periodic);
    boundary_conditions.max_corner.Fill(BASIS_INTEGRATION_CUTTING<TV,2>::periodic);

    SYSTEM_MATRIX_HELPER<T> helper_uu[TV::m][TV::m],helper_p[TV::m],helper_q[TV::m];
    BASIS_INTEGRATION_CUTTING<TV,2> bic(boundary_conditions,grid,coarse_grid,phi.phi);

    // Diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++){
            bic.Add_Block(helper_uu[i][i],*udx_stencil[i][j],*udx_stencil[i][j],*index_map_u[i],*index_map_u[i],mu*(1+(i==j)));}

    // Off-diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=i+1;j<TV::m;j++)
            bic.Add_Block(helper_uu[i][j],*udx_stencil[i][j],*udx_stencil[j][i],*index_map_u[i],*index_map_u[j],mu);

    // Pressure blocks
    for(int i=0;i<TV::m;i++)
        bic.Add_Block(helper_p[i],*udx_stencil[i][i],p_stencil,*index_map_u[i],*index_map_p,VECTOR<T,2>(1,1));

    // Traction blocks
    for(int i=0;i<TV::m;i++)
        bic.Add_Block(helper_q[i],*u_stencil[i],*index_map_u[i],VECTOR<T,2>(1,1));

    // Rhs traction blocks
    for(int i=0;i<TV::m;i++){
        helper_rhs_q[i]=new SYSTEM_MATRIX_HELPER<T>;
        bic.Add_Block(*helper_rhs_q[i],*u_stencil[i],*index_map_u[i],-(T).5*VECTOR<T,2>(1,-1));}

    // Rhs pressure blocks
    for(int i=0;i<TV::m;i++){
        helper_rhs_p[i]=new SYSTEM_MATRIX_HELPER<T>;
        bic.Add_Block(*helper_rhs_p[i],*u_stencil[i],p_stencil,*index_map_u[i],*index_map_p,VECTOR<T,2>(1,1));}

    bic.Compute();

    index_range_u[0]=INTERVAL<int>(0,index_map_u[0]->next_index);
    for(int i=1;i<TV::m;i++) index_range_u[i]=INTERVAL<int>(index_range_u[i-1].max_corner,index_range_u[i-1].max_corner+index_map_u[i]->next_index);
    index_range_p=INTERVAL<int>(index_range_u[TV::m-1].max_corner,index_range_u[TV::m-1].max_corner+index_map_p->next_index);
    index_range_q[0]=INTERVAL<int>(index_range_p.max_corner,index_range_p.max_corner+object.mesh.elements.m);
    for(int i=1;i<TV::m;i++) index_range_q[i]=INTERVAL<int>(index_range_q[i-1].max_corner,index_range_q[i-1].max_corner+object.mesh.elements.m);
    system_size=index_range_q[TV::m-1].max_corner;

    for(int i=0;i<TV::m;i++) index_map_u[i]->shift=index_range_u[i].min_corner;
    index_map_p->shift=index_range_p.min_corner;

    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++){
            helper_uu[i][j].Shift(index_range_u[i].min_corner,index_range_u[j].min_corner);
            helper.Add_Helper(helper_uu[i][j]);
            if(i!=j) helper.Add_Transpose();}

    for(int i=0;i<TV::m;i++){
        helper_p[i].Shift(index_range_u[i].min_corner,index_range_p.min_corner);
        helper.Add_Helper(helper_p[i]);
        helper.Add_Transpose();}

    for(int i=0;i<TV::m;i++){
        helper_q[i].Shift(index_range_u[i].min_corner,index_range_q[i].min_corner);
        helper.Add_Helper(helper_q[i]);
        helper.Add_Transpose();}

    for(int i=0;i<TV::m;i++){
        helper_rhs_p[i]->Shift(index_range_u[i].min_corner,index_range_p.min_corner);
        helper_rhs_q[i]->Shift(index_range_u[i].min_corner,index_range_q[i].min_corner);}

    helper.Set_Matrix(system_size,system_size,matrix,1e-14);

    for(int i=0;i<TV::m;i++){
        null_u[i].Resize(system_size);
        for(int j=index_range_u[i].min_corner;j<index_range_u[i].max_corner;j++)
            null_u[i](j)=1;
        null_u[i].Normalize();}

    null_p.Resize(system_size);
    for(int i=index_range_p.min_corner;i<index_range_p.max_corner;i++)
        null_p(i)=1;
    for(int i=0;i<TV::m;i++)
        for(int j=index_range_q[i].min_corner;j<index_range_q[i].max_corner;j++)
            null_p(j)=-object.Get_Element(j-index_range_q[i].min_corner).Normal()(i);
    null_p.Normalize();
}
//#####################################################################
// Function Set_RHS
//#####################################################################
template<class TV> void INTERFACE_FLUID_SYSTEM<TV>::
Set_RHS(VECTOR_T& rhs,const ARRAY<TV,TV_INT> f_body[2],const ARRAY<TV>& f_interface)
{
    VECTOR_ND<T> f(system_size);
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<f_interface.m;j++)
            f(j+index_range_q[i].min_corner)=f_interface(j)(i);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        for(int i=0;i<TV::m;i++){
            for(int s=0;s<2;s++){
                int index=index_map_u[i]->Get_Shifted_Index_Fixed(it.index,s);
                if(index>=0)
                    f(index)=f_body[s](it.index)(i);}}
    
    rhs.v.Resize(system_size);
    for(int i=0;i<TV::m;i++){
        for(int j=0;j<helper_rhs_q[i]->data.m;j++)
            rhs.v(helper_rhs_q[i]->data(j).x)+=helper_rhs_q[i]->data(j).z*f(helper_rhs_q[i]->data(j).y);
        delete helper_rhs_q[i];}
    for(int i=0;i<TV::m;i++){
        for(int j=0;j<helper_rhs_p[i]->data.m;j++)
            rhs.v(helper_rhs_p[i]->data(j).x)+=helper_rhs_p[i]->data(j).z*f(helper_rhs_p[i]->data(j).y);
        delete helper_rhs_p[i];}
}
//#####################################################################
// Function Resize_Vector
//#####################################################################
template<class TV> void INTERFACE_FLUID_SYSTEM<TV>::
Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const
{
    debug_cast<VECTOR_T&>(x).v.Resize(system_size);
}
//#####################################################################
// Function Get_U_Part
//#####################################################################
template<class TV> void INTERFACE_FLUID_SYSTEM<TV>::
Get_U_Part(const VECTOR_ND<T>& x,ARRAY<T,FACE_INDEX<TV::m> >& u) const
{
    u.Resize(grid);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        int s=phi.Phi(it.Location())<0;
        int index=index_map_u[it.Axis()]->Get_Shifted_Index_Fixed(it.index,s);
        assert(index>=0);
        u(it.Full_Index())=x(index);}
}
//#####################################################################
// Function Get_P_Part
//#####################################################################
template<class TV> void INTERFACE_FLUID_SYSTEM<TV>::
Get_P_Part(const VECTOR_ND<T>& x,ARRAY<T,TV_INT>& p) const
{
    p.Resize(grid.Domain_Indices());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int s=phi.Phi(it.Location())<0;
        int index=index_map_p->Get_Shifted_Index_Fixed(it.index,s);
        assert(index>=0);
        p(it.index)=x(index);}
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void INTERFACE_FLUID_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    matrix.Times(debug_cast<const VECTOR_T&>(x).v,debug_cast<VECTOR_T&>(result).v);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double INTERFACE_FLUID_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    return debug_cast<const VECTOR_T&>(x).v.Dot(debug_cast<const VECTOR_T&>(y).v);
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_FLUID_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    return debug_cast<const VECTOR_T&>(x).v.Max_Abs();
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void INTERFACE_FLUID_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    VECTOR_ND<T>& v=debug_cast<VECTOR_T&>(x).v;
    v-=v.Dot(null_p)*null_p;
    for(int i=0;i<TV::m;i++)
        v-=v.Dot(null_u[i])*null_u[i];
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void INTERFACE_FLUID_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void INTERFACE_FLUID_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project(x);
}
template class INTERFACE_FLUID_SYSTEM<VECTOR<float,2> >;
template class INTERFACE_FLUID_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class INTERFACE_FLUID_SYSTEM<VECTOR<double,2> >;
template class INTERFACE_FLUID_SYSTEM<VECTOR<double,3> >;
#endif
