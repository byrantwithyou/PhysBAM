//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_POISSON_SYSTEM_COLOR<TV>::
INTERFACE_POISSON_SYSTEM_COLOR(const GRID<TV>& grid_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input)
    :BASE(false,false),grid(grid_input),phi_grid(grid.counts*2,grid.domain,true),phi_value(phi_grid.Node_Indices()),phi_color(phi_grid.Node_Indices())
{
    CELL_DOMAIN_INTERFACE_COLOR<TV>::Interpolate_Level_Set_To_Double_Fine_Grid(grid_input,phi_value_input,phi_color_input,phi_grid,phi_value,phi_color);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_POISSON_SYSTEM_COLOR<TV>::
~INTERFACE_POISSON_SYSTEM_COLOR()
{
    delete cm_u;
    delete cdi;
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Set_Matrix(const ARRAY<T>& mu,bool wrap,BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>* abc,bool double_fine)
{
    // SET LEVELSET EXACTLY ON DOUBLE FINE GRID

    if(double_fine){
        T tol=1e-2;
        T panic_threshold=phi_grid.dX.Min()*tol;
        ANALYTIC_TEST<TV> *at=debug_cast<ANALYTIC_TEST<TV>*>(abc);
        for(UNIFORM_GRID_ITERATOR_NODE<TV> it(phi_grid);it.Valid();it.Next()){
            phi_value(it.index)=max(at->phi_value(it.Location()),panic_threshold);
            phi_color(it.index)=at->phi_color(it.Location());}}

    // SET UP STENCILS

    BASIS_STENCIL_UNIFORM<TV,1> u_stencil(grid.dX);
    VECTOR<BASIS_STENCIL_UNIFORM<TV,1>*,TV::m> udx_stencil;
    u_stencil.Set_Center();
    u_stencil.Set_Multilinear_Stencil();
    u_stencil.Dice_Stencil();
    for(int i=0;i<TV::m;i++){
        udx_stencil(i)=new BASIS_STENCIL_UNIFORM<TV,1>(u_stencil);
        udx_stencil(i)->Differentiate(i);
        udx_stencil(i)->Dice_Stencil();}
    
    // GATHER CELL DOMAIN & INTERFACE INFO 

    int padding;
    if(wrap) padding=u_stencil.Overlap_Padding(u_stencil);
    else padding=u_stencil.Padding();

    cdi=new CELL_DOMAIN_INTERFACE_COLOR<TV>(grid,padding,mu.m,wrap); 
    cm_u=new CELL_MANAGER_COLOR<TV>(*cdi);

    // STENCILS INTEGRATION 
    
    BASIS_INTEGRATION_UNIFORM_COLOR<TV,2> biu(grid,phi_grid,phi_value,phi_color,*cdi);
    SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV> helper_uu,helper_rhs_uu;
    SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV> helper_qu;

    helper_uu.Initialize(u_stencil,u_stencil,*cm_u,*cm_u,*cdi);
    helper_rhs_uu.Initialize(u_stencil,u_stencil,*cm_u,*cm_u,*cdi);
    helper_qu.Initialize(u_stencil,*cm_u,*cdi);

    rhs_surface.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++)
        rhs_surface(c).Resize(cdi->flat_size);

    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_uu,*udx_stencil(i),*udx_stencil(i),mu);
    biu.Add_Surface_Block_Scalar(helper_qu,u_stencil,abc,rhs_surface,1);
    biu.Add_Volume_Block(helper_rhs_uu,u_stencil,u_stencil,ARRAY<T>(CONSTANT_ARRAY<T>(mu.m,(T)1)));

    biu.Compute_Entries(double_fine);
        
    // BUILD SYSTEM MATRIX BLOCKS
    
    helper_uu.Mark_Active_Cells();
    helper_qu.Mark_Active_Cells();
    helper_rhs_uu.Mark_Active_Cells();
    
    cm_u->Compress_Indices();

    helper_uu.Build_Matrix(matrix_uu);
    helper_qu.Build_Matrix(matrix_qu,rhs_constraint);
    helper_rhs_uu.Build_Matrix(matrix_rhs_uu);

    // FILL IN THE NULL MODES

    inactive_u.Resize(cdi->colors);

    if(this->use_preconditioner) Set_Jacobi_Preconditioner();

    Resize_Vector(null_u);
    for(int c=0;c<cdi->colors;c++){
        VECTOR_ND<T>& u=null_u.u(c);
        const ARRAY<int>& inactive=inactive_u(c);
        u.Fill(1);
        for(int k=0;k<inactive.m;k++) u(inactive(k))=0;}
    null_u.Normalize();

    for(int i=0;i<TV::m;i++) delete udx_stencil(i);
}
//#####################################################################
// Function Set_RHS
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Set_RHS(VECTOR_T& rhs,const ARRAY<ARRAY<T,TV_INT> >& f_volume,const ARRAY<ARRAY<T,TV_INT> >& u)
{
    ARRAY<VECTOR_ND<T> > F_volume,U;
    
    F_volume.Resize(cdi->colors);
    U.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++){
        F_volume(c).Resize(cm_u->dofs(c));
        U(c).Resize(cm_u->dofs(c));}

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        for(int c=0;c<cdi->colors;c++){
            int k=cm_u->Get_Index(it.index,c);
            if(k>=0)
                for(int i=0;i<TV::m;i++){
                    F_volume(c)(k)=f_volume(c)(it.index);
                    U(c)(k)=u(c)(it.index);}}

    Resize_Vector(rhs); // assumes rhs was 0
    rhs.q=rhs_constraint;

    for(int c=0;c<cdi->colors;c++)
        for(int j=0;j<cdi->flat_size;j++){
            int k=cm_u->Get_Index(j,c);
            if(k>=0) rhs.u(c)(k)+=rhs_surface(c)(j);}
    
    for(int c=0;c<cdi->colors;c++){
        matrix_rhs_uu(c).Transpose_Times_Add(F_volume(c),rhs.u(c));
        matrix_qu(c).Times_Add(U(c),rhs.q);}

    matrix_rhs_uu.Clean_Memory();
    rhs_surface.Clean_Memory();
}
//#####################################################################
// Function Set_Jacobi_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Set_Jacobi_Preconditioner()
{
    Resize_Vector(J);
    for(int c=0;c<cdi->colors;c++){
        int u_dofs=cm_u->dofs(c);
        SPARSE_MATRIX_FLAT_MXN<T>& m_uu=matrix_uu(c);
        for(int k=0;k<u_dofs;k++){
            T d=abs(m_uu(k,k));
            if(d<1e-13){
                inactive_u(c).Append(k);
                LOG::cout<<"WARNING: small diagonal entry in the UU block."<<std::endl;}
            else J.u(c)(k)=1/abs(m_uu(k,k));}}
    for(int k=0;k<J.q.n;k++){
        T sum=0;
        for(int c=0;c<cdi->colors;c++){
            SPARSE_MATRIX_FLAT_MXN<T>& m_qu=matrix_qu(c);
            int start=m_qu.offsets(k);
            int end=m_qu.offsets(k+1);
            for(int j=start;j<end;j++)
                sum+=sqr(m_qu.A(j).a)*J.u(c)(m_qu.A(j).j);}
        if(sum<1e-13){
            inactive_q.Append(k);
            LOG::cout<<"WARNING: small row sum in the QU block."<<std::endl;}
        else J.q(k)=1/sum;}
}
//#####################################################################
// Function Resize_Vector
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const
{
    VECTOR_T& v=debug_cast<VECTOR_T&>(x);
    v.colors=cdi->colors;
    v.u.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++)
        v.u(c).Resize(cm_u->dofs(c));
    v.q.Resize(cdi->constraint_base_scalar);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    const VECTOR_T& xc=debug_cast<const VECTOR_T&>(x);
    VECTOR_T& rc=debug_cast<VECTOR_T&>(result);
    for(int c=0;c<cdi->colors;c++){
        matrix_uu(c).Times(xc.u(c),rc.u(c));
        matrix_qu(c).Transpose_Times_Add(xc.q,rc.u(c));}
    rc.q.Fill(0);
    for(int c=0;c<cdi->colors;c++)
        matrix_qu(c).Times_Add(xc.u(c),rc.q);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    return debug_cast<const VECTOR_T&>(x).Dot(debug_cast<const VECTOR_T&>(y));
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    return debug_cast<const VECTOR_T&>(x).Max_Abs();
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    // TODO: This needs to change for N/D BC.
    VECTOR_T& v=debug_cast<VECTOR_T&>(x);
    if(!cdi->dc_present) v.Copy(-v.Dot(null_u),null_u,v);

    for(int c=0;c<cdi->colors;c++){
        VECTOR_ND<T>& u=v.u(c);
        const ARRAY<int>& inactive=inactive_u(c);
        for(int k=0;k<inactive.m;k++)
            u(inactive(k))=0;}
    for(int k=0;k<inactive_q.m;k++)
        v.q(inactive_q(k))=0;
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project(x);
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    debug_cast<VECTOR_T&>(z).Scale(debug_cast<const VECTOR_T&>(r),debug_cast<const VECTOR_T&>(J));
}
template class INTERFACE_POISSON_SYSTEM_COLOR<VECTOR<float,2> >;
template class INTERFACE_POISSON_SYSTEM_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INTERFACE_POISSON_SYSTEM_COLOR<VECTOR<double,2> >;
template class INTERFACE_POISSON_SYSTEM_COLOR<VECTOR<double,3> >;
#endif
