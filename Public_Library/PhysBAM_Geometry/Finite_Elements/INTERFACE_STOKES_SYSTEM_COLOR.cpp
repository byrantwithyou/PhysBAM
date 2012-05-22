//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
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
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_HELPER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_COLOR<TV>::
INTERFACE_STOKES_SYSTEM_COLOR(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_value_input,ARRAY<int,TV_INT>& phi_color_input)
    :BASE(false,false),grid(grid_input),phi_grid(grid.counts*2,grid.domain,true),phi_value(phi_grid.Node_Indices()),phi_color(phi_grid.Node_Indices())
{
    // DEFINE DOUBLE FINE LEVELSET

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        phi_value(it.index*2)=phi_value_input(it.index);
        phi_color(it.index*2)=phi_color_input(it.index);}

    TV_INT counts(grid.counts+1);
    TV_INT scale(TV_INT()+2);
    for(int axis=0;axis<TV::m;axis++){
        counts(axis)=grid.counts(axis);
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),counts));it.Valid();it.Next()){
            TV_INT a=scale*it.index;
            TV_INT m(a); m+=TV_INT::Axis_Vector(axis);
            TV_INT b(m); b+=TV_INT::Axis_Vector(axis);
            const T& phi_value_a=phi_value(a);
            const T& phi_value_b=phi_value(b);
            phi_value(m)=(T).5*abs(phi_value_a-phi_value_b);
            phi_color(m)=(phi_value_a<phi_value_b)?phi_color(a):phi_color(b);}
        counts(axis)=phi_grid.counts(axis)+1;
        scale(axis)=1;}
    
    // PERTURB LEVELSET
    
    T panic_threshold=phi_grid.dX.Min()*1e-2;
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(phi_grid);it.Valid();it.Next()){
        T& value=phi_value(it.index);
        if(value<panic_threshold) value=panic_threshold;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_COLOR<TV>::
~INTERFACE_STOKES_SYSTEM_COLOR()
{
    for(int i=0;i<TV::m;i++) delete cm_u(i);
    delete cm_p;
    delete cdi;
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Set_Matrix(const ARRAY<T>& mu,bool wrap,ANALYTIC_BOUNDARY_CONDITIONS_COLOR<TV>* abc)
{
    // SET UP STENCILS

    BASIS_STENCIL_UNIFORM<TV,0> p_stencil(grid.dX);
    VECTOR<BASIS_STENCIL_UNIFORM<TV,1>*,TV::m> u_stencil;
    VECTOR<VECTOR<BASIS_STENCIL_UNIFORM<TV,1>*,TV::m>,TV::m> udx_stencil;
    p_stencil.Set_Center();
    p_stencil.Set_Constant_Stencil();
    p_stencil.Dice_Stencil();
    for(int i=0;i<TV::m;i++){
        u_stencil(i)=new BASIS_STENCIL_UNIFORM<TV,1>(grid.dX);
        u_stencil(i)->Set_Face(i);
        u_stencil(i)->Set_Multilinear_Stencil();
        u_stencil(i)->Dice_Stencil();
        for(int j=0;j<TV::m;j++){
            udx_stencil(i)(j)=new BASIS_STENCIL_UNIFORM<TV,1>(*u_stencil(i));
            udx_stencil(i)(j)->Differentiate(j);
            udx_stencil(i)(j)->Dice_Stencil();}}
    
    // GATHER CELL DOMAIN & INTERFACE INFO 

    int padding;
    if(wrap){
        padding=0;
        for(int i=0;i<TV::m;i++)
            for(int j=i;j<TV::m;j++)
                padding=max(u_stencil(i)->Overlap_Padding(*u_stencil(j)),padding);
        for(int i=0;i<TV::m;i++)
            padding=max(p_stencil.Overlap_Padding(*u_stencil(i)),padding);}
    else{
        padding=p_stencil.Padding();
        for(int i=0;i<TV::m;i++) padding=max(u_stencil(i)->Padding(),padding);}

    cdi=new CELL_DOMAIN_INTERFACE_COLOR<TV>(grid,padding,mu.m,wrap); 

    cm_p=new CELL_MANAGER_COLOR<TV>(*cdi);
    for(int i=0;i<TV::m;i++) cm_u(i)=new CELL_MANAGER_COLOR<TV>(*cdi);

    // STENCILS INTEGRATION 
    
    BASIS_INTEGRATION_UNIFORM_COLOR<TV,2> biu(grid,phi_grid,phi_value,phi_color,*cdi);
    VECTOR<VECTOR<SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>,TV::m>,TV::m> helper_uu;
    VECTOR<SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>,TV::m> helper_pu,helper_rhs_pu;
    VECTOR<SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>,TV::m> helper_qu;

    for(int i=0;i<TV::m;i++){
        for(int j=i;j<TV::m;j++)
            helper_uu(i)(j).Initialize(*u_stencil(i),*u_stencil(j),*cm_u(i),*cm_u(j),*cdi);
        helper_pu(i).Initialize(p_stencil,*u_stencil(i),*cm_p,*cm_u(i),*cdi);
        helper_qu(i).Initialize(*u_stencil(i),*cm_u(i),*cdi);
        helper_rhs_pu(i).Initialize(p_stencil,*u_stencil(i),*cm_p,*cm_u(i),*cdi);}

    ARRAY<T> double_mu(mu),ones(mu.m),minus_ones(mu.m);
    for(int i=0;i<double_mu.m;i++) double_mu(i)*=2;
    ones.Fill(1); minus_ones.Fill(-1);

    // Diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
            biu.Add_Volume_Block(helper_uu(i)(i),*udx_stencil(i)(j),*udx_stencil(i)(j),(i==j)?double_mu:mu);
    // Off-diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=i+1;j<TV::m;j++)
            biu.Add_Volume_Block(helper_uu(i)(j),*udx_stencil(i)(j),*udx_stencil(j)(i),mu);
    // Pressure blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_pu(i),p_stencil,*udx_stencil(i)(i),minus_ones);
    // Traction blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Surface_Block(helper_qu(i),*u_stencil(i),abc,i,1);
    // RHS pressure blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_rhs_pu(i),p_stencil,*u_stencil(i),ones);

    rhs_surface_forces=new VECTOR<ARRAY<VECTOR_ND<T> >,TV::m>;
    for(int i=0;i<TV::m;i++){
        (*rhs_surface_forces)(i).Resize(cdi->colors);
        for(int c=0;c<cdi->colors;c++)
            (*rhs_surface_forces)(i)(c).Resize(cdi->flat_size);}

    biu.Compute_Entries(rhs_surface_forces);
        
    // BUILD SYSTEM MATRIX BLOCKS
    
    for(int i=0;i<TV::m;i++){
        for(int j=i;j<TV::m;j++)
            helper_uu(i)(j).Mark_Active_Cells();
        helper_pu(i).Mark_Active_Cells();
        helper_qu(i).Mark_Active_Cells();
        helper_rhs_pu(i).Mark_Active_Cells();}
    
    cm_p->Compress_Indices();
    for(int i=0;i<TV::m;i++) cm_u(i)->Compress_Indices();

    // UU Block
    for(int i=0;i<TV::m;i++)
        for(int j=i;j<TV::m;j++)
            helper_uu(i)(j).Build_Matrix(matrix_uu(i)(j));
    // PU Block
    for(int i=0;i<TV::m;i++)
        helper_pu(i).Build_Matrix(matrix_pu(i));
    // QU Block
    for(int i=0;i<TV::m;i++)
        helper_qu(i).Build_Matrix(matrix_qu(i));
    // RHS PU Block 
    matrix_rhs_pu=new VECTOR<ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >,TV::m>;
    for(int i=0;i<TV::m;i++)
        helper_rhs_pu(i).Build_Matrix((*matrix_rhs_pu)(i));

    // FILL IN THE NULL MODES

/*    Resize_Vector(active_dofs);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            active_dofs.u(i)[s].Fill(1);
    for(int s=0;s<2;s++)
        active_dofs.p[s].Fill(1);
    active_dofs.q.Fill(1);

    if(BASE::use_preconditioner) Set_Jacobi_Preconditioner();

    Resize_Vector(null_p);
    for(int s=0;s<2;s++){
        null_p.p[s].Fill(1);
        int start=cdi->interface_dofs*(TV::m-1);
        int end=cdi->interface_dofs*TV::m;
        for(int k=start;k<end;k++)
            null_p.q(k)=-1;}
    null_p.Scale(active_dofs);
    null_p.Normalize();

    for(int i=0;i<TV::m;i++){
        Resize_Vector(null_u(i));
        for(int s=0;s<2;s++) null_u(i).u(i)[s].Fill(1);
        null_u(i).Scale(active_dofs);
        null_u(i).Normalize();}
    
    for(int i=0;i<TV::m;i++){
        delete u_stencil(i);
        for(int j=0;j<TV::m;j++)
        delete udx_stencil(i)(j);}*/
}
//#####################################################################
// Function Set_RHS
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Set_RHS(VECTOR_T& rhs,const ARRAY<ARRAY<TV,TV_INT> > f_body,const ARRAY<ARRAY<T,FACE_INDEX<TV::m> > >& u)
{
    /*VECTOR<VECTOR<VECTOR_ND<T>,2>,TV::m> F_body;
    VECTOR<VECTOR<VECTOR_ND<T>,2>,TV::m> U;
    
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            F_body(i)[s].Resize(cm_p->dofs[s]);
            U(i)[s].Resize(cm_u(i)->dofs[s]);}

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        for(int s=0;s<2;s++){
            int k=cm_p->Get_Index(it.index,s);
            if(k>=0)
                for(int i=0;i<TV::m;i++)
                    F_body(i)[s](k)=f_body[s](it.index)(i);}

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face(it.Full_Index()); 
        for(int s=0;s<2;s++){
            int k=cm_u(face.axis)->Get_Index(it.index,s);
            if(k>=0) U(face.axis)[s](k)=u[s](face);}}

    Resize_Vector(rhs); // assumes rhs was 0

    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            for(int j=0;j<cdi->flat_size;j++){
                int k=cm_u(i)->Get_Index(j,s);
                if(k>=0) rhs.u(i)[s](k)+=(*rhs_interface)(i)[s](j);}
    
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            (*matrix_rhs_pu)(i)[s].Transpose_Times_Add(F_body(i)[s],rhs.u(i)[s]);
            matrix_qu(i)[s].Times_Add(U(i)[s],rhs.q);
            matrix_pu(i)[s].Times_Add(U(i)[s],rhs.p[s]);}*/

    delete matrix_rhs_pu;
}
//#####################################################################
// Function Set_Jacobi_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Set_Jacobi_Preconditioner()
{
/*    Resize_Vector(J);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            int u_dofs=cm_u(i)->dofs[s];
            SPARSE_MATRIX_FLAT_MXN<T>& m_uu=matrix_uu(i)(i)[s];
            for(int k=0;k<u_dofs;k++){
                T d=abs(m_uu(k,k));
                if(d<1e-12) {active_dofs.u(i)[s](k)=0;LOG::cout<<"WARNING: small diagonal entry in the UU block."<<std::endl;}
                else J.u(i)[s](k)=1/abs(m_uu(k,k));}}
    for(int s=0;s<2;s++){
        for(int k=0;k<cm_p->dofs[s];k++){
            T sum=0;
            for(int i=0;i<TV::m;i++){
                SPARSE_MATRIX_FLAT_MXN<T>& m_pu=matrix_pu(i)[s];
                int start=m_pu.offsets(k);
                int end=m_pu.offsets(k+1);
                for(int e=start;e<end;e++)
                    sum+=sqr(m_pu.A(e).a)*J.u(i)[s](m_pu.A(e).j);}
            if(sum<1e-12) {active_dofs.p[s](k)=0;LOG::cout<<"WARNING: small row sum in the PU block."<<std::endl;}
            else J.p[s](k)=1/sum;}}
    for(int k=0;k<J.q.n;k++){
        T sum=0;
        for(int i=0;i<TV::m;i++)
            for(int s=0;s<2;s++){
                SPARSE_MATRIX_FLAT_MXN<T>& m_qu=matrix_qu(i)[s];
                int start=m_qu.offsets(k);
                int end=m_qu.offsets(k+1);
                for(int j=start;j<end;j++)
                    sum+=sqr(m_qu.A(j).a)*J.u(i)[s](m_qu.A(j).j);}
        if(sum<1e-12) {active_dofs.q(k)=0;LOG::cout<<"WARNING: small row sum in the QU block."<<std::endl;}
        else J.q(k)=1/sum;}*/
}
//#####################################################################
// Function Resize_Vector
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const
{
    VECTOR_T& v=debug_cast<VECTOR_T&>(x);
    for(int i=0;i<TV::m;i++){
        v.u(i).Resize(cdi->colors);
        for(int c=0;c<cdi->colors;c++)
            v.u(i)(c).Resize(cm_u(i)->dofs(c));}
    v.p.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++)
        v.p(c).Resize(cm_p->dofs(c));
    v.q.Resize(cdi->total_number_of_surface_constraints);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
/*    const VECTOR_T& xc=debug_cast<const VECTOR_T&>(x);
    VECTOR_T& rc=debug_cast<VECTOR_T&>(result);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            matrix_pu(i)[s].Transpose_Times(xc.p[s],rc.u(i)[s]);
            matrix_qu(i)[s].Transpose_Times_Add(xc.q,rc.u(i)[s]);
            for(int j=0;j<TV::m;j++){
                if(j>=i) matrix_uu(i)(j)[s].Times_Add(xc.u(j)[s],rc.u(i)[s]);
                else matrix_uu(j)(i)[s].Transpose_Times_Add(xc.u(j)[s],rc.u(i)[s]);}}
    for(int s=0;s<2;s++){
        rc.p[s].Fill(0);
        for(int i=0;i<TV::m;i++)
            matrix_pu(i)[s].Times_Add(xc.u(i)[s],rc.p[s]);}
    rc.q.Fill(0);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
        matrix_qu(i)[s].Times_Add(xc.u(i)[s],rc.q);*/
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    return debug_cast<const VECTOR_T&>(x).Dot(debug_cast<const VECTOR_T&>(y));
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    return debug_cast<const VECTOR_T&>(x).Max_Abs();
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    VECTOR_T& v=debug_cast<VECTOR_T&>(x);
    v.Copy(-v.Dot(null_p),null_p,v);
    for(int i=0;i<TV::m;i++)
        v.Copy(-v.Dot(null_u(i)),null_u(i),v);
    v.Scale(active_dofs);
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project(x);
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    debug_cast<VECTOR_T&>(z).Scale(debug_cast<const VECTOR_T&>(r),debug_cast<const VECTOR_T&>(J));
}
template class INTERFACE_STOKES_SYSTEM_COLOR<VECTOR<float,2> >;
template class INTERFACE_STOKES_SYSTEM_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INTERFACE_STOKES_SYSTEM_COLOR<VECTOR<double,2> >;
template class INTERFACE_STOKES_SYSTEM_COLOR<VECTOR<double,3> >;
#endif
