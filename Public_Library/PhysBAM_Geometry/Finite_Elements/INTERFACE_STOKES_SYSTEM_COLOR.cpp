//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
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
INTERFACE_STOKES_SYSTEM_COLOR(const GRID<TV>& grid_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input)
    :BASE(false,false),grid(grid_input),phi_grid(grid.counts*2,grid.domain,true),phi_value(phi_grid.Node_Indices()),phi_color(phi_grid.Node_Indices()),use_p_null_mode(false),use_u_null_mode(false)
{
    CELL_DOMAIN_INTERFACE_COLOR<TV>::Interpolate_Level_Set_To_Double_Fine_Grid(grid_input,phi_value_input,phi_color_input,phi_grid,phi_value,phi_color);
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
    null_modes.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Set_Matrix(const ARRAY<T>& mu,bool wrap,BOUNDARY_CONDITIONS_COLOR<TV>* abc,ARRAY<T>* system_inertia,ARRAY<T>* rhs_inertia)
{
    PHYSBAM_ASSERT((bool)system_inertia==(bool)rhs_inertia);
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
    VECTOR<SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>,TV::m> helper_pu,helper_rhs_pu,helper_inertial_rhs;
    VECTOR<SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>,TV::m> helper_qu;

    for(int i=0;i<TV::m;i++){
        for(int j=i;j<TV::m;j++)
            helper_uu(i)(j).Initialize(*u_stencil(i),*u_stencil(j),*cm_u(i),*cm_u(j),*cdi);
        helper_pu(i).Initialize(p_stencil,*u_stencil(i),*cm_p,*cm_u(i),*cdi);
        helper_qu(i).Initialize(*u_stencil(i),*cm_u(i),*cdi);
        helper_rhs_pu(i).Initialize(p_stencil,*u_stencil(i),*cm_p,*cm_u(i),*cdi);
        if(system_inertia)
            helper_inertial_rhs(i).Initialize(*u_stencil(i),*u_stencil(i),*cm_u(i),*cm_u(i),*cdi);}

    ARRAY<T> double_mu(mu*(T)2),ones(CONSTANT_ARRAY<T>(mu.m,(T)1)),minus_ones(CONSTANT_ARRAY<T>(mu.m,-(T)1));

    for(int i=0;i<TV::m;i++){
        rhs_surface(i).Resize(cdi->colors);
        for(int c=0;c<cdi->colors;c++)
            rhs_surface(i)(c).Resize(cdi->flat_size);}

    // Diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
            biu.Add_Volume_Block(helper_uu(i)(i),*udx_stencil(i)(j),*udx_stencil(i)(j),(i==j)?double_mu:mu);
    // Off-diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=i+1;j<TV::m;j++)
            biu.Add_Volume_Block(helper_uu(i)(j),*udx_stencil(i)(j),*udx_stencil(j)(i),mu);
    // Diagonal inertial term
    if(system_inertia){
        PHYSBAM_ASSERT(system_inertia->m==mu.m && system_inertia->m==rhs_inertia->m);
        for(int i=0;i<TV::m;i++){
            biu.Add_Volume_Block(helper_uu(i)(i),*u_stencil(i),*u_stencil(i),*system_inertia);
            biu.Add_Volume_Block(helper_inertial_rhs(i),*u_stencil(i),*u_stencil(i),*rhs_inertia);}}
    // Pressure blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_pu(i),p_stencil,*udx_stencil(i)(i),minus_ones);
    // Traction blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Surface_Block(helper_qu(i),*u_stencil(i),abc,rhs_surface(i),i,1);
    // RHS pressure blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_rhs_pu(i),p_stencil,*u_stencil(i),ones);

    biu.Compute_Entries();
        
    // BUILD SYSTEM MATRIX BLOCKS
    
    for(int i=0;i<TV::m;i++){
        for(int j=i;j<TV::m;j++)
            helper_uu(i)(j).Mark_Active_Cells();
        helper_pu(i).Mark_Active_Cells();
        helper_qu(i).Mark_Active_Cells();
        helper_rhs_pu(i).Mark_Active_Cells();
        if(system_inertia)
            helper_inertial_rhs(i).Mark_Active_Cells();}
    
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
        helper_qu(i).Build_Matrix(matrix_qu(i),q_rhs);
    // RHS PU Block 
    for(int i=0;i<TV::m;i++)
        helper_rhs_pu(i).Build_Matrix(matrix_rhs_pu(i));
    // RHS Inertial Block 
    if(system_inertia)
        for(int i=0;i<TV::m;i++)
            helper_inertial_rhs(i).Build_Matrix(matrix_inertial_rhs(i));
    // FILL IN THE NULL MODES

    for(int i=0;i<TV::m;i++) inactive_u(i).Resize(cdi->colors);
    inactive_p.Resize(cdi->colors);

    if(this->use_preconditioner) Set_Jacobi_Preconditioner();

    if(use_p_null_mode){
        null_modes.Append(new VECTOR_T);
        Resize_Vector(*null_modes.Last());
        for(int c=0;c<cdi->colors;c++){
            null_modes.Last()->p(c).Fill(1);
            int start=cdi->constraint_base_t*(TV::m-1);
            int end=start+cdi->constraint_base_n;
            for(int k=start;k<end;k++)
                null_modes.Last()->q(k)=-1;}}

    if(!system_inertia && use_u_null_mode)
        for(int i=0;i<TV::m;i++){
            null_modes.Append(new VECTOR_T);
            Resize_Vector(*null_modes.Last());
            for(int c=0;c<cdi->colors;c++)
                null_modes.Last()->u(i)(c).Fill(1);}

    for(int i=0;i<null_modes.m;i++){
        Clear_Unused_Entries(*null_modes(i));
        null_modes(i)->Normalize();}

    for(int i=0;i<TV::m;i++){
        delete u_stencil(i);
        for(int j=0;j<TV::m;j++)
            delete udx_stencil(i)(j);}
}
//#####################################################################
// Function Set_RHS
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Set_RHS(VECTOR_T& rhs,const ARRAY<ARRAY<TV,TV_INT> >& f_volume,const ARRAY<ARRAY<T,FACE_INDEX<TV::m> > >* u,bool analytic_velocity_correction)
{
    VECTOR<ARRAY<VECTOR_ND<T> >,TV::m> F_volume;
    
    Resize_Vector(rhs); // assumes rhs was 0
    rhs.q=q_rhs;

    for(int i=0;i<TV::m;i++){
        F_volume(i).Resize(cdi->colors);
        for(int c=0;c<cdi->colors;c++)
            F_volume(i)(c).Resize(cm_p->dofs(c));}

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        for(int c=0;c<cdi->colors;c++){
            int k=cm_p->Get_Index(it.index,c);
            if(k>=0)
                for(int i=0;i<TV::m;i++)
                    F_volume(i)(c)(k)=f_volume(c)(it.index)(i);}

    for(int i=0;i<TV::m;i++)
        for(int c=0;c<cdi->colors;c++)
            matrix_rhs_pu(i)(c).Transpose_Times_Add(F_volume(i)(c),rhs.u(i)(c));

    if(u){
        VECTOR<ARRAY<VECTOR_ND<T> >,TV::m> U;
        for(int i=0;i<TV::m;i++){
            U(i).Resize(cdi->colors);
            for(int c=0;c<cdi->colors;c++){
                U(i)(c).Resize(cm_u(i)->dofs(c));}}

        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
            FACE_INDEX<TV::m> face(it.Full_Index()); 
            for(int c=0;c<cdi->colors;c++){
                int k=cm_u(face.axis)->Get_Index(it.index,c);
                if(k>=0) U(face.axis)(c)(k)=(*u)(c)(face);}}

        for(int i=0;i<TV::m;i++)
            if(matrix_inertial_rhs(i).m)
                for(int c=0;c<cdi->colors;c++)
                    matrix_inertial_rhs(i)(c).Times_Add(U(i)(c),rhs.u(i)(c));

        if(analytic_velocity_correction)
            for(int i=0;i<TV::m;i++)
                for(int c=0;c<cdi->colors;c++){
                    matrix_qu(i)(c).Times_Add(U(i)(c),rhs.q);
                    matrix_pu(i)(c).Times_Add(U(i)(c),rhs.p(c));}}

    for(int i=0;i<TV::m;i++)
        for(int c=0;c<cdi->colors;c++)
            for(int j=0;j<cdi->flat_size;j++){
                int k=cm_u(i)->Get_Index(j,c);
                if(k>=0) rhs.u(i)(c)(k)+=rhs_surface(i)(c)(j);}

    for(int i=0;i<TV::m;i++) matrix_rhs_pu(i).Clean_Memory();
    for(int i=0;i<TV::m;i++) rhs_surface(i).Clean_Memory();
}
//#####################################################################
// Function Set_Jacobi_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Set_Jacobi_Preconditioner()
{
    Resize_Vector(J);
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<cdi->colors;c++){
            int u_dofs=cm_u(i)->dofs(c);
            SPARSE_MATRIX_FLAT_MXN<T>& m_uu=matrix_uu(i)(i)(c);
            for(int k=0;k<u_dofs;k++){
                T d=abs(m_uu(k,k));
                if(d<1e-13){
                    inactive_u(i)(c).Append(k);
                    LOG::cout<<"WARNING: small diagonal entry in the UU block."<<std::endl;}
                else J.u(i)(c)(k)=1/abs(m_uu(k,k));}}
    for(int c=0;c<cdi->colors;c++){
        for(int k=0;k<cm_p->dofs(c);k++){
            T sum=0;
            for(int i=0;i<TV::m;i++){
                SPARSE_MATRIX_FLAT_MXN<T>& m_pu=matrix_pu(i)(c);
                int start=m_pu.offsets(k);
                int end=m_pu.offsets(k+1);
                for(int e=start;e<end;e++)
                    sum+=sqr(m_pu.A(e).a)*J.u(i)(c)(m_pu.A(e).j);}
            if(sum<1e-13){
                inactive_p(c).Append(k);
                LOG::cout<<"WARNING: small row sum in the PU block."<<std::endl;}
            else J.p(c)(k)=1/sum;}}
    for(int k=0;k<J.q.n;k++){
        T sum=0;
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<cdi->colors;c++){
                SPARSE_MATRIX_FLAT_MXN<T>& m_qu=matrix_qu(i)(c);
                int start=m_qu.offsets(k);
                int end=m_qu.offsets(k+1);
                for(int j=start;j<end;j++)
                    sum+=sqr(m_qu.A(j).a)*J.u(i)(c)(m_qu.A(j).j);}
        if(sum<1e-13){
            inactive_q.Append(k);
            LOG::cout<<"WARNING: small row sum in the QU block."<<std::endl;}
        else J.q(k)=1/sum;}
}
//#####################################################################
// Function Resize_Vector
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const
{
    VECTOR_T& v=debug_cast<VECTOR_T&>(x);
    v.colors=cdi->colors;
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
    const VECTOR_T& xc=debug_cast<const VECTOR_T&>(x);
    VECTOR_T& rc=debug_cast<VECTOR_T&>(result);
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<cdi->colors;c++){
            matrix_pu(i)(c).Transpose_Times(xc.p(c),rc.u(i)(c));
            matrix_qu(i)(c).Transpose_Times_Add(xc.q,rc.u(i)(c));
            for(int j=0;j<TV::m;j++){
                if(j>=i) matrix_uu(i)(j)(c).Times_Add(xc.u(j)(c),rc.u(i)(c));
                else matrix_uu(j)(i)(c).Transpose_Times_Add(xc.u(j)(c),rc.u(i)(c));}}
    for(int c=0;c<cdi->colors;c++){
        rc.p(c).Fill(0);
        for(int i=0;i<TV::m;i++)
            matrix_pu(i)(c).Times_Add(xc.u(i)(c),rc.p(c));}
    rc.q.Fill(0);
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<cdi->colors;c++)
        matrix_qu(i)(c).Times_Add(xc.u(i)(c),rc.q);
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
// Function Clear_Unused_Entries
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Clear_Unused_Entries(VECTOR_T& v) const
{
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<cdi->colors;c++){
            VECTOR_ND<T>& u=v.u(i)(c);
            const ARRAY<int>& inactive=inactive_u(i)(c);
            for(int k=0;k<inactive.m;k++)
                u(inactive(k))=0;}

    for(int c=0;c<cdi->colors;c++)
        for(int k=0;k<inactive_p(c).m;k++)
            v.p(c)(inactive_p(c)(k))=0;

    for(int k=0;k<inactive_q.m;k++)
        v.q(inactive_q(k))=0;
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    // TODO: This needs to change for N/D BC.
    VECTOR_T& v=debug_cast<VECTOR_T&>(x);
    for(int i=0;i<null_modes.m;i++)
        v.Copy(-v.Dot(*null_modes(i)),*null_modes(i),v);

    Clear_Unused_Entries(v);
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
