//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/CONSTANT_ARRAY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Utilities/DEBUG_CAST.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM_COLOR.h>
#include <Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_HELPER_COLOR.h>
#include <Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_COLOR<TV>::
INTERFACE_STOKES_SYSTEM_COLOR(const GRID<TV>& grid_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,bool mac_phi)
    :BASE(false,false),grid(grid_input),phi_grid(grid.counts*2,grid.domain,true),phi_value(phi_grid.Node_Indices(mac_phi)),phi_color(phi_grid.Node_Indices(mac_phi)),
    use_p_null_mode(false),use_u_null_mode(false),use_polymer_stress(false)
{
    T tol=(grid.dX/TV(grid.counts)).Min();
    if(mac_phi) CELL_DOMAIN_INTERFACE_COLOR<TV>::Interpolate_Mac_Level_Set_To_Double_Fine_Grid(grid_input,phi_value_input,phi_color_input,phi_grid,phi_value,phi_color,tol);
    else CELL_DOMAIN_INTERFACE_COLOR<TV>::Interpolate_Level_Set_To_Double_Fine_Grid(grid_input,phi_value_input,phi_color_input,phi_grid,phi_value,phi_color,tol);
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
Set_Matrix(const ARRAY<T>& mu,bool use_discontinuous_velocity,boost::function<TV(const TV& X,int color0,int color1)> u_jump,
    boost::function<TV(const TV& X,int color0,int color1)> j_surface,ARRAY<T>* inertia,
    bool use_rhs)
{
    // SET UP STENCILS

    BASIS_STENCIL_UNIFORM<TV,0> p_stencil;
    VECTOR<BASIS_STENCIL_UNIFORM<TV,1>,TV::m> u_stencil;
    VECTOR<VECTOR<BASIS_STENCIL_UNIFORM<TV,1>,TV::m>,TV::m> udx_stencil;
    BASIS_STENCIL_UNIFORM<TV,1> polymer_stress_stencil;
    p_stencil.Set_Center();
    p_stencil.Set_Constant_Stencil(grid.dX);
    p_stencil.Dice_Stencil(grid.dX);
    polymer_stress_stencil.Set_Center();
    polymer_stress_stencil.Set_Multilinear_Stencil(grid.dX);
    polymer_stress_stencil.Dice_Stencil(grid.dX);
    for(int i=0;i<TV::m;i++){
        u_stencil(i).Set_Face(i);
        u_stencil(i).Set_Multilinear_Stencil(grid.dX);
        u_stencil(i).Dice_Stencil(grid.dX);
        for(int j=0;j<TV::m;j++){
            udx_stencil(i)(j)=u_stencil(i);
            udx_stencil(i)(j).Differentiate(j);
            udx_stencil(i)(j).Dice_Stencil(grid.dX);}}

    // GATHER CELL DOMAIN & INTERFACE INFO

    int padding;
    padding=0;
    for(int i=0;i<TV::m;i++)
        for(int j=i;j<TV::m;j++)
            padding=max(u_stencil(i).Overlap_Padding(u_stencil(j)),padding);
    for(int i=0;i<TV::m;i++)
        padding=max(p_stencil.Overlap_Padding(u_stencil(i)),padding);

    cdi=new CELL_DOMAIN_INTERFACE_COLOR<TV>(grid,padding,mu.m);
    cdi->Construct_Surface_Meshes(phi_grid,phi_value,phi_color);

    cm_p=new CELL_MANAGER_COLOR<TV>(*cdi);
    for(int i=0;i<TV::m;i++) cm_u(i)=new CELL_MANAGER_COLOR<TV>(*cdi);

    // STENCILS INTEGRATION

    BASIS_INTEGRATION_UNIFORM_COLOR<TV,2> biu(grid,phi_grid,phi_color,*cdi);
    VECTOR<VECTOR<SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>,TV::m>,TV::m> helper_uu;
    VECTOR<SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>,TV::m> helper_pu,helper_rhs_pu,helper_inertial_rhs;
    VECTOR<SYSTEM_SURFACE_BLOCK_HELPER_COLOR<TV>,TV::m> helper_qu;
    VECTOR<VECTOR<SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV>,TV::m>,TV::m> helper_polymer_stress_rhs; // u, then diff

    for(int i=0;i<TV::m;i++){
        for(int j=i;j<TV::m;j++)
            helper_uu(i)(j).Initialize(u_stencil(i),u_stencil(j),*cm_u(i),*cm_u(j),*cdi);
        helper_pu(i).Initialize(p_stencil,u_stencil(i),*cm_p,*cm_u(i),*cdi);
        helper_qu(i).Initialize(u_stencil(i),*cm_u(i),*cdi);
        if(use_rhs) helper_rhs_pu(i).Initialize(p_stencil,u_stencil(i),*cm_p,*cm_u(i),*cdi);
        if(inertia)
            if(use_rhs) helper_inertial_rhs(i).Initialize(u_stencil(i),u_stencil(i),*cm_u(i),*cm_u(i),*cdi);}
    if(use_rhs && use_polymer_stress)
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++)
                helper_polymer_stress_rhs(i)(j).Initialize(u_stencil(i),polymer_stress_stencil,*cm_u(i),*cm_p,*cdi);

    ARRAY<T> double_mu(mu*(T)2),ones(CONSTANT_ARRAY<T>(mu.m,(T)1)),minus_ones(CONSTANT_ARRAY<T>(mu.m,-(T)1));

    for(int i=0;i<TV::m;i++){
        rhs_surface(i).Resize(cdi->colors);
        for(int c=0;c<cdi->colors;c++)
            rhs_surface(i)(c).Resize(cdi->flat_size);}

    // Diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
            biu.Add_Volume_Block(helper_uu(i)(i),udx_stencil(i)(j),udx_stencil(i)(j),(i==j)?double_mu:mu);
    // Off-diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=i+1;j<TV::m;j++)
            biu.Add_Volume_Block(helper_uu(i)(j),udx_stencil(i)(j),udx_stencil(j)(i),mu);
    // Diagonal inertial term
    if(inertia){
        for(int i=0;i<TV::m;i++){
            biu.Add_Volume_Block(helper_uu(i)(i),u_stencil(i),u_stencil(i),*inertia);
            if(use_rhs) biu.Add_Volume_Block(helper_inertial_rhs(i),u_stencil(i),u_stencil(i),*inertia);}}
    // Pressure blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_pu(i),p_stencil,udx_stencil(i)(i),minus_ones);
    // Traction blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Surface_Block(helper_qu(i),u_stencil(i),use_discontinuous_velocity,u_jump,j_surface,rhs_surface(i),i,1);
    // RHS pressure blocks
    if(use_rhs)
        for(int i=0;i<TV::m;i++)
            biu.Add_Volume_Block(helper_rhs_pu(i),p_stencil,u_stencil(i),ones);
    // RHS polymer stress blocks
    if(use_rhs && use_polymer_stress)
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++)
                biu.Add_Volume_Block(helper_polymer_stress_rhs(i)(j),udx_stencil(i)(j),polymer_stress_stencil,ones);

    biu.Compute_Entries();

    // BUILD SYSTEM MATRIX BLOCKS

    for(int i=0;i<TV::m;i++){
        for(int j=i;j<TV::m;j++)
            helper_uu(i)(j).Mark_Active_Cells();
        helper_pu(i).Mark_Active_Cells();
        helper_qu(i).Mark_Active_Cells();
        if(use_rhs) helper_rhs_pu(i).Mark_Active_Cells();
        if(inertia && use_rhs)
            helper_inertial_rhs(i).Mark_Active_Cells();}
    if(use_polymer_stress && use_rhs)
        for(int i=0;i<TV::m;i++)
            for(int j=i;j<TV::m;j++)
                helper_polymer_stress_rhs(i)(j).Mark_Active_Cells();

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
    if(use_rhs){
        // RHS PU Block
        for(int i=0;i<TV::m;i++)
            helper_rhs_pu(i).Build_Matrix(matrix_rhs_pu(i));
        // RHS Inertial Block
        if(inertia)
            for(int i=0;i<TV::m;i++)
                helper_inertial_rhs(i).Build_Matrix(matrix_inertial_rhs(i));
        if(use_polymer_stress)
            for(int i=0;i<TV::m;i++)
                for(int j=0;j<TV::m;j++)
                    helper_polymer_stress_rhs(i)(j).Build_Matrix(matrix_polymer_stress_rhs(i)(j));}

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

    if(!inertia && use_u_null_mode)
        for(int i=0;i<TV::m;i++){
            null_modes.Append(new VECTOR_T);
            Resize_Vector(*null_modes.Last());
            for(int c=0;c<cdi->colors;c++)
                null_modes.Last()->u(i)(c).Fill(1);}

    for(int i=0;i<null_modes.m;i++){
        Clear_Unused_Entries(*null_modes(i));
        null_modes(i)->Normalize();}
}
//#####################################################################
// Function Set_RHS
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Set_RHS(VECTOR_T& rhs,boost::function<TV(const TV& X,int color)> body_force,const ARRAY<ARRAY<T,FACE_INDEX<TV::m> > >* u,bool analytic_velocity_correction)
{
    VECTOR<ARRAY<ARRAY<T> >,TV::m> F_volume;

    Resize_Vector(rhs); // assumes rhs was 0
    rhs.q=q_rhs;

    for(int i=0;i<TV::m;i++){
        F_volume(i).Resize(cdi->colors);
        for(int c=0;c<cdi->colors;c++)
            F_volume(i)(c).Resize(cm_p->dofs(c));}

    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next())
        for(int c=0;c<cdi->colors;c++){
            int k=cm_p->Get_Index(it.index,c);
            if(k>=0)
                for(int i=0;i<TV::m;i++)
                    F_volume(i)(c)(k)=body_force(it.Location(),c)(i);}

    for(int i=0;i<TV::m;i++)
        for(int c=0;c<cdi->colors;c++)
            matrix_rhs_pu(i)(c).Transpose_Times_Add(F_volume(i)(c),rhs.u(i)(c));

    if(u){
        VECTOR<ARRAY<ARRAY<T> >,TV::m> U;
        Pack(*u,U);

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

    // TODO: polymer_stress term here

    for(int i=0;i<TV::m;i++) matrix_rhs_pu(i).Clean_Memory();
    for(int i=0;i<TV::m;i++) rhs_surface(i).Clean_Memory();
}
//#####################################################################
// Function Add_Polymer_Stress_RHS
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Add_Polymer_Stress_RHS(VECTOR_T& rhs,const ARRAY<ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> >& polymer_stress)
{
    VECTOR<VECTOR<ARRAY<ARRAY<T> >,TV::m>,TV::m> S;
    Pack(polymer_stress,S);

    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
            for(int c=0;c<cdi->colors;c++)
                matrix_polymer_stress_rhs(i)(j)(c).Times_Add(S(i)(j)(c),rhs.u(i)(c));
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
                if(!d){
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
            if(!sum){
                inactive_p(c).Append(k);
                LOG::cout<<"WARNING: small row sum in the PU block."<<std::endl;}
            else J.p(c)(k)=1/sum;}}
    for(int k=0;k<J.q.m;k++){
        T sum=0;
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<cdi->colors;c++){
                SPARSE_MATRIX_FLAT_MXN<T>& m_qu=matrix_qu(i)(c);
                int start=m_qu.offsets(k);
                int end=m_qu.offsets(k+1);
                for(int j=start;j<end;j++)
                    sum+=sqr(m_qu.A(j).a)*J.u(i)(c)(m_qu.A(j).j);}
        if(!sum){
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
// Function Pack
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Pack(const ARRAY<ARRAY<T,FACE_INDEX<TV::m> > >& u,VECTOR<ARRAY<ARRAY<T> >,TV::m>& v) const
{
    for(int i=0;i<TV::m;i++){
        v(i).Resize(cdi->colors);
        for(int c=0;c<cdi->colors;c++){
            v(i)(c).Resize(cm_u(i)->dofs(c));}}

    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face(it.Full_Index());
        for(int c=0;c<cdi->colors;c++){
            int k=cm_u(face.axis)->Get_Index(it.index,c);
            if(k>=0) v(face.axis)(c)(k)=u(c)(face);}}
}
//#####################################################################
// Function Pack
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Pack(const ARRAY<ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> >& polymer_stress,VECTOR<VECTOR<ARRAY<ARRAY<T> >,TV::m>,TV::m>& S) const
{
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++){
            S(i)(j).Resize(cdi->colors);
            for(int c=0;c<cdi->colors;c++)
                S(i)(j)(c).Resize(cm_p->dofs(c));}

    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next())
        for(int c=0;c<cdi->colors;c++){
            int k=cm_p->Get_Index(it.index,c);
            if(k>=0){
                SYMMETRIC_MATRIX<T,TV::m> SM=polymer_stress(c)(it.index);
                for(int i=0;i<TV::m;i++)
                    for(int j=0;j<TV::m;j++)
                        S(i)(j)(c)(k)=SM(i,j);}}
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    const VECTOR_T& xc=debug_cast<const VECTOR_T&>(x);
    VECTOR_T& rc=debug_cast<VECTOR_T&>(result);

#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<cdi->colors;c++)
#pragma omp task
            {
                matrix_pu(i)(c).Transpose_Times(xc.p(c),rc.u(i)(c));
                matrix_qu(i)(c).Transpose_Times_Add(xc.q,rc.u(i)(c));
                for(int j=0;j<TV::m;j++){
                    if(j>=i) matrix_uu(i)(j)(c).Times_Add(xc.u(j)(c),rc.u(i)(c));
                    else matrix_uu(j)(i)(c).Transpose_Times_Add(xc.u(j)(c),rc.u(i)(c));}
            }
        for(int c=0;c<cdi->colors;c++)
#pragma omp task
        {
            rc.p(c).Fill(0);
            for(int i=0;i<TV::m;i++)
                matrix_pu(i)(c).Times_Add(xc.u(i)(c),rc.p(c));
        }
#pragma omp task
        {
            rc.q.Fill(0);
            for(int i=0;i<TV::m;i++)
                for(int c=0;c<cdi->colors;c++)
                    matrix_qu(i)(c).Times_Add(xc.u(i)(c),rc.q);
        }
    }
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
            ARRAY<T>& u=v.u(i)(c);
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
//#####################################################################
// Function Get_Sparse_Matrix
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_COLOR<TV>::
Get_Sparse_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& M) const
{
    int size=0;
    VECTOR<ARRAY<int>,TV::m> first_row_u;
    ARRAY<int> first_row_p;
    int first_row_q;
    for(int i=0;i<TV::m;i++)
        for(int k=0;k<matrix_uu(i)(i).m;k++){
            first_row_u(i).Append(size);
            size+=matrix_uu(i)(i)(k).m;}
    for(int k=0;k<matrix_pu(0).m;k++){
        first_row_p.Append(size);
        size+=matrix_pu(0)(k).m;}
    first_row_q=size;
    size+=matrix_qu(0).m?matrix_qu(0)(0).m:0;
    M.m=size;
    M.n=size;

    ARRAY<int> row_entries(size),&next_entry=row_entries;

    for(int i=0;i<TV::m;i++){
        for(int k=0;k<matrix_uu(i)(i).m;k++){
            const SPARSE_MATRIX_FLAT_MXN<T>& mat=matrix_uu(i)(i)(k);
            int first_row=first_row_u(i)(k);
            for(int r=0;r<mat.m;r++)
                row_entries(first_row+r)+=mat.offsets(r+1)-mat.offsets(r);}

        for(int j=i+1;j<TV::m;j++)
            for(int k=0;k<matrix_uu(i)(j).m;k++){
                const SPARSE_MATRIX_FLAT_MXN<T>& mat=matrix_uu(i)(j)(k);
                int first_row=first_row_u(i)(k);
                int first_col=first_row_u(j)(k);
                for(int r=0;r<mat.m;r++){
                    row_entries(first_row+r)+=mat.offsets(r+1)-mat.offsets(r);
                    for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++)
                        row_entries(first_col+mat.A(e).j)++;}}}

    for(int i=0;i<TV::m;i++)
        for(int k=0;k<matrix_pu(i).m;k++){
            const SPARSE_MATRIX_FLAT_MXN<T>& mat=matrix_pu(i)(k);
            int first_row=first_row_p(k);
            int first_col=first_row_u(i)(k);
            for(int r=0;r<mat.m;r++){
                row_entries(first_row+r)+=mat.offsets(r+1)-mat.offsets(r);
                for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++)
                    row_entries(first_col+mat.A(e).j)++;}}

    for(int i=0;i<TV::m;i++)
        for(int k=0;k<matrix_qu(i).m;k++){
            const SPARSE_MATRIX_FLAT_MXN<T>& mat=matrix_qu(i)(k);
            int first_row=first_row_q;
            int first_col=first_row_u(i)(k);
            for(int r=0;r<mat.m;r++){
                row_entries(first_row+r)+=mat.offsets(r+1)-mat.offsets(r);
                for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++)
                    row_entries(first_col+mat.A(e).j)++;}}

    M.offsets.Append(0);
    for(int i=0;i<row_entries.m;i++) M.offsets.Append(M.offsets.Last()+row_entries(i));
    M.A.Resize(M.offsets.Last());
    next_entry=M.offsets;

    for(int i=0;i<TV::m;i++){
        for(int k=0;k<matrix_uu(i)(i).m;k++){
            const SPARSE_MATRIX_FLAT_MXN<T>& mat=matrix_uu(i)(i)(k);
            int first_row=first_row_u(i)(k);
            int first_col=first_row_u(i)(k);
            for(int r=0;r<mat.m;r++)
                for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++)
                    M.A(next_entry(first_row+r)++)=SPARSE_MATRIX_ENTRY<T>(first_col+mat.A(e).j,mat.A(e).a);}

        for(int j=i+1;j<TV::m;j++)
            for(int k=0;k<matrix_uu(i)(j).m;k++){
                const SPARSE_MATRIX_FLAT_MXN<T>& mat=matrix_uu(i)(j)(k);
                int first_row=first_row_u(i)(k);
                int first_col=first_row_u(j)(k);
                for(int r=0;r<mat.m;r++)
                    for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++){
                        M.A(next_entry(first_row+r)++)=SPARSE_MATRIX_ENTRY<T>(first_col+mat.A(e).j,mat.A(e).a);
                        M.A(next_entry(first_col+mat.A(e).j)++)=SPARSE_MATRIX_ENTRY<T>(first_row+r,mat.A(e).a);}}}

    for(int i=0;i<TV::m;i++)
        for(int k=0;k<matrix_pu(i).m;k++){
            const SPARSE_MATRIX_FLAT_MXN<T>& mat=matrix_pu(i)(k);
            int first_row=first_row_p(k);
            int first_col=first_row_u(i)(k);
            for(int r=0;r<mat.m;r++)
                for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++){
                    M.A(next_entry(first_row+r)++)=SPARSE_MATRIX_ENTRY<T>(first_col+mat.A(e).j,mat.A(e).a);
                    M.A(next_entry(first_col+mat.A(e).j)++)=SPARSE_MATRIX_ENTRY<T>(first_row+r,mat.A(e).a);}}

    for(int i=0;i<TV::m;i++)
        for(int k=0;k<matrix_qu(i).m;k++){
            const SPARSE_MATRIX_FLAT_MXN<T>& mat=matrix_qu(i)(k);
            int first_row=first_row_q;
            int first_col=first_row_u(i)(k);
            for(int r=0;r<mat.m;r++)
                for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++){
                    M.A(next_entry(first_row+r)++)=SPARSE_MATRIX_ENTRY<T>(first_col+mat.A(e).j,mat.A(e).a);
                    M.A(next_entry(first_col+mat.A(e).j)++)=SPARSE_MATRIX_ENTRY<T>(first_row+r,mat.A(e).a);}}
}
namespace PhysBAM{
template class INTERFACE_STOKES_SYSTEM_COLOR<VECTOR<float,2> >;
template class INTERFACE_STOKES_SYSTEM_COLOR<VECTOR<float,3> >;
template class INTERFACE_STOKES_SYSTEM_COLOR<VECTOR<double,2> >;
template class INTERFACE_STOKES_SYSTEM_COLOR<VECTOR<double,3> >;
}
