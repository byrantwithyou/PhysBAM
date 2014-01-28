//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/CONSTANT_ARRAY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
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
#include <Geometry/Finite_Elements/CONSTRAINT_AGGREGATION_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_COLOR_NEW.h>
#include <Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR.h>
#include <Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
INTERFACE_POISSON_SYSTEM_COLOR_NEW(const GRID<TV>& grid_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input)
    :BASE(false,false),grid(grid_input),phi_grid(grid.counts*2,grid.domain,true),phi_value(phi_grid.Node_Indices()),phi_color(phi_grid.Node_Indices())
{
    T tol=(grid.dX/TV(grid.counts)).Min();
    CELL_DOMAIN_INTERFACE_COLOR<TV>::Interpolate_Level_Set_To_Double_Fine_Grid(grid_input,phi_value_input,phi_color_input,phi_grid,phi_value,phi_color,tol);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
~INTERFACE_POISSON_SYSTEM_COLOR_NEW()
{
    delete cm_u;
    delete cdi;
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Set_Matrix(const ARRAY<T>& mu,bool use_discontinuous_scalar_field,
    boost::function<T(const TV& X,int color0,int color1)> u_jump,
    boost::function<T(const TV& X,int color0,int color1)> j_surface,
    bool aggregated_constraints,bool cell_centered_u,bool eliminate_nullspace_input)
{
    // SET UP STENCILS

    BASIS_STENCIL_UNIFORM<TV,1> u_stencil;
    VECTOR<BASIS_STENCIL_UNIFORM<TV,1>,TV::m> udx_stencil;
    if(cell_centered_u)u_stencil.Set_Center(); else u_stencil.Set_Node();
    u_stencil.Set_Multilinear_Stencil(grid.dX);
    u_stencil.Dice_Stencil(grid.dX);
    for(int i=0;i<TV::m;i++){
        udx_stencil(i)=u_stencil;
        udx_stencil(i).Differentiate(i);
        udx_stencil(i).Dice_Stencil(grid.dX);}
    
    // GATHER CELL DOMAIN & INTERFACE INFO 
    eliminate_nullspace=eliminate_nullspace_input;

    int padding;
    padding=u_stencil.Overlap_Padding(u_stencil);

    cdi=new CELL_DOMAIN_INTERFACE_COLOR<TV>(grid,padding,mu.m); 
    cm_u=new CELL_MANAGER_COLOR<TV>(*cdi);
    cdi->Construct_Surface_Meshes(phi_grid,phi_value,phi_color);

    // STENCILS INTEGRATION 
    BASIS_INTEGRATION_UNIFORM_COLOR<TV,2> biu(grid,phi_grid,phi_color,*cdi);
    SYSTEM_VOLUME_BLOCK_HELPER_COLOR<TV> helper_uu,helper_rhs_uu;
    SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV> helper_qu;

    helper_uu.Initialize(u_stencil,u_stencil,*cm_u,*cm_u,*cdi);
    helper_rhs_uu.Initialize(u_stencil,u_stencil,*cm_u,*cm_u,*cdi);
    helper_qu.Initialize(u_stencil,*cm_u,*cdi);

    rhs_surface.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++)
        rhs_surface(c).Resize(cdi->flat_size);

    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_uu,udx_stencil(i),udx_stencil(i),mu);
    biu.Add_Surface_Block_Scalar(helper_qu,u_stencil,use_discontinuous_scalar_field,u_jump,j_surface,rhs_surface,1);
    biu.Add_Volume_Block(helper_rhs_uu,u_stencil,u_stencil,ARRAY<T>(CONSTANT_ARRAY<T>(mu.m,(T)1)));

    biu.Compute_Entries();//This is where the magic happens
    
//    for(int i=0;i<biu.all_constraint_offsets.m;i++){
//        for(int j=0;j<biu.all_constraint_offsets(i).m;j++){
//            std::cout<<biu.all_constraint_offsets(i)(j)<<" ";
//        }std::cout<<std::endl;
//    }

    // BUILD SYSTEM MATRIX BLOCKS
    
    helper_uu.Mark_Active_Cells();
    helper_qu.Mark_Active_Cells();
    helper_rhs_uu.Mark_Active_Cells();
    cm_u->Compress_Indices();

    helper_uu.Build_Matrix(matrix_uu);
    helper_qu.Build_Matrix(matrix_qu_full,rhs_constraint_full);
    helper_rhs_uu.Build_Matrix(matrix_rhs_uu);
    
    //This is where we get matrix_qu_agg from matrix_qu_full and decide which one becomes matrix_qu
    if(!aggregated_constraints){
        matrix_qu=matrix_qu_full; //For now
        rhs_constraint=rhs_constraint_full;
    }
    else
    {
        CONSTRAINT_AGGREGATION_COLOR<TV> cac(grid,cdi,cm_u,biu,matrix_qu_full,rhs_constraint_full);
        cac.Aggregate_Constraints(matrix_uu,phi_color);
        cac.Build_Condensed_Constraint_Matrix(matrix_qu_agg);
        matrix_qu=matrix_qu_agg;
        rhs_constraint=cac.RHS_Constraint_Agg();
        if(eliminate_nullspace){
            std::cout<<"\nNullspace Eliminated"<<std::endl;
            Resize_Vector(specific_solution_to_constraints);
            cac.Build_Nullspace_And_Specific_Constraint_Solution(matrix_uu,matrix_z,matrix_z_t,z_t_a,z_t_a_z,specific_solution_to_constraints);
            spd_system_size=z_t_a_z.m;

        }
    }
    //This is where it ends

    matrix_qu_t.Resize(cdi->colors);
#pragma omp parallel for
    for(int c=0;c<cdi->colors;c++)
        matrix_qu(c).Transpose(matrix_qu_t(c));

    // FILL IN THE NULL MODES

    inactive_u.Resize(cdi->colors);

    if(this->use_preconditioner) Set_Jacobi_Preconditioner();

    Resize_Vector(null_u);
    for(int c=0;c<cdi->colors;c++){
        ARRAY<T>& u=null_u.u(c);
        const ARRAY<int>& inactive=inactive_u(c);
        u.Fill(1);
        for(int k=0;k<inactive.m;k++) u(inactive(k))=0;}
    null_u.Normalize();
    if(eliminate_nullspace){
        null_u_condensed.v.Resize(spd_system_size);
        null_u_condensed.v.Fill(1);
        null_u_condensed.v.Normalize();
    }
}
//#####################################################################
// Function Build_Full_Solution_From_Condensed
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Build_Full_Solution_From_Condensed(KRYLOV_VECTOR_BASE<T>& small,KRYLOV_VECTOR_BASE<T>& big)
{
    //u=c+Zv
    VECTOR_T& sol=debug_cast<VECTOR_T&>(big);
    CONDENSED_VECTOR_T& condensed_sol=debug_cast<CONDENSED_VECTOR_T&>(small);
    Resize_Vector(sol);
    sol.Zero_Out();
    for(int c=0;c<cdi->colors;c++){
        matrix_z(c).Times(condensed_sol.v,sol.u(c)); 
    }
    sol+=specific_solution_to_constraints;
}
//#####################################################################
// Function Set_RHS
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Set_RHS(VECTOR_T& rhs,boost::function<T(const TV& X,int color)> body_force,bool cell_centered_u)
{
    ARRAY<ARRAY<T> > F_volume;

    Resize_Vector(rhs); // assumes rhs was 0
    rhs.q=rhs_constraint;
    
    F_volume.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++)
        F_volume(c).Resize(cm_u->dofs(c));
    if(cell_centered_u){
        for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next())
            for(int c=0;c<cdi->colors;c++){
                int k=cm_u->Get_Index(it.index,c);
                if(k>=0) F_volume(c)(k)=body_force(it.Location(),c);}}
    else{
        for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next())
            for(int c=0;c<cdi->colors;c++){
                int k=cm_u->Get_Index(it.index,c);
                if(k>=0) F_volume(c)(k)=body_force(it.Location(),c);}}
    for(int c=0;c<cdi->colors;c++)
        matrix_rhs_uu(c).Transpose_Times_Add(F_volume(c),rhs.u(c));

    for(int c=0;c<cdi->colors;c++)
        for(int j=0;j<cdi->flat_size;j++){
            int k=cm_u->Get_Index(j,c);
            if(k>=0) rhs.u(c)(k)+=rhs_surface(c)(j);}

    matrix_rhs_uu.Clean_Memory();
    rhs_surface.Clean_Memory();
    
    if(eliminate_nullspace){Create_Condensed_RHS(rhs);}
}
//#####################################################################
// Function Create_Condensed_RHS
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Create_Condensed_RHS(VECTOR_T& rhs)
{
    //condensed RHS is Z'*(rhs-Ac)
    condensed_rhs.v.Resize(spd_system_size);
    ARRAY<T> carrier;carrier.Resize(spd_system_size);
    for(int c=0;c<cdi->colors;c++){
        matrix_z_t(c).Times(rhs.u(c),carrier);condensed_rhs.v+=carrier;
        z_t_a(c).Times(specific_solution_to_constraints.u(c),carrier);condensed_rhs.v-=carrier;
    }
}
//#####################################################################
// Function Set_Jacobi_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Set_Jacobi_Preconditioner()
{
    if(eliminate_nullspace){
        J0.v.Resize(spd_system_size);
        for(int i=0;i<spd_system_size;i++){
            T d=1/abs(z_t_a_z(i,i));
            if(d<1e-13){LOG::cout<<"WARNING: small diagonal entry in Z'AZ."<<std::endl;
        }J0.v(i)=d;
        }
    }
    else{
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
    for(int k=0;k<J.q.m;k++){
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
}
//#####################################################################
// Function Resize_Vector
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const
{
    VECTOR_T& v=debug_cast<VECTOR_T&>(x);
    v.colors=cdi->colors;
    v.u.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++)
        v.u(c).Resize(cm_u->dofs(c));
    v.q.Resize(rhs_constraint.m);
}
//#####################################################################
// Function Resize_Condensed_Vector
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Resize_Condensed_Vector(KRYLOV_VECTOR_BASE<T>& x) const
{
    CONDENSED_VECTOR_T& v=debug_cast<CONDENSED_VECTOR_T&>(x);
    v.v.Resize(spd_system_size);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    
    if(eliminate_nullspace){
        const CONDENSED_VECTOR_T& xc=debug_cast<const CONDENSED_VECTOR_T&>(x);
        CONDENSED_VECTOR_T& rc=debug_cast<CONDENSED_VECTOR_T&>(result);
        rc.v.Fill(0);
#pragma omp parallel for
        for(int i=0;i<spd_system_size;i++)
            z_t_a_z.Times_Add_Row(xc.v,rc.v,i);
    }
    else{
        const VECTOR_T& xc=debug_cast<const VECTOR_T&>(x);
        VECTOR_T& rc=debug_cast<VECTOR_T&>(result);
        rc.Zero_Out();
        for(int c=0;c<cdi->colors;c++)
#pragma omp parallel for
            for(int i=0;i<cm_u->dofs(c);i++){
                matrix_uu(c).Times_Add_Row(xc.u(c),rc.u(c),i);
                matrix_qu_t(c).Times_Add_Row(xc.q,rc.u(c),i);}
        rc.q.Fill(0);
        for(int c=0;c<cdi->colors;c++)
#pragma omp parallel for
            for(int i=0;i<rhs_constraint.m;i++)
                matrix_qu(c).Times_Add_Row(xc.u(c),rc.q,i);
    }
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    if(eliminate_nullspace) return debug_cast<const CONDENSED_VECTOR_T&>(x).v.Dot(debug_cast<const CONDENSED_VECTOR_T&>(y).v);
    else return debug_cast<const VECTOR_T&>(x).Dot(debug_cast<const VECTOR_T&>(y));
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    if(eliminate_nullspace){return debug_cast<const CONDENSED_VECTOR_T&>(x).v.Max_Abs();}
    else return debug_cast<const VECTOR_T&>(x).Max_Abs();
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    
    if(eliminate_nullspace){
    //This will also need to change at some point. Also nodes aren't deactivated for this.
        CONDENSED_VECTOR_T& v=debug_cast<CONDENSED_VECTOR_T&>(x);
        if(!cdi->dc_present) v.Copy(-v.v.Dot(null_u_condensed.v),null_u_condensed,v);
        //This is where we sould then deal with inactive nodes
    }
    else{
        // TODO: This needs to change for N/D BC.
        VECTOR_T& v=debug_cast<VECTOR_T&>(x);
        if(!cdi->dc_present) v.Copy(-v.Dot(null_u),null_u,v);

        for(int c=0;c<cdi->colors;c++){
            ARRAY<T>& u=v.u(c);
            const ARRAY<int>& inactive=inactive_u(c);
#pragma omp parallel for
            for(int k=0;k<inactive.m;k++)
                u(inactive(k))=0;}
#pragma omp parallel for
        for(int k=0;k<inactive_q.m;k++)
            v.q(inactive_q(k))=0;
    }
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project(x);
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR_NEW<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    if(eliminate_nullspace){
       CONDENSED_VECTOR_T& zz=debug_cast<CONDENSED_VECTOR_T&>(z);
       const CONDENSED_VECTOR_T& rr=debug_cast<const CONDENSED_VECTOR_T&>(r);
#pragma omp parallel for
        for(int i=0;i<spd_system_size;i++)zz.v(i)=rr.v(i)*J0.v(i);
    }
    else debug_cast<VECTOR_T&>(z).Scale(debug_cast<const VECTOR_T&>(r),debug_cast<const VECTOR_T&>(J));
}
namespace PhysBAM{
template class INTERFACE_POISSON_SYSTEM_COLOR_NEW<VECTOR<float,2> >;
template class INTERFACE_POISSON_SYSTEM_COLOR_NEW<VECTOR<float,3> >;
template class INTERFACE_POISSON_SYSTEM_COLOR_NEW<VECTOR<double,2> >;
template class INTERFACE_POISSON_SYSTEM_COLOR_NEW<VECTOR<double,3> >;
}
