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
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK_HELPER_NEW.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_NEW.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_NEW<TV>::
INTERFACE_STOKES_SYSTEM_NEW(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,bool periodic_bc_input)
    :BASE(false,false),grid(grid_input),periodic_bc(periodic_bc_input),phi_grid(grid.counts*2,grid.domain,true),phi(phi_grid.Node_Indices())
{
    // DEFINE DOUBLE FINE LEVELSET

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()) phi(it.index*2)=phi_input(it.index);
    TV_INT counts(grid.counts+1);
    TV_INT scale(TV_INT()+2);
    for(int axis=0;axis<TV::m;axis++){
        counts(axis)=grid.counts(axis);
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),counts));it.Valid();it.Next()){
            TV_INT scaled_index=scale*it.index;
            phi(scaled_index+TV_INT::Axis_Vector(axis))=(T).5*(phi(scaled_index)+phi(scaled_index+TV_INT::Axis_Vector(axis)*2));}
        counts(axis)=phi_grid.counts(axis)+1;
        scale(axis)=1;}
    
    // PERTURB LEVELSET
    
    T panic_threshold=phi_grid.dX.Min()*1e-2;
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(phi_grid);it.Valid();it.Next()){
        T& phi_value=phi(it.index);
        if(fabs(phi_value)<panic_threshold) phi_value=phi_value<0?-panic_threshold:panic_threshold;}

    // GET THE NUMBER OF THE CUT CELLS

    cut_cells=0;
    int all_positive=(1<<(1<<TV::m))-1;
    const VECTOR<TV_INT,(1<<TV::m)>& phi_offsets=GRID<TV>::Binary_Counts(TV_INT(),2);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        TV_INT phi_base=it.index*2;
        int signs=0;
        for(int b=0;b<(1<<TV::m);b++) signs|=(phi(phi_base+phi_offsets(b))<0)<<b;
        if(!(signs==0||signs==all_positive)) cut_cells++;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_NEW<TV>::
~INTERFACE_STOKES_SYSTEM_NEW()
{
    for(int i=0;i<TV::m;i++) delete cm_u(i);
    delete cm_p;
    delete cdi;
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_NEW<TV>::
Set_Matrix(const VECTOR<T,2>& mu)
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
    if(periodic_bc){
        padding=0;
        for(int i=0;i<TV::m;i++)
            for(int j=i;j<TV::m;j++)
                padding=max(u_stencil(i)->Overlap_Padding(*u_stencil(j)),padding);
        for(int i=0;i<TV::m;i++)
            padding=max(p_stencil.Overlap_Padding(*u_stencil(i)),padding);}
    else{
        padding=p_stencil.Padding();
        for(int i=0;i<TV::m;i++) padding=max(u_stencil(i)->Padding(),padding);}
    
    // cdi=new CELL_DOMAIN_INTERFACE_NEW<TV>(grid,padding,cut_cells,periodic_bc); 

    cm_p=new CELL_MANAGER_NEW<TV>(*cdi);
    for(int i=0;i<TV::m;i++) cm_u(i)=new CELL_MANAGER_NEW<TV>(*cdi);

    // STENCILS INTEGRATION 
    
    BASIS_INTEGRATION_UNIFORM_NEW<TV,2> biu(grid,phi_grid,phi,*cdi);
    VECTOR<VECTOR<SYSTEM_VOLUME_BLOCK_HELPER_NEW<TV>,TV::m>,TV::m> helper_uu;
    VECTOR<SYSTEM_VOLUME_BLOCK_HELPER_NEW<TV>,TV::m> helper_pu,helper_rhs_pu;
    VECTOR<VECTOR<SYSTEM_INTERFACE_BLOCK_HELPER_NEW<TV>,TV::m>,TV::m> helper_qu;

    for(int i=0;i<TV::m;i++){
        for(int j=i;j<TV::m;j++){
            helper_uu(i)(j).Initialize(*u_stencil(i),*u_stencil(j),*cm_u(i),*cm_u(j),*cdi);
            helper_qu(i)(j).Initialize(*u_stencil(j),*cm_u(j),*cdi);}
        helper_pu(i).Initialize(p_stencil,*u_stencil(i),*cm_p,*cm_u(i),*cdi);
        helper_rhs_pu(i).Initialize(p_stencil,*u_stencil(i),*cm_p,*cm_u(i),*cdi);}
    
    // Diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
            biu.Add_Volume_Block(helper_uu(i)(i),*udx_stencil(i)(j),*udx_stencil(i)(j),mu*(1+(i==j)));
    // Off-diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=i+1;j<TV::m;j++)
            biu.Add_Volume_Block(helper_uu(i)(j),*udx_stencil(i)(j),*udx_stencil(j)(i),mu);
    // Pressure blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_pu(i),p_stencil,*udx_stencil(i)(i),VECTOR<T,2>(-1,-1));
    // Traction blocks
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
            biu.Add_Interface_Block(helper_qu(i)(j),*u_stencil(j),1,false);
    // RHS pressure blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_rhs_pu(i),p_stencil,*u_stencil(i),VECTOR<T,2>(1,1));

    biu.Compute_Entries();

    // BUILD SYSTEM MATRIX BLOCKS

    for(int i=0;i<TV::m;i++){
        for(int j=i;j<TV::m;j++){
            helper_uu(i)(j).Mark_Active_Cells();
            helper_qu(i)(j).Mark_Active_Cells();}
        helper_pu(i).Mark_Active_Cells();
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
        for(int j=i;j<TV::m;j++)
            helper_qu(i)(j).Build_Matrix(matrix_qu(i)(j));
    // RHS PU Block 
    for(int i=0;i<TV::m;i++)
        helper_rhs_pu(i).Build_Matrix(matrix_f_pu(i));

    // FILL IN THE NULL MODES

    Resize_Vector(active_dofs);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            active_dofs.u(i)[s].Fill(1);
    for(int s=0;s<2;s++)
        active_dofs.p[s].Fill(1);
    for(int i=0;i<TV::m;i++)
        active_dofs.q(i).Fill(1);

    if(BASE::use_preconditioner) Set_Jacobi_Preconditioner();

    Resize_Vector(null_p);
    for(int s=0;s<2;s++){
        null_p.p[s].Fill(1);
        null_p.q(0).Fill(-1);}
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
            delete udx_stencil(i)(j);}

    MARCHING_CUBES<TV>::Create_Surface(object,phi_grid,phi);
}
//#####################################################################
// Function Set_RHS
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_NEW<TV>::
Set_RHS(VECTOR_T& rhs,const VECTOR<ARRAY<TV,TV_INT>,2> f_body,const ARRAY<TV>& f_interface,const VECTOR<ARRAY<T,FACE_INDEX<TV::m> >,2>& u)
{
    // VECTOR<VECTOR_ND<T>,TV::m> F_interface;
    // VECTOR<VECTOR<VECTOR_ND<T>,2>,TV::m> F_body;
    // VECTOR<VECTOR<VECTOR_ND<T>,2>,TV::m> U;
    
    // for(int i=0;i<TV::m;i++){
        // F_interface(i).Resize(object.mesh.elements.m);
        // for(int s=0;s<2;s++){
            // F_body(i)[s].Resize(cm_p->dofs[s]);
            // U(i)[s].Resize(cm_u(i)->dofs[s]);}}

    // for(int i=0;i<TV::m;i++)
        // for(int k=0;k<object.mesh.elements.m;k++)
            // F_interface(i)(k)=f_interface(k)(i);
    // for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        // for(int s=0;s<2;s++){
            // int k=cm_p->Get_Index(it.index,s);
            // if(k>=0)
                // for(int i=0;i<TV::m;i++)
                    // F_body(i)[s](k)=f_body[s](it.index)(i);}

    // for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        // FACE_INDEX<TV::m> face(it.Full_Index()); 
        // for(int s=0;s<2;s++){
            // int k=cm_u(face.axis)->Get_Index(it.index,s);
            // if(k>=0) U(face.axis)[s](k)=u[s](face);}}

    // Resize_Vector(rhs); // assumes rhs was 0

    // for(int i=0;i<TV::m;i++)
        // for(int s=0;s<2;s++){
            // matrix_f_qu(i)[s].Transpose_Times_Add(F_interface(i),rhs.u(i)[s]);
            // matrix_f_pu(i)[s].Transpose_Times_Add(F_body(i)[s],rhs.u(i)[s]);
            // matrix_qu(i)[s].Times_Add(U(i)[s],rhs.q(i));
            // matrix_pu(i)[s].Times_Add(U(i)[s],rhs.p[s]);}
}
//#####################################################################
// Function Set_Jacobi_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_NEW<TV>::
Set_Jacobi_Preconditioner()
{
    // Resize_Vector(J);
    // for(int i=0;i<TV::m;i++)
        // for(int s=0;s<2;s++){
            // int u_dofs=cm_u(i)->dofs[s];
            // SPARSE_MATRIX_FLAT_MXN<T>& m_uu=matrix_uu(i)(i)[s];
            // for(int k=0;k<u_dofs;k++){
                // T d=abs(m_uu(k,k));
                // if(d<1e-12) {active_dofs.u(i)[s](k)=0;LOG::cout<<"WARNING: small diagonal entry in the UU block."<<std::endl;}
                // else J.u(i)[s](k)=1/abs(m_uu(k,k));}}
    // for(int s=0;s<2;s++){
        // for(int k=0;k<cm_p->dofs[s];k++){
            // T sum=0;
            // for(int i=0;i<TV::m;i++){
                // SPARSE_MATRIX_FLAT_MXN<T>& m_pu=matrix_pu(i)[s];
                // int start=m_pu.offsets(k);
                // int end=m_pu.offsets(k+1);
                // for(int j=start;j<end;j++)
                    // sum+=sqr(m_pu.A(j).a)*J.u(i)[s](m_pu.A(j).j);}
            // if(sum<1e-12) {active_dofs.p[s](k)=0;LOG::cout<<"WARNING: small row sum in the PU block."<<std::endl;}
            // else J.p[s](k)=1/sum;}}
    // for(int i=0;i<TV::m;i++)
        // for(int k=0;k<object.mesh.elements.m;k++){
            // T sum=0;
            // for(int s=0;s<2;s++){
                // SPARSE_MATRIX_FLAT_MXN<T>& m_qu=matrix_qu(i)[s];
                // int start=m_qu.offsets(k);
                // int end=m_qu.offsets(k+1);
                // for(int j=start;j<end;j++)
                    // sum+=sqr(m_qu.A(j).a)*J.u(i)[s](m_qu.A(j).j);}
            // if(sum<1e-12) {active_dofs.q(i)(k)=0;LOG::cout<<"WARNING: small row sum in the QU block."<<std::endl;}
            // else J.q(i)(k)=1/sum;}
}
//#####################################################################
// Function Resize_Vector
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_NEW<TV>::
Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const
{
    VECTOR_T& v=debug_cast<VECTOR_T&>(x);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            v.u(i)[s].Resize(cm_u(i)->dofs[s]);
    for(int s=0;s<2;s++)
        v.p[s].Resize(cm_p->dofs[s]);
    for(int i=0;i<TV::m;i++)
        v.q(i).Resize(cut_cells);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_NEW<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    const VECTOR_T& xc=debug_cast<const VECTOR_T&>(x);
    VECTOR_T& rc=debug_cast<VECTOR_T&>(result);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            matrix_pu(i)[s].Transpose_Times(xc.p[s],rc.u(i)[s]);
            for(int j=0;j<TV::m;j++){
                if(j>=i) matrix_uu(i)(j)[s].Times_Add(xc.u(j)[s],rc.u(i)[s]);
                else matrix_uu(j)(i)[s].Transpose_Times_Add(xc.u(j)[s],rc.u(i)[s]);
                matrix_qu(i)(j)[s].Transpose_Times_Add(xc.q(i),rc.u(j)[s]);}}
    for(int s=0;s<2;s++){
        rc.p[s].Fill(0);
        for(int i=0;i<TV::m;i++)
            matrix_pu(i)[s].Times_Add(xc.u(i)[s],rc.p[s]);}
    for(int i=0;i<TV::m;i++){
        rc.q(i).Fill(0);
        for(int j=0;j<TV::m;j++)
            for(int s=0;s<2;s++)
                matrix_qu(i)(j)[s].Times_Add(xc.u(j)[s],rc.q(i));}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double INTERFACE_STOKES_SYSTEM_NEW<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    return debug_cast<const VECTOR_T&>(x).Dot(debug_cast<const VECTOR_T&>(y));
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_NEW<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    return debug_cast<const VECTOR_T&>(x).Max_Abs();
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_NEW<TV>::
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
template<class TV> void INTERFACE_STOKES_SYSTEM_NEW<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_NEW<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project(x);
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_NEW<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    debug_cast<VECTOR_T&>(z).Scale(debug_cast<const VECTOR_T&>(r),debug_cast<const VECTOR_T&>(J));
}
template class INTERFACE_STOKES_SYSTEM_NEW<VECTOR<float,2> >;
template class INTERFACE_STOKES_SYSTEM_NEW<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INTERFACE_STOKES_SYSTEM_NEW<VECTOR<double,2> >;
template class INTERFACE_STOKES_SYSTEM_NEW<VECTOR<double,3> >;
#endif
