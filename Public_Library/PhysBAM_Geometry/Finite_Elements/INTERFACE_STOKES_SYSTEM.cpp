//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_INTERFACE_BLOCK_HELPER.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM<TV>::
INTERFACE_STOKES_SYSTEM(const GRID<TV>& grid_input,GRID<TV>& coarse_grid_input,ARRAY<T,TV_INT>& phi_input)
    :BASE(false,false),grid(grid_input),coarse_grid(coarse_grid_input),phi_grid(coarse_grid_input.Get_Regular_Grid().Get_MAC_Grid_At_Regular_Positions())
{
    phi=new LEVELSET_UNIFORM<GRID<TV> >(phi_grid,phi_input,0);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM<TV>::
~INTERFACE_STOKES_SYSTEM()
{
    for(int i=0;i<TV::m;i++) delete cm_u(i);
    delete cm_p;
    delete phi;
    delete cdi;
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
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

    int padding=p_stencil.Padding();
    for(int i=0;i<TV::m;i++) padding=max(padding,u_stencil(i)->Padding());

    MARCHING_CUBES<TV>::Create_Surface(object,coarse_grid,phi->phi);
    cdi=new CELL_DOMAIN_INTERFACE<TV>(grid,padding,coarse_grid.counts.x/grid.counts.x,object.mesh.elements.m,true); 

    cm_p=new CELL_MANAGER<TV>(*cdi);
    for(int i=0;i<TV::m;i++) cm_u(i)=new CELL_MANAGER<TV>(*cdi);

    // STENCILS INTEGRATION 
    
    BASIS_INTEGRATION_UNIFORM<TV,2> biu(grid,coarse_grid,phi->phi,*cdi);
    VECTOR<VECTOR<SYSTEM_VOLUME_BLOCK_HELPER<TV>,TV::m>,TV::m> helper_uu;
    VECTOR<SYSTEM_VOLUME_BLOCK_HELPER<TV>,TV::m> helper_pu;
    VECTOR<SYSTEM_INTERFACE_BLOCK_HELPER<TV>,TV::m> helper_qu,helper_rhs_qu;

    for(int i=0;i<TV::m;i++){
        for(int j=i;j<TV::m;j++)
            helper_uu(i)(j).Initialize(*u_stencil(i),*u_stencil(j),*cm_u(i),*cm_u(j),*cdi);
        helper_pu(i).Initialize(p_stencil,*u_stencil(i),*cm_p,*cm_u(i),*cdi);
        helper_qu(i).Initialize(*u_stencil(i),*cm_u(i),*cdi);
        helper_rhs_qu(i).Initialize(*u_stencil(i),*cm_u(i),*cdi);}

    // Diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++){
            biu.Add_Volume_Block(helper_uu(i)(i),*udx_stencil(i)(j),*udx_stencil(i)(j),mu*(1+(i==j)));}
    // Off-diagonal blocks
    for(int i=0;i<TV::m;i++)
        for(int j=i+1;j<TV::m;j++)
            biu.Add_Volume_Block(helper_uu(i)(j),*udx_stencil(i)(j),*udx_stencil(j)(i),mu);
    // Pressure blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Volume_Block(helper_pu(i),p_stencil,*udx_stencil(i)(i),VECTOR<T,2>(-1,-1));
    // Traction blocks
    for(int i=0;i<TV::m;i++)
        biu.Add_Interface_Block(helper_qu(i),*u_stencil(i),1,false);
    // RHS traction blocks
    for(int i=0;i<TV::m;i++){
        biu.Add_Interface_Block(helper_rhs_qu(i),*u_stencil(i),-0.5,true);}

    biu.Compute_Entries();

    // BUILD SYSTEM MATRIX BLOCKS

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
    // RHS QU Block 
    for(int i=0;i<TV::m;i++)
        helper_rhs_qu(i).Build_Matrix(matrix_rhs_qu(i));

    // FILL IN THE NULL MODES

    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            null_p.u(i)[s].Resize(cm_u(i)->dofs[s]);
    for(int s=0;s<2;s++)
        null_p.p[s].Resize(cm_p->dofs[s]);
    for(int i=0;i<TV::m;i++)
        null_p.q(i).Resize(object.mesh.elements.m);

    for(int s=0;s<2;s++)
        null_p.p[s].Fill(1);
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<object.mesh.elements.m;j++)
            null_p.q(i)(j)=object.Get_Element(j).Normal()(i);
    null_p.Normalize();

    for(int i=0;i<TV::m;i++){
        Resize_Vector(null_u(i));
        for(int s=0;s<2;s++) null_u(i).u(i)[s].Fill(1);
        null_u(i).Normalize();}
    
    if(BASE::use_preconditioner) Set_Jacobi_Preconditioner();

    for(int i=0;i<TV::m;i++){
        delete u_stencil(i);
        for(int j=0;j<TV::m;j++)
            delete udx_stencil(i)(j);}
}
//#####################################################################
// Function Set_RHS
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Set_RHS(VECTOR_T& rhs,const VECTOR<ARRAY<TV,TV_INT>,2> f_body,const ARRAY<TV>& f_interface)
{
    VECTOR<VECTOR_ND<T>,TV::m> F_interface;
    VECTOR<VECTOR<VECTOR_ND<T>,2>,TV::m> F_body;
    
    for(int i=0;i<TV::m;i++){
        F_interface(i).Resize(object.mesh.elements.m);
        for(int s=0;s<2;s++) F_body(i)[s].Resize(cm_p->dofs[s]);}

    for(int i=0;i<TV::m;i++)
        for(int k=0;k<object.mesh.elements.m;k++)
            F_interface(i)(k)=f_interface(k)(i);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        for(int s=0;s<2;s++){
            int k=cm_p->Get_Index(it.index,s);
            if(k>=0)
                for(int i=0;i<TV::m;i++)
                    F_body(i)[s](k)=f_body[s](it.index)(i);}

    Resize_Vector(rhs); // assumes rhs was 0
    
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            matrix_rhs_qu(i)[s].Transpose_Times(F_interface(i),rhs.u(i)[s]);
            matrix_pu(i)[s].Transpose_Times_Subtract(F_body(i)[s],rhs.u(i)[s]);}
}
//#####################################################################
// Function Set_Jacobi_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Set_Jacobi_Preconditioner()
{
    Resize_Vector(J);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            int u_dofs=cm_u(i)->dofs[s];
            SPARSE_MATRIX_FLAT_MXN<T>& m_uu=matrix_uu(i)(i)[s];
            for(int k=0;k<u_dofs;k++) J.u(i)[s](k)=1/abs(m_uu(k,k));}
    for(int s=0;s<2;s++){
        for(int k=0;k<cm_p->dofs[s];k++){
            T sum=0;
            for(int i=0;i<TV::m;i++){
                SPARSE_MATRIX_FLAT_MXN<T>& m_pu=matrix_pu(i)[s];
                int start=m_pu.offsets(k);
                int end=m_pu.offsets(k+1);
                for(int j=start;j<end;j++)
                    sum+=sqr(m_pu.A(j).a)*J.u(i)[s](m_pu.A(j).j);}
            J.p[s](k)=1/sum;}}
    for(int i=0;i<TV::m;i++)
        for(int k=0;k<object.mesh.elements.m;k++){
            T sum=0;
            for(int s=0;s<2;s++){
                SPARSE_MATRIX_FLAT_MXN<T>& m_qu=matrix_qu(i)[s];
                int start=m_qu.offsets(k);
                int end=m_qu.offsets(k+1);
                for(int j=start;j<end;j++)
                    sum+=sqr(m_qu.A(j).a)*J.u(i)[s](m_qu.A(j).j);}
            J.q(i)(k)=1/sum;}
}
//#####################################################################
// Function Resize_Vector
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Resize_Vector(KRYLOV_VECTOR_BASE<T>& x) const
{
    debug_cast<VECTOR_T&>(x).Resize(null_p);
}
//#####################################################################
// Function Get_U_Part
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Get_U_Part(const VECTOR_T& x,ARRAY<T,FACE_INDEX<TV::m> >& u) const
{
    u.Resize(grid);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        int i=it.Axis();
        int s=phi->Phi(it.Location())<0;
        int k=cm_u(it.Axis())->Get_Index(it.index,s);
        assert(k>=0);
        u(it.Full_Index())=x.u(i)[s](k);}
}
//#####################################################################
// Function Get_P_Part
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Get_P_Part(const VECTOR_T& x,ARRAY<T,TV_INT>& p) const
{
    p.Resize(grid.Domain_Indices());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int s=phi->Phi(it.Location())<0;
        int k=cm_p->Get_Index(it.index,s);
        assert(k>=0);
        p(it.index)=x.p[s](k);}
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    const VECTOR_T& xc=debug_cast<const VECTOR_T&>(x);
    VECTOR_T& rc=debug_cast<VECTOR_T&>(result);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            matrix_pu(i)[s].Transpose_Times(xc.p[s],rc.u(i)[s]);
            matrix_qu(i)[s].Transpose_Times_Add(xc.q(i),rc.u(i)[s]);
            for(int j=0;j<TV::m;j++)
                if(j>=i) matrix_uu(i)(j)[s].Times_Add(xc.u(j)[s],rc.u(i)[s]);
                else matrix_uu(j)(i)[s].Transpose_Times_Add(xc.u(j)[s],rc.u(i)[s]);}
    for(int s=0;s<2;s++){
        rc.p[s].Fill(0);
        for(int i=0;i<TV::m;i++)
            matrix_pu(i)[s].Times_Add(xc.u(i)[s],rc.p[s]);}
    for(int i=0;i<TV::m;i++){
        rc.q(i).Fill(0);
        for(int s=0;s<2;s++)
            matrix_qu(i)[s].Times_Add(xc.u(i)[s],rc.q(i));}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double INTERFACE_STOKES_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    return debug_cast<const VECTOR_T&>(x).Dot(debug_cast<const VECTOR_T&>(y));
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    return debug_cast<const VECTOR_T&>(x).Max_Abs();
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    // VECTOR_T& v=debug_cast<VECTOR_T&>(x);
    // v-=null_p*v.Dot(null_p);
    // for(int i=0;i<TV::m;i++)
        // v-=null_u(i)*v.Dot(null_u(i));
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project(x);
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    debug_cast<VECTOR_T&>(z).Scale(debug_cast<const VECTOR_T&>(r));
}
template class INTERFACE_STOKES_SYSTEM<VECTOR<float,2> >;
template class INTERFACE_STOKES_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class INTERFACE_STOKES_SYSTEM<VECTOR<double,2> >;
template class INTERFACE_STOKES_SYSTEM<VECTOR<double,3> >;
#endif
