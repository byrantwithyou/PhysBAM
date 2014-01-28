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
#include <Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_COLOR.h>
#include <Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR.h>
#include <Geometry/Finite_Elements/SYSTEM_VOLUME_BLOCK_HELPER_COLOR.h>
#include <Geometry/Finite_Elements/TRIPLE_JUNCTION_CORRECTION.h>
#include <Geometry/Finite_Elements/VOLUME_FORCE_SCALAR_COLOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_POISSON_SYSTEM_COLOR<TV>::
INTERFACE_POISSON_SYSTEM_COLOR(const GRID<TV>& grid_input,const ARRAY<ARRAY<T,TV_INT> >& color_phi_input)
    :BASE(false,false),grid(grid_input),phi_grid(grid.counts*2,grid.domain,true),color_phi(color_phi_input.m),ghost(3)
{
    T tol=(grid.dX/TV(grid.counts)).Min();
    for(int c=0;c<color_phi.m;c++){
        RANGE<TV_INT> domain(color_phi_input(c).domain*2);
        domain.max_corner-=1;
        color_phi(c).Resize(domain);
        CELL_DOMAIN_INTERFACE_COLOR<TV>::Interpolate_Level_Set_To_Double_Fine_Grid(color_phi_input(c).domain,color_phi_input(c),domain,color_phi(c),tol);}
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

template<class TV,class CELL_ELEMENTS> void Dump(const GRID<TV>& grid,const HASHTABLE<VECTOR<int,2>,CELL_ELEMENTS>& index_to_cell_elements)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,2> TV_INT;
    typedef typename MARCHING_CUBES_COLOR<TV>::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    typedef typename MARCHING_CUBES_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    VECTOR<T,3> color_map[4]={VECTOR<T,3>(0,0.7,0),VECTOR<T,3>(0.8,0.8,0),VECTOR<T,3>(0,0.4,1),VECTOR<T,3>(0.8,0.2,0)};
    T sep=(T).02;
    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(index_to_cell_elements);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<INTERFACE_ELEMENT>& interface_elements=cell_elements.interface;
        for(int i=0;i<interface_elements.m;i++){
            const INTERFACE_ELEMENT& V=interface_elements(i);
            if(V.color_pair.y>=0){
                if(V.color_pair.y>=0) Add_Debug_Object(V.face.X-V.face.Normal()*sep*grid.dX.Min(),color_map[V.color_pair.y]);
                if("Alexey was here") Add_Debug_Object(V.face.X+V.face.Normal()*sep*grid.dX.Min(),color_map[V.color_pair.x]);}
            else if(V.color_pair.x>=0) Add_Debug_Object(V.face.X-V.face.Normal()*sep*grid.dX.Min(),color_map[V.color_pair.x]);}}

    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(index_to_cell_elements);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<BOUNDARY_ELEMENT>& boundary_elements=cell_elements.boundary;
        for(int i=0;i<boundary_elements.m;i++){
            const BOUNDARY_ELEMENT& V=boundary_elements(i);
            Add_Debug_Object(V.face.X-V.face.Normal()*sep*grid.dX.Min(),color_map[V.color]);}}
}
template<class TV,class CELL_ELEMENTS> void Dump(const GRID<TV>& grid,const HASHTABLE<VECTOR<int,3>,CELL_ELEMENTS>& index_to_cell_elements){}

//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Set_Matrix(const ARRAY<T>& mu,BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>* abc)
{
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
    padding=u_stencil.Overlap_Padding(u_stencil);

    cdi=new CELL_DOMAIN_INTERFACE_COLOR<TV>(grid,padding,mu.m); 
    cm_u=new CELL_MANAGER_COLOR<TV>(*cdi);

    TRIPLE_JUNCTION_CORRECTION<TV> tjc(phi_grid.Get_Regular_Grid(),color_phi,ghost*2);
    tjc.Compute_Pairwise_Level_Set_Data();
    tjc.Cut_Interface(cdi->index_to_cell_elements);
    Dump(grid,cdi->index_to_cell_elements);
    Flush_Frame<TV>("cutting");

    // STENCILS INTEGRATION 
    
    BASIS_INTEGRATION_UNIFORM_COLOR<TV,2> biu(grid,phi_grid,tjc.combined_color,*cdi);
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

    biu.Compute_Entries();
        
    // BUILD SYSTEM MATRIX BLOCKS
    
    helper_uu.Mark_Active_Cells();
    helper_qu.Mark_Active_Cells();
    helper_rhs_uu.Mark_Active_Cells();
    
    cm_u->Compress_Indices();

    helper_uu.Build_Matrix(matrix_uu);
    helper_qu.Build_Matrix(matrix_qu,rhs_constraint);
    helper_rhs_uu.Build_Matrix(matrix_rhs_uu);

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

    for(int i=0;i<TV::m;i++) delete udx_stencil(i);
}
//#####################################################################
// Function Set_RHS
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_COLOR<TV>::
Set_RHS(VECTOR_T& rhs,VOLUME_FORCE_SCALAR_COLOR<TV>* vfsc)
{
    ARRAY<ARRAY<T> > F_volume;

    Resize_Vector(rhs); // assumes rhs was 0
    rhs.q=rhs_constraint;
    
    F_volume.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++)
        F_volume(c).Resize(cm_u->dofs(c));

    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next())
        for(int c=0;c<cdi->colors;c++){
            int k=cm_u->Get_Index(it.index,c);
            if(k>=0) F_volume(c)(k)=vfsc->F(it.Location(),c);}

    for(int c=0;c<cdi->colors;c++)
        matrix_rhs_uu(c).Transpose_Times_Add(F_volume(c),rhs.u(c));

    for(int c=0;c<cdi->colors;c++)
        for(int j=0;j<cdi->flat_size;j++){
            int k=cm_u->Get_Index(j,c);
            if(k>=0) rhs.u(c)(k)+=rhs_surface(c)(j);}

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
    rc.Zero_Out();
    for(int c=0;c<cdi->colors;c++)
#pragma omp parallel for
        for(int i=0;i<cm_u->dofs(c);i++){
            matrix_uu(c).Times_Add_Row(xc.u(c),rc.u(c),i);
            matrix_qu_t(c).Times_Add_Row(xc.q,rc.u(c),i);}
    rc.q.Fill(0);
    for(int c=0;c<cdi->colors;c++)
#pragma omp parallel for
        for(int i=0;i<cdi->constraint_base_scalar;i++)
            matrix_qu(c).Times_Add_Row(xc.u(c),rc.q,i);
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
        ARRAY<T>& u=v.u(c);
        const ARRAY<int>& inactive=inactive_u(c);
#pragma omp parallel for
        for(int k=0;k<inactive.m;k++)
            u(inactive(k))=0;}
#pragma omp parallel for
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
namespace PhysBAM{
template class INTERFACE_POISSON_SYSTEM_COLOR<VECTOR<float,2> >;
template class INTERFACE_POISSON_SYSTEM_COLOR<VECTOR<float,3> >;
template class INTERFACE_POISSON_SYSTEM_COLOR<VECTOR<double,2> >;
template class INTERFACE_POISSON_SYSTEM_COLOR<VECTOR<double,3> >;
}
