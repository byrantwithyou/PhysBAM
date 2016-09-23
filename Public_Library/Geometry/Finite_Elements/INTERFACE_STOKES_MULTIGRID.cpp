//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/SPARSE_MATRIX_ROW.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_MULTIGRID.h>
#include <climits>
#ifdef USE_UMFPACK
#include <suitesparse/umfpack.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_MULTIGRID<TV>::
INTERFACE_STOKES_MULTIGRID(int num_levels,const GRID<TV>& grid_input,const ARRAY<T,TV_INT>& phi_value_input,
    const ARRAY<int,TV_INT>& phi_color_input,bool mac_phi,const ARRAY<ARRAY<T,TV_INT> >& phi_per_color_input,
    const VECTOR<ARRAY<T,TV_INT>,num_bc>& phi_boundary_input,int number_of_ghost_cells)
    :INTERFACE_STOKES_SYSTEM_COLOR<TV>(grid_input,phi_value_input,phi_color_input,mac_phi),
    levels(num_levels),boundary_smoother_iterations(5),number_of_ghost_cells(number_of_ghost_cells)
{
    levels(0).iss=this;
    levels(0).phi_per_color=phi_per_color_input;
    levels(0).phi_boundary=phi_boundary_input;

    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+2));it.Valid();it.Next())
        p_restriction_stencil.Append(it.index);

    for(int i=0;i<TV::m;i++){
        RANGE<TV_INT> range(TV_INT(),TV_INT()+2);
        range.min_corner(i)--;
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            u_restriction_stencil(i).Append(it.index);}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_STOKES_MULTIGRID<TV>::
~INTERFACE_STOKES_MULTIGRID()
{
    for(int i=1;i<levels.m;i++) delete levels(i).iss;
}
//#####################################################################
// Function Construct_Level
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Construct_Level(int l,const ARRAY<T>& mu,ARRAY<T>* inertia,T dt)
{
    ARRAY<ARRAY<T,TV_INT> >& cr_phis=levels(l).phi_per_color;
    ARRAY<ARRAY<T,TV_INT> >& bc_phis=levels(l).phi_boundary;
    if(l){
        const GRID<TV>& fine_grid=levels(l-1).iss->grid;
        TV_INT coarse_counts=(fine_grid.counts+1)/2;
        TV coarse_dX=fine_grid.dX*2;
        RANGE<TV> coarse_range(fine_grid.domain.min_corner,fine_grid.domain.min_corner+TV(coarse_counts)*coarse_dX);
        GRID<TV> coarse_grid(coarse_counts,coarse_range,true);

        Coarsen_Levelset(coarse_grid,levels(l-1).phi_per_color,cr_phis);
        Coarsen_Levelset(coarse_grid,levels(l-1).phi_boundary,bc_phis);
        Fill_Color_Levelset(coarse_grid,cr_phis,bc_phis,levels(l).color_levelset_phi,levels(l).color_levelset_color);

        levels(l).iss=new INTERFACE_STOKES_SYSTEM_COLOR<TV>(coarse_grid,levels(l).color_levelset_phi,
            levels(l).color_levelset_color,true);
        levels(l).iss->Set_Matrix(mu,false,0,0,inertia,false,dt);}

    levels(l).Initialize();
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Apply_Preconditioner(T_VECTOR& z,const T_VECTOR& x,bool initial_guess) const
{
    if(initial_guess){
        levels(0).iss->Multiply(z,levels(0).tmp1);
        levels(0).tmp2.Copy(-1,levels(0).tmp1,x);}
    else levels(0).tmp2=x;

    for(int i=0;i<levels.m-1;i++){
        levels(i).tmp0*=0;
        levels(i).Interior_Smoother(levels(i).tmp0,levels(i).tmp2);
        levels(i).Boundary_Smoother(levels(i).tmp0,levels(i).tmp2,boundary_smoother_iterations);
        levels(i).Interior_Smoother(levels(i).tmp0,levels(i).tmp2);
        levels(i).iss->Multiply(levels(i).tmp0,levels(i).tmp1);
        levels(i).tmp0.Copy(-1,levels(i).tmp1,levels(i).tmp2);
        levels(i).iss->Project(levels(i).tmp0);
        Restriction(levels(i+1).tmp2,levels(i).tmp2,i);}

    levels.Last().iss->Project(levels.Last().tmp0);
    levels.Last().Exact_Solve(levels.Last().tmp0,levels.Last().tmp2);

    for(int i=levels.m-2;i>=0;i--){
        levels.Last().iss->Project(levels.Last().tmp0);
        Prolongation(levels(i).tmp1,levels(i+1).tmp0,i);
        levels(i).tmp0.Copy(1,levels(i).tmp1,levels(i).tmp0);
        levels(i).Interior_Smoother(levels(i).tmp0,levels(i).tmp2);
        levels(i).Boundary_Smoother(levels(i).tmp0,levels(i).tmp2,boundary_smoother_iterations);
        levels(i).Interior_Smoother(levels(i).tmp0,levels(i).tmp2);}

        //z=levels(0).tmp0;
    if(initial_guess) z.Copy(1,levels(0).tmp0,z);
    else z=levels(0).tmp0;
}
//#####################################################################
// Function Restriction
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Restriction(T_VECTOR& z,const T_VECTOR& x,int fine_level) const
{
    INTERFACE_STOKES_SYSTEM_COLOR<TV>* coarse_iss=levels(fine_level+1).iss;
    INTERFACE_STOKES_SYSTEM_COLOR<TV>* fine_iss=levels(fine_level).iss;

    for(CELL_ITERATOR<TV> it(coarse_iss->grid);it.Valid();it.Next())
        for(int c=0;c<coarse_iss->cdi->colors;c++){
            int k=coarse_iss->cm_p->Get_Index(it.index,c); //fine_level + 1 = coarse_level
            if(k>=0){
                int weight=0;
                T value=0;
                for(int i=0;i<p_restriction_stencil.m;i++){
                    TV_INT index=it.index*2+p_restriction_stencil(i);
                    int l=fine_iss->cm_p->Get_Index(index,c);
                    if(l>=0){
                        weight++;
                        value+=x.p(c)(l);}}
                z.p(c)(k)=value/weight;}}

    for(FACE_ITERATOR<TV> it(coarse_iss->grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face(it.Full_Index());
        for(int c=0;c<coarse_iss->cdi->colors;c++){
            int k=coarse_iss->cm_u(face.axis)->Get_Index(it.index,c);
            if(k>=0){
                int weight=0;
                T value=0;
                for(int i=0;i<u_restriction_stencil(face.axis).m;i++){
                    TV_INT offset=u_restriction_stencil(face.axis)(i);
                    FACE_INDEX<TV::m> f(face.axis,face.index*2+offset);
                    int l=fine_iss->cm_u(f.axis)->Get_Index(f.index,c);
                    if(l>=0){
                        int w=1+!offset(f.axis);
                        weight+=w;
                        value+=x.u(f.axis)(c)(l)*w;}}
                z.u(face.axis)(c)(k)=value/weight;}}}
}
//#####################################################################
// Function Prolongation
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Prolongation(T_VECTOR& z,const T_VECTOR& x,int fine_level) const
{
    INTERFACE_STOKES_SYSTEM_COLOR<TV>* coarse_iss=levels(fine_level+1).iss;
    INTERFACE_STOKES_SYSTEM_COLOR<TV>* fine_iss=levels(fine_level).iss;

    for(CELL_ITERATOR<TV> it(coarse_iss->grid);it.Valid();it.Next())
        for(int c=0;c<coarse_iss->cdi->colors;c++){
            int k=coarse_iss->cm_p->Get_Index(it.index,c); //fine_level + 1 = coarse_level
            if(k>=0){
                VECTOR<const T*,1<<TV::m> from;
                VECTOR<T*,1<<TV::m> to;
                int n=0;
                for(int i=0;i<p_restriction_stencil.m;i++){
                    TV_INT index=it.index*2+p_restriction_stencil(i);
                    int l=fine_iss->cm_p->Get_Index(index,c);
                    if(l>=0){
                        to(n)=&z.p(c)(l);
                        from(n++)=&x.p(c)(k);}}
                for(int i=0;i<n;i++) *to(i)+=*from(i)*(1<<TV::m)/n;}}

    for(FACE_ITERATOR<TV> it(coarse_iss->grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face(it.Full_Index());
        for(int c=0;c<coarse_iss->cdi->colors;c++){
            int k=coarse_iss->cm_u(face.axis)->Get_Index(it.index,c);
            if(k>=0){
                VECTOR<const T*,3<<(TV::m-1)> from;
                VECTOR<T*,3<<(TV::m-1)> to;
                VECTOR<int,3<<(TV::m-1)> weights;
                int n=0;
                for(int i=0;i<u_restriction_stencil(face.axis).m;i++){
                    TV_INT offset=u_restriction_stencil(face.axis)(i);
                    FACE_INDEX<TV::m> f(face.axis,face.index*2+offset);
                    int l=fine_iss->cm_u(f.axis)->Get_Index(f.index,c);
                    if(l>=0){
                        weights(n)=1+!offset(f.axis);
                        to(n)=&z.u(f.axis)(c)(l);
                        from(n++)=&x.u(f.axis)(c)(k);}}
                int s=weights.Sum();
                for(int i=0;i<n;i++) *to(i)+=*from(i)*weights(i)*(2<<TV::m)/s;}}}
}
//#####################################################################
// Function Exact Solve
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::LEVEL::
Exact_Solve(T_VECTOR& z,const T_VECTOR& rhs) const
{
#ifdef USE_UMFPACK
    SPARSE_MATRIX_FLAT_MXN<T> M;
    iss->Get_Sparse_Matrix(M);
    ARRAY<double> rhs_umf(M.m);

    int track=0;
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<iss->cdi->colors;c++)
            for(int k=0;k<rhs.u(i)(c).m;k++)
                rhs_umf(track++)=rhs.u(i)(c)(k);
    for(int c=0;c<iss->cdi->colors;c++)
        for(int k=0;k<rhs.p(c).m;k++)
                rhs_umf(track++)=rhs.p(c)(k);
    for(int k=0;k<rhs.q.m;k++)
        rhs_umf(track++)=rhs.q(k);

    ARRAY<int> Mp;
    Mp.Append_Elements(M.offsets);
    ARRAY<int> Mi;
    ARRAY<double> Mx;

    for(int r=0;r<M.m;r++)
        for(int e=M.offsets(r);e<M.offsets(r+1);e++){
            Mx.Append(M.A(e).a);Mi.Append(M.A(e).j);}

        //for(int i=0;i<M.m;i++){
        //  for(int j=0;j<M.n;j++){
        //    if(M.Element_Present(i,j)){
        //       Mx.Append(M(i,j)); Mi.Append(j);}}}

    ARRAY<double> x_umf(M.m);
    //T *null_umf = (T *) NULL ;
    void *Symbolic_umf,*Numeric_umf;
    int status;
    double Control[UMFPACK_CONTROL],Info[UMFPACK_INFO];
    umfpack_di_defaults(Control);
    Control[UMFPACK_STRATEGY]=UMFPACK_STRATEGY_UNSYMMETRIC;
    status=umfpack_di_symbolic(M.n,M.n,Mp.Get_Array_Pointer(),Mi.Get_Array_Pointer(),Mx.Get_Array_Pointer(),
        &Symbolic_umf,Control,Info);
    PHYSBAM_ASSERT(status>=0);
    status=umfpack_di_numeric(Mp.Get_Array_Pointer(),Mi.Get_Array_Pointer(),Mx.Get_Array_Pointer(),Symbolic_umf,
        &Numeric_umf,Control,Info);
    PHYSBAM_ASSERT(status>=0);
    umfpack_di_free_symbolic (&Symbolic_umf);
    status=umfpack_di_solve(UMFPACK_A,Mp.Get_Array_Pointer(),Mi.Get_Array_Pointer(),Mx.Get_Array_Pointer(),
        x_umf.Get_Array_Pointer(),rhs_umf.Get_Array_Pointer(),Numeric_umf,Control,Info);
    PHYSBAM_ASSERT(status>=0);
    umfpack_di_free_numeric(&Numeric_umf);
    track=0;
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<iss->cdi->colors;c++)
            for(int k=0;k<rhs.u(i)(c).m;k++)
             z.u(i)(c)(k)=x_umf(track++);
    for(int c=0;c<iss->cdi->colors;c++)
        for(int k=0;k<rhs.p(c).m;k++)
            z.p(c)(k)=x_umf(track++);
    for(int k=0;k<rhs.q.m;k++)
        z.q(k)=x_umf(track++);
#else
    MINRES<T> mr;
    mr.Solve(*iss,z,rhs,av,1e-10,0,100);
#endif
}
//#####################################################################
// Function Interior_Smoother
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::LEVEL::
Interior_Smoother(T_VECTOR& z,const T_VECTOR& x) const
{
    ARRAY<T> x_flat(L.m),z_flat(L.m);
    
    int track=0;
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<iss->cdi->colors;c++)
            for(int k=0;k<x.u(i)(c).m;k++){
                x_flat(track)=x.u(i)(c)(k);
                z_flat(track++)=z.u(i)(c)(k);}
    for(int c=0;c<iss->cdi->colors;c++)
        for(int k=0;k<x.p(c).m;k++){
                x_flat(track)=x.p(c)(k);
                z_flat(track++)=z.p(c)(k);}
    for(int k=0;k<x.q.m;k++){
        x_flat(track)=x.q(k);
        z_flat(track++)=z.q(k);}

    int iterations=1;
    for(int i=0;i<iterations;i++)
        for(int k=0;k<interior_indices.m;k++){
            int index=interior_indices(k);
            SPARSE_MATRIX_ROW<T> LR(L,index),MR(M,index);
            T t=LR.Dot(MR),r=x_flat(index)-LR.Dot(z_flat);
            if(abs(t)>1e-7) MR.Add_Multiple(r/t,MR,z_flat);}

    track=0;
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<iss->cdi->colors;c++)
            for(int k=0;k<x.u(i)(c).m;k++)
                z.u(i)(c)(k)=z_flat(track++);
    for(int c=0;c<iss->cdi->colors;c++)
        for(int k=0;k<x.p(c).m;k++)
            z.p(c)(k)=z_flat(track++);
    for(int k=0;k<x.q.m;k++)
        z.q(k)=z_flat(track++);
}
//#####################################################################
// Function Boundary_Smoother
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::LEVEL::
Boundary_Smoother(T_VECTOR& z,const T_VECTOR& x,int iterations) const
{
    PHYSBAM_WARNING("Boundary_Smoother is not implemented.");
}
//#####################################################################
// Function Fill_Ghost
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Fill_Ghost(const GRID<TV>& grid,ARRAY<ARRAY<T,TV_INT> >& phi) const
{
    for(int i=0;i<phi.m;i++)
        Fill_Ghost(grid,phi(i));
}
//#####################################################################
// Function Fill_Ghost
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Fill_Ghost(const GRID<TV>& grid,ARRAY<T,TV_INT>& phi) const
{
    BOUNDARY_MAC_GRID_PERIODIC<TV,T>().Fill_Ghost_Cells(grid,phi,phi,0,0,number_of_ghost_cells);
}
//#####################################################################
// Function Coarsen_Levelset
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Coarsen_Levelset(const GRID<TV>& coarse_grid,const ARRAY<ARRAY<T,TV_INT> >& fine_phi,ARRAY<ARRAY<T,TV_INT> >& phi) const
{
    phi.Resize(fine_phi.m);
    for(int c=0;c<fine_phi.m;c++)
        Coarsen_Levelset(coarse_grid,fine_phi(c),phi(c));
}
//#####################################################################
// Function Coarsen_Levelset
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Coarsen_Levelset(const GRID<TV>& coarse_grid,const ARRAY<T,TV_INT>& fine_phi,ARRAY<T,TV_INT>& phi) const
{
    phi.Resize(coarse_grid.Cell_Indices(number_of_ghost_cells));
    for(CELL_ITERATOR<TV> it(coarse_grid);it.Valid();it.Next()){
        T value=0;
        for(int i=0;i<p_restriction_stencil.m;i++){
            TV_INT fine_index=it.index*2+p_restriction_stencil(i);
            value+=fine_phi(fine_index);}
        phi(it.index)=value/(1<<TV::m);}

    Fill_Ghost(coarse_grid,phi);
}
//#####################################################################
// Function Fill_Color_Levelset
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Fill_Color_Levelset(const GRID<TV>& grid,const ARRAY<ARRAY<T,TV_INT> >& cr_phis,const ARRAY<ARRAY<T,TV_INT> >& bc_phis,
    ARRAY<T,TV_INT>& color_phi,ARRAY<int,TV_INT>& colors) const
{
    color_phi.Resize(grid.Cell_Indices(number_of_ghost_cells));
    colors.Resize(grid.Cell_Indices(number_of_ghost_cells));

    for(CELL_ITERATOR<TV> it(grid,number_of_ghost_cells);it.Valid();it.Next()){
        T phi=0;
        int color=INT_MAX;
        for(int i=0;i<bc_phis.m;i++)
            if(bc_phis(i)(it.index)<=0){
                phi=-bc_phis(i)(it.index);
                color=~i;
                break;}

        if(color==INT_MAX)
            for(int k=0;k<cr_phis.m-1;k++){
                phi=cr_phis(k)(it.index);
                if(phi<=0){
                    color=k;
                    break;}}

        if(color==INT_MAX){
            color=cr_phis.m-1;
            phi=cr_phis(color)(it.index);}

        color_phi(it.index)=phi;
        colors(it.index)=color;}

    BOUNDARY_MAC_GRID_PERIODIC<TV,T>().Fill_Ghost_Cells(grid,color_phi,color_phi,0,0,number_of_ghost_cells);
    BOUNDARY_MAC_GRID_PERIODIC<TV,int>().Fill_Ghost_Cells(grid,colors,colors,0,0,number_of_ghost_cells);
}
//#####################################################################
// Function Get_Change_Of_Variables_Matrix
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::LEVEL::
Get_Change_Of_Variables_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& M) const
{
    const int colors=iss->cdi->colors;
    SPARSE_MATRIX_FLAT_MXN<T> temp;
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > gtg_poisson(colors); //per color

    for(int c=0;c<colors;c++){
        gtg_poisson(c).m=iss->matrix_pu(0)(c).m;
        gtg_poisson(c).n=iss->matrix_pu(0)(c).m;
        gtg_poisson(c).offsets.Resize(iss->matrix_pu(0)(c).m+1);
        temp.m=iss->matrix_pu(1)(c).n;
        temp.n=iss->matrix_pu(1)(c).m;
        for(int i=0;i<TV::m;i++){
            iss->matrix_pu(i)(c).Transpose(temp);
            temp*=(T)2;
            gtg_poisson(c)=gtg_poisson(c)+iss->matrix_pu(i)(c)*temp;}}
    
    int size=0;
    VECTOR<ARRAY<int>,TV::m> first_row_u;
    ARRAY<int> first_row_p;
    
    for(int i=0;i<TV::m;i++)
        for(int k=0;k<colors;k++){
            first_row_u(i).Append(size);
            size+=iss->matrix_uu(i)(i)(k).m;}
    int u_size=size;
    for(int k=0;k<colors;k++){
        first_row_p.Append(size);
        size+=iss->matrix_pu(0)(k).m;}
    
    M.m=size;
    M.n=size;

    ARRAY<int> row_entries(size),&next_entry=row_entries;

    for(int i=0;i<TV::m;i++)
        for(int k=0;k<colors;k++){
            const SPARSE_MATRIX_FLAT_MXN<T>& mat=iss->matrix_pu(i)(k);
            int first_row=first_row_p(k);
            int first_col=first_row_u(i)(k);
            for(int r=0;r<mat.m;r++){
                row_entries(first_row+r)+=mat.offsets(r+1)-mat.offsets(r);
                for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++)
                    row_entries(first_col+mat.A(e).j)++;}}

    for(int k=0;k<colors;k++){
        const SPARSE_MATRIX_FLAT_MXN<T>& mat=gtg_poisson(k);
        int first_row=first_row_p(k);
        for(int r=0;r<mat.m;r++)
            row_entries(first_row+r)+=mat.offsets(r+1)-mat.offsets(r);}
    
    M.offsets.Append(0);
    for(int i=0;i<row_entries.m;i++) M.offsets.Append(M.offsets.Last()+row_entries(i)+(int)(i<u_size)); //+1 for the identity block
    M.A.Resize(M.offsets.Last());
    next_entry=M.offsets;

    for(int i=0;i<u_size;i++)
        M.A(next_entry(i)++)=SPARSE_MATRIX_ENTRY<T>(i,1); //identity block

    for(int i=0;i<TV::m;i++)
        for(int k=0;k<colors;k++){
            const SPARSE_MATRIX_FLAT_MXN<T>& mat=iss->matrix_pu(i)(k);
            int first_row=first_row_p(k);
            int first_col=first_row_u(i)(k);
            for(int r=0;r<mat.m;r++)
                for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++){
                    M.A(next_entry(first_row+r)++)=SPARSE_MATRIX_ENTRY<T>(first_col+mat.A(e).j,mat.A(e).a);
                    M.A(next_entry(first_col+mat.A(e).j)++)=SPARSE_MATRIX_ENTRY<T>(first_row+r,mat.A(e).a);}}

    int p_offset=0;
    for(int k=0;k<colors;k++){
        const SPARSE_MATRIX_FLAT_MXN<T>& mat=gtg_poisson(k);
        int first_row=first_row_p(k);
        for(int r=0;r<mat.m;r++)
            for(int e=mat.offsets(r),end=mat.offsets(r+1);e<end;e++)
                M.A(next_entry(first_row+r)++)=SPARSE_MATRIX_ENTRY<T>(u_size+p_offset+mat.A(e).j,mat.A(e).a);
        p_offset+=mat.m;}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::LEVEL::
Initialize()
{
    iss->Resize_Vector(tmp0);
    iss->Resize_Vector(tmp1);
    iss->Resize_Vector(tmp2);
    iss->Get_Sparse_Matrix(L);
    Get_Change_Of_Variables_Matrix(M);
    interior_indices=IDENTITY_ARRAY<>(L.m);
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    Apply_Preconditioner(debug_cast<T_VECTOR&>(z),debug_cast<const T_VECTOR&>(r),false);
}
//#####################################################################
// Function Set_Matrix
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Set_Matrix(const ARRAY<T>& mu,bool use_discontinuous_velocity,std::function<TV(const TV& X,int color0,int color1)> u_jump,
    std::function<TV(const TV& X,int color0,int color1)> j_surface,ARRAY<T>* inertia,bool use_rhs,T dt)
{
    INTERFACE_STOKES_SYSTEM_COLOR<TV>::Set_Matrix(mu,use_discontinuous_velocity,u_jump,j_surface,inertia,use_rhs,dt);
    for(int i=0;i<levels.m;i++) Construct_Level(i,mu,inertia,dt);
}
namespace PhysBAM{
template class INTERFACE_STOKES_MULTIGRID<VECTOR<float,2> >;
template class INTERFACE_STOKES_MULTIGRID<VECTOR<float,3> >;
template class INTERFACE_STOKES_MULTIGRID<VECTOR<double,2> >;
template class INTERFACE_STOKES_MULTIGRID<VECTOR<double,3> >;
}
