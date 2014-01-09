//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_MULTIGRID.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_MULTIGRID<TV>::
INTERFACE_STOKES_MULTIGRID(int num_levels,INTERFACE_STOKES_SYSTEM_COLOR<TV>* iss,const ARRAY<ARRAY<T,TV_INT> >& phi_per_color_input,const ARRAY<ARRAY<T,TV_INT> >& phi_boundary_input)
    :levels(num_levels),boundary_smoother_iterations(5)
{
    levels(0).iss=iss;
    levels(0).phi_per_color=phi_per_color_input;
    levels(0).phi_boundary=phi_boundary_input;

    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+2));it.Valid();it.Next())
        p_restriction_stencil.Append(it.index);

    for(int i=0;i<TV::m;i++){
        RANGE<TV_INT> range(TV_INT(),TV_INT()+2);
        range.min_corner(i)--;
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next())
            u_restriction_stencil(i).Append(it.index);}

    for(int i=0;i<levels.m;i++){
        Construct_Level(i);}
                
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
Construct_Level(int l)
{    
    
    if(l){
        const GRID<TV>& fine_grid=levels(l-1).iss->grid;
        TV_INT fine_counts=fine_grid.counts;
        RANGE<TV> fine_domain=fine_grid.domain;
        TV fine_dX=fine_grid.dX;
        TV_INT coarse_counts=(fine_counts+1)/2;
        TV coarse_dX=fine_grid.dX/2;
        //RANGE<TV> coarse_range(fine_domain.min_corner,coarse_counts.Last()*coarse_dX.Last());
        RANGE<TV> coarse_range(fine_domain.min_corner,TV(coarse_counts)*coarse_dX);
        GRID<TV> coarse_grid;
        coarse_grid.Initialize(coarse_counts,coarse_range,true);
        
        const ARRAY<ARRAY<T,TV_INT> >& fine_phi=levels(l-1).phi_per_color;
        levels(l).phi_per_color.Resize(fine_phi.m);
        for(CELL_ITERATOR<TV> it(coarse_grid);it.Valid();it.Next())
            for(int c=0;c<fine_phi.m;c++){
                TV_INT coarse_index=it.index;
                //Insert(const T& element,const ID index)
                //levels(l).phi_per_color(c).Append(coarse_index);
                //levels(l).phi_per_color(c).Insert(0,coarse_index);
                T value=0;
                for(int i=0;i<p_restriction_stencil.m;i++){
                    TV_INT fine_index=it.index*2+p_restriction_stencil(i);
                    value+=fine_phi(c)(fine_index);}
                levels(l).phi_per_color(c)(coarse_index)=value/4;}

        const ARRAY<ARRAY<T,TV_INT> >& fine_phi_bdry=levels(l-1).phi_boundary;
        levels(l).phi_boundary.Resize(fine_phi_bdry.m);
        for(CELL_ITERATOR<TV> it(coarse_grid);it.Valid();it.Next())
            for(int c=0;c<fine_phi_bdry.m;c++){
                TV_INT coarse_index=it.index;
                //levels(l).phi_boundary(c).Append(coarse_index);
                T value=0;
                for(int i=0;i<p_restriction_stencil.m;i++){
                    TV_INT fine_index=it.index*2+p_restriction_stencil(i);
                    value+=fine_phi_bdry(c)(fine_index);}
                levels(l).phi_boundary(c)(coarse_index)=value/4;}

        const ARRAY<ARRAY<T,TV_INT> >& cr_phis=levels(l).phi_per_color;
        const ARRAY<ARRAY<T,TV_INT> >& bc_phis=levels(l).phi_boundary;
        T phi=0; //=cr_phis(cr_phis.m-1)(it.index);
        int color=cr_phis.m-1;
        for(CELL_ITERATOR<TV> it(coarse_grid,5);it.Valid();it.Next()){
            for(int i=0;i<bc_phis.m;i++){
                if(bc_phis(i)(it.index)<=0){
                    phi=-bc_phis(i)(it.index);
                    color=~i;
                    break;}
                for(int k=0;k<cr_phis.m-1;k++){
                    phi=cr_phis(k)(it.index);
                    if(bc_phis(k)(it.index)<=0){
                            color=k;
                            break;}
                    phi=cr_phis(cr_phis.m-1)(it.index);
                    color=cr_phis.m-1;}}
            levels(l).second_representation_phi_value(it.index)=phi;
            levels(l).second_representation_phi_color(it.index)=color;}
        
        // INTERFACE_STOKES_SYSTEM_COLOR<TV> temp_iss(coarse_grid,levels(l).second_representation_phi_value,levels(l).second_representation_phi_color,true);
        levels(l).iss=new INTERFACE_STOKES_SYSTEM_COLOR<TV>(coarse_grid,levels(l).second_representation_phi_value,levels(l).second_representation_phi_color,true);}
        //INTERFACE_STOKES_SYSTEM_COLOR<TV> *levels(l).iss(coarse_grid,levels(l).second_representation_phi_value,levels(l).second_representation_phi_color,true);}
    SPARSE_MATRIX_FLAT_MXN<T> M;
    levels(l).iss->Get_Sparse_Matrix(M);
    levels(l).interior_indices.Resize(M.m);
    for(int i=0;i<M.m;i++)
        levels(l).interior_indices(i)=i;

}

/*
template<class TV_input> void PLS_FC_EXAMPLE<TV_input>::
Rebuild_Levelset_Color()
{
    Make_Levelsets_Consistent();
    for(CELL_ITERATOR<TV> it(grid,number_of_ghost_cells);it.Valid();it.Next())
        levelset_color.color(it.index)=Color_At_Cell(it.index,levelset_color.phi(it.index));
    boundary.Fill_Ghost_Cells(grid,levelset_color.phi,levelset_color.phi,0,0,number_of_ghost_cells);
    boundary_int.Fill_Ghost_Cells(grid,levelset_color.color,levelset_color.color,0,0,number_of_ghost_cells);
}
        
template<class TV_input> int PLS_FC_EXAMPLE<TV_input>::
Color_At_Cell(const TV_INT& index,T& phi) const
{
    for(int i=0;i<bc_phis.m;i++)
        if(bc_phis(i)(index)<=0){
            phi=-bc_phis(i)(index);
            return ~i;}
    int c=particle_levelset_evolution_multiple.particle_levelset_multiple.levelset_multiple.Inside_Region(index,phi);
    phi=-phi;
    return c;
}    

template<class TV> int LEVELSET_MULTIPLE<TV>::
Inside_Region(const TV_INT& index,T& phi) const // assumes exactly one Phi<0 on a node
{
    for(int k=0;k<phis.m-1;k++){phi=Phi(k,index);if(phi<=0) return k;}
    phi=Phi(phis.m-1,index);
    assert(phi<=0);
    return phis.m-1;
}
*/
    // NEED TO GENERATE:
    // INTERFACE_STOKES_SYSTEM_COLOR<TV>* iss;
    // ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > pressure_poisson; // per color
    // ARRAY<int> interior_indices;
    // ARRAY<int> boundary_indices;
    //////////////////////////////////////////////////////////////////
 
//#####################################################################
// Function Update
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Update()
{
    for(int i=0;i<levels.m;i++){
        LEVEL& l=levels(i);
        l.iss->Resize_Vector(l.tmp0);
        l.iss->Resize_Vector(l.tmp1);
        l.iss->Resize_Vector(l.tmp2);}
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Apply_Preconditioner(T_VECTOR& z,const T_VECTOR& x,bool initial_guess)
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
    Exact_Solve(levels.Last().tmp0,levels.Last().tmp2);

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
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::
Exact_Solve(T_VECTOR& z,const T_VECTOR& rhs) const
{
    SPARSE_MATRIX_FLAT_MXN<T> M;
    levels.Last().iss->Get_Sparse_Matrix(M);
    ARRAY<double,int> rhs_umf(M.m,true,0);
    int track=0;
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<levels.Last().iss->cdi->colors;c++)
            for(int k=0;k<rhs.u(i)(c).m;k++)
                rhs_umf(track++)=rhs.u(i)(c)(k);
    for(int c=0;c<levels.Last().iss->cdi->colors;c++)
        for(int k=0;k<rhs.p(c).m;k++)
                rhs_umf(track++)=rhs.p(c)(k);
    for(int k=0;k<rhs.q.m;k++)
        rhs_umf(track++)=rhs.q(k);
    ARRAY<int> Mp;
    Mp.Append_Elements(M.offsets);
    ARRAY<int> Mi;
    ARRAY<double,int> Mx;
    for(int i=0;i<M.m;i++){
        for(int j=0;j<M.n;j++){
            if(M.Element_Present(i,j)){
                Mx.Append(M(i,j)); Mi.Append(j);}}}
    ARRAY<double,int> x_umf(M.m,true,0);
    //T *null_umf = (T *) NULL ;
    void *Symbolic_umf,*Numeric_umf;
    int status;
    double Control[UMFPACK_CONTROL],Info[UMFPACK_INFO];
    umfpack_di_defaults(Control);
    Control[UMFPACK_STRATEGY]=UMFPACK_STRATEGY_UNSYMMETRIC;
    status=umfpack_di_symbolic(M.n,M.n,Mp.Get_Array_Pointer(),Mi.Get_Array_Pointer(),Mx.Get_Array_Pointer(),&Symbolic_umf,
Control,Info);//std::cout <<"\n UMFPACK_Symbolic status = " << status << std::endl;
    PHYSBAM_ASSERT(status>=0);
    status=umfpack_di_numeric(Mp.Get_Array_Pointer(),Mi.Get_Array_Pointer(),Mx.Get_Array_Pointer(),Symbolic_umf,&Numeric_umf,Control,Info);//std::cout <<" UMFPACK_Numeric status = " <<  status << std::endl;
    PHYSBAM_ASSERT(status>=0);
    umfpack_di_free_symbolic (&Symbolic_umf);
    status=umfpack_di_solve(UMFPACK_A,Mp.Get_Array_Pointer(),Mi.Get_Array_Pointer(),Mx.Get_Array_Pointer(),x_umf.Get_Array_Pointer(),rhs_umf.Get_Array_Pointer(),Numeric_umf,Control,Info);//std::cout <<" UMFPACK_Solve status = " <<  status << std::endl;
    PHYSBAM_ASSERT(status>=0);
    umfpack_di_free_numeric(&Numeric_umf);
    track=0;
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<levels.Last().iss->cdi->colors;c++)
            for(int k=0;k<rhs.u(i)(c).m;k++)
             z.u(i)(c)(k)=x_umf(track++);
    for(int c=0;c<levels.Last().iss->cdi->colors;c++)
        for(int k=0;k<rhs.p(c).m;k++)
            z.p(c)(k)=x_umf(track++);
    for(int k=0;k<rhs.q.m;k++)
        z.q(k)=x_umf(track++);
}
//#####################################################################
// Function Interior_Smoother
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::LEVEL::
Interior_Smoother(T_VECTOR& z,const T_VECTOR& x) const
{
    
    SPARSE_MATRIX_FLAT_MXN<T> L,M; // L is iss operator, M is the change of basis matrix; M = [id -G; -G' -2*pressure_poisson]
        
    iss->Get_Sparse_Matrix(L);
    
    M.m=L.m;
    M.n=L.n;
    
    ARRAY<T> x_flat(L.m,true,0),z_flat(L.m,true,0);

    int track=0;
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<iss->cdi->colors;c++)
            for(int k=0;k<x.u(i)(c).m;k++){
                M(track,track)=1;
                x_flat(track)=x.u(i)(c)(k);
                z_flat(track++)=z.u(i)(c)(k);}
    int u_offset=track;
    for(int c=0;c<iss->cdi->colors;c++){
        int p_offset=0;
        for(int i=0;i<TV::m;i++){
            for(int k=0;k<iss->matrix_pu(i)(c).m;k++){
                for(int l=0;l<iss->matrix_pu(i)(c).n;l++){
                    M(u_offset+k,l+p_offset)=iss->matrix_pu(i)(c)(k,l);
                    M(l+p_offset,u_offset+k)=iss->matrix_pu(i)(c)(k,l);}
                x_flat(track)=x.p(c)(k);
                z_flat(track++)=z.p(c)(k);}
            p_offset+=iss->matrix_pu(i)(c).n;}}
    for(int k=0;k<x.q.m;k++){
        x_flat(track)=x.q(k);
        z_flat(track++)=z.q(k);}
    
    SPARSE_MATRIX_FLAT_MXN<T> temp;
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > gtg_poisson; //per color
    
    int last_offset=0;
    for(int c=0;c<iss->cdi->colors;c++){
        gtg_poisson(c).m=iss->matrix_pu(1)(c).m;
        gtg_poisson(c).n=iss->matrix_pu(1)(c).m;
        gtg_poisson(c)*=0;
        temp.m=iss->matrix_pu(1)(c).n;
        temp.n=iss->matrix_pu(1)(c).m;
        for(int i=0;i<TV::m;i++){
            iss->matrix_pu(i)(c).Transpose(temp);
            temp*=(T)2;
            gtg_poisson(c)=gtg_poisson(c)+temp*iss->matrix_pu(i)(c);}
        for(int i=0;i<gtg_poisson(c).m;i++){
            for(int j=0;j<gtg_poisson(c).n;j++){
                M(u_offset+last_offset+i,u_offset+last_offset+j)=gtg_poisson(c)(i,j);}}
        last_offset+=gtg_poisson(c).m;}

    // Smoother below //

    for(int k=0;k<interior_indices.m;k++){
        int index=interior_indices(k);
        T temp1=0;
        T temp2=0;
        for(int i=0;i<L.m;i++){
            temp1+=L(index,i)*M(i,index);
            temp2+=L(index,i)*z_flat(i);}
        T delta=(x_flat(index)-temp2)/temp1;
        for(int i=0;i<L.m;i++){
            z_flat(i)+=delta*M(i,index);}}

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

/*
    x = initial_guess;
    for i = 1:ITERATIONS
        for index = INTERIOR_INDICES
            temp = L(index,:)*M(:,index);
            r = RHS(index) - L(index,:)*x;
            delta  = r/temp;
            x += delta*M(:,index)*(abs(temp)>1e-7);
        endfor
    endfor */

//    PHYSBAM_FATAL_ERROR("TODO");
}
//#####################################################################
// Function Boundary_Smoother
//#####################################################################
template<class TV> void INTERFACE_STOKES_MULTIGRID<TV>::LEVEL::
Boundary_Smoother(T_VECTOR& z,const T_VECTOR& x,int iterations) const
{
// PHYSBAM_FATAL_ERROR("TODO");
}
namespace PhysBAM{
template class INTERFACE_STOKES_MULTIGRID<VECTOR<float,2> >;
template class INTERFACE_STOKES_MULTIGRID<VECTOR<float,3> >;
template class INTERFACE_STOKES_MULTIGRID<VECTOR<double,2> >;
template class INTERFACE_STOKES_MULTIGRID<VECTOR<double,3> >;
}
