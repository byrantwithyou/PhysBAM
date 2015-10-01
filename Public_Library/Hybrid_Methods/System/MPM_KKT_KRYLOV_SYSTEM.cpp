//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar, Chenfanfu Jiang, Chuyuan Fu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_KKT_EXAMPLE.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/MPM_KKT_KRYLOV_SYSTEM.h>
#include <Hybrid_Methods/System/MPM_KKT_KRYLOV_VECTOR.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_KKT_KRYLOV_SYSTEM<TV>::
MPM_KKT_KRYLOV_SYSTEM(const MPM_KKT_EXAMPLE<TV>& example)
    :KRYLOV_SYSTEM_BASE<T>(false,false),example(example)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_KKT_KRYLOV_SYSTEM<TV>::
~MPM_KKT_KRYLOV_SYSTEM()
{
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MPM_KKT_KRYLOV_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    PHYSBAM_ASSERT(&BV!=&BF);
    const MPM_KKT_KRYLOV_VECTOR<TV>& V=debug_cast<const MPM_KKT_KRYLOV_VECTOR<TV>&>(BV);
    MPM_KKT_KRYLOV_VECTOR<TV>& F=debug_cast<MPM_KKT_KRYLOV_VECTOR<TV>&>(BF);
    F.u.Fill(TV());F.p.Fill((T)0);
    // -D^T*p
    int index=D.offsets(0);
    for(int row=0;row<D.m;row++){
        int end=D.offsets(row+1);T p_value=V.p.array(example.valid_pressure_indices(row));
        for(;index<end;index++){
            int col=D.A(index).j;
            F.u.array.Flattened()(col)-=D.A(index).a*p_value;}}
    // (1/dt)*M*u
    T one_over_dt=(T)1/example.dt;
    for(int t=0;t<example.valid_velocity_indices.m;t++){
        int cell=example.valid_velocity_indices(t);
        F.u.array(cell)+=example.mass.array(cell)*V.u.array(cell)*one_over_dt;}
    // -D*u
    index=D.offsets(0);
    for(int row=0;row<D.m;row++){
        int end=D.offsets(row+1);T sum=(T)0;
        for(;index<end;index++){
            int col=D.A(index).j;
            sum-=D.A(index).a*V.u.array.Flattened()(col);}
        F.p.array(example.valid_pressure_indices(row))+=sum;}
    //-(one_over_lambda/(J*dt))*p
    const T pressure_volume_scale=example.coarse_grid.DX().Product();
    for(int t=0;t<example.valid_pressure_indices.m;t++){
        TV_INT cell=example.valid_pressure_cell_indices(t);
        if(example.mass_coarse(cell))
            F.p(cell)-=pressure_volume_scale*one_over_dt*example.one_over_lambda(cell)/example.J(cell)*V.p(cell);}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MPM_KKT_KRYLOV_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const MPM_KKT_KRYLOV_VECTOR<TV>& X=debug_cast<const MPM_KKT_KRYLOV_VECTOR<TV>&>(x);
    const MPM_KKT_KRYLOV_VECTOR<TV>& Y=debug_cast<const MPM_KKT_KRYLOV_VECTOR<TV>&>(y);
    T r=0;
#pragma omp parallel for reduction(+:r)
    for(int k=0;k<example.valid_velocity_indices.m;k++){
        int i=example.valid_velocity_indices(k);
        r+=X.u.array(i).Dot(Y.u.array(i));}
#pragma omp parallel for reduction(+:r)
    for(int k=0;k<example.valid_pressure_indices.m;k++){
        int i=example.valid_pressure_indices(k);
        r+=X.p.array(i)*Y.p.array(i);}
    return r;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MPM_KKT_KRYLOV_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    return sqrt(Inner_Product(BR,BR));
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void MPM_KKT_KRYLOV_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void MPM_KKT_KRYLOV_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void MPM_KKT_KRYLOV_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BR) const
{
} 
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void MPM_KKT_KRYLOV_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Get_Appropriate_DF
//#####################################################################
template<class TV> ARRAY<ARRAY<TV,VECTOR<int,TV::m> >,VECTOR<int,TV::m> >
Get_Appropriate_DF(VECTOR<int,TV::m> element_index)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    static const T dwt1[][3]={{-(5.0/96),1.0/24,1.0/96},{-(7.0/16),0,7.0/16},{-(1.0/96),-(1.0/24),5.0/96}};
    static const T wt1[][3]={{17.0/768,5.0/128,1.0/768},{55.0/384,113.0/192,55.0/384},{1.0/768,5.0/128,17.0/768}};
    static const T dwt2[][3]={{-(7.0/24),1.0/12,5.0/24},{-(5.0/24),-(1.0/12),7.0/24}};
    static const T wt2[][3]={{5.0/48,1.0/3,1.0/16},{1.0/16,1.0/3,5.0/48}};
    TV_INT pressure_max_corner=TV_INT::Constant_Vector(2);
    for(int axis=0;axis<TV::m;axis++)
        if(element_index(axis)%2==0) pressure_max_corner(axis)++;
    RANGE<TV_INT> pressure_range(TV_INT(),pressure_max_corner);
    RANGE<TV_INT> velocity_range(TV_INT(),TV_INT::Constant_Vector(3));
    ARRAY<ARRAY<TV,TV_INT>,TV_INT> wt(pressure_range);
    for(RANGE_ITERATOR<TV::m> pit(pressure_range);pit.Valid();pit.Next()){
        wt(pit.index).Resize(velocity_range);
        wt(pit.index).Fill(TV::All_Ones_Vector());
        for(RANGE_ITERATOR<TV::m> vit(velocity_range);vit.Valid();vit.Next())
            for(int xyz=0;xyz<TV::m;xyz++)
                for(int paxis=0;paxis<TV::m;paxis++){
                    if(element_index(xyz)%2==0)
                        wt(pit.index)(vit.index)(xyz)*=xyz==paxis?dwt1[pit.index(paxis)][vit.index(paxis)]:wt1[pit.index(paxis)][vit.index(paxis)];
                    else
                        wt(pit.index)(vit.index)(xyz)*=xyz==paxis?dwt2[pit.index(paxis)][vit.index(paxis)]:wt2[pit.index(paxis)][vit.index(paxis)];}}
    return wt;
}
//#####################################################################
// Function Precompute_Div_Coefficients
//#####################################################################
template<class TV> ARRAY<TV,VECTOR<int,TV::m> > 
Precompute_Div_Coefficients()
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    static const T w[]={1.0/768,13.0/128,383.0/768,51.0/64,383.0/768,13.0/128,1.0/768};
    static const T dw[]={-1.0/96,-1.0/4,-15.0/32,0,15.0/32,1.0/4,1.0/96};
    RANGE<TV_INT> local_range(-3*TV_INT::All_Ones_Vector(),4*TV_INT::All_Ones_Vector());
    ARRAY<TV,TV_INT> wt(local_range,true,TV::All_Ones_Vector());
    for(RANGE_ITERATOR<TV::m> it(local_range);it.Valid();it.Next()){
        for(int axis=0;axis<TV::m;axis++)
            for(int paxis=0;paxis<TV::m;paxis++)
                wt(it.index)(axis)*=axis==paxis?dw[it.index(paxis)+3]:w[it.index(paxis)+3];}
    return wt;
}
//#####################################################################
// Function Build_Div_Matrix
//#####################################################################
template<class TV> void MPM_KKT_KRYLOV_SYSTEM<TV>::
Build_Div_Matrix()
{
    PHYSBAM_ASSERT(example.weights->Order()==2);
    static const ARRAY<TV,TV_INT> wt=Precompute_Div_Coefficients<TV>();
    const T dx_scale=pow(example.grid.DX()(0),TV::m-1);
    D.Reset(TV::m*example.velocity.array.m);
    for(int row=0;row<example.valid_pressure_indices.m;row++){
        TV_INT pid=example.valid_pressure_cell_indices(row);
        TV_INT vid=2*pid;
        RANGE<TV_INT> local_range(-4*TV_INT::All_Ones_Vector(),5*TV_INT::All_Ones_Vector());
        for(RANGE_ITERATOR<TV::m> it(local_range);it.Valid();it.Next())
            for(int axis=0;axis<TV::m;axis++){
                D.Append_Entry_To_Current_Row(TV::m*example.velocity.Standard_Index(vid+it.index)+axis,0);}
        D.Finish_Row();}
    D.Sort_Entries();
    int skip=0;
    for(int element=0;element<example.valid_velocity_indices.m;element++){
        TV_INT vid_cell=example.valid_velocity_cell_indices(element);
        ARRAY<ARRAY<TV,TV_INT>,TV_INT> wt=Get_Appropriate_DF<TV>(vid_cell);
        TV_INT min_corner(vid_cell/2);
        for(int axis=0;axis<TV::m;axis++)
            if(vid_cell(axis)%2==0) min_corner(axis)--;
        for(RANGE_ITERATOR<TV::m> pit(wt.domain);pit.Valid();pit.Next()){
            TV_INT pressure_index=pit.index+min_corner;
            int row=example.inv_valid_pressure_cell(pressure_index);
            if(row==-1){
                skip++;
                //LOG::printf("skip=%P\n",skip);
                continue;}
            for(RANGE_ITERATOR<TV::m> vit(wt(pit.index).domain);vit.Valid();vit.Next()){
                if(example.inv_valid_velocity_cell(vit.index+vid_cell-1)!=-1)
                for(int axis=0;axis<TV::m;axis++){
                    int col=TV::m*example.velocity.Standard_Index(vit.index+vid_cell-1)+axis;
                        D(row,col)+=wt(pit.index)(vit.index)(axis);}}}}
} 
template class MPM_KKT_KRYLOV_SYSTEM<VECTOR<float,1> >;
template class MPM_KKT_KRYLOV_SYSTEM<VECTOR<float,2> >;
template class MPM_KKT_KRYLOV_SYSTEM<VECTOR<float,3> >;
template class MPM_KKT_KRYLOV_SYSTEM<VECTOR<double,1> >;
template class MPM_KKT_KRYLOV_SYSTEM<VECTOR<double,2> >;
template class MPM_KKT_KRYLOV_SYSTEM<VECTOR<double,3> >;
}
