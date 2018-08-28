//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Poisson/BASIC_DISCRETIZATIONS.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <fstream>
#include <map>
#include <string>
#include "COMMON.h"
#include "FLUID_LAYOUT.h"
namespace PhysBAM{

template<class T,class TV>
void Compute_Full_Matrix(const GRID<TV>& grid,ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    ARRAY<T,CODE_ID>& code_values,ARRAY<T,DOF_ID>& rhs_vector,const FLUID_LAYOUT<TV>& fl,T mu)
{
    typedef VECTOR<int,TV::m> TV_INT;

    rhs_vector.Resize(fl.Total_Dofs(),init_all,0);

    TV de=mu*grid.one_over_dX*grid.one_over_dX;
    for(auto t:grid.one_over_dX) code_values.Append(t);
    for(auto t:-grid.one_over_dX) code_values.Append(t);
    for(auto t:-de) code_values.Append(t);
    for(int i=0;i<(1<<2*TV::m);i++){
        T x=0;
        for(int a=0;a<TV::m;a++)
            x+=((i>>2*a)&3)*de(a);
        code_values.Append(x);}
    
    for(FACE_RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()){
        int diag_code=3*TV::m;
        auto& uf=fl.used_faces(it.face);
        if(uf.type!=fluid) continue;

        for(int a=0;a<TV::m;a++){
            for(int s=0,sn=-1;s<2;s++,sn+=2){
                FACE_INDEX<TV::m> n(it.face);
                n.index(a)+=sn;
                auto& un=fl.used_faces(n);

                if(un.type==nodof){
                    if(it.face.axis==a) continue;
                    FACE_INDEX<TV::m> g(it.face);
                    g.index(a)+=s;
                    if(fl.used_faces(g).type!=wall){
                        FACE_INDEX<TV::m> h(g);
                        h.index(it.face.axis)--;
                        if(fl.used_faces(h).type!=wall)
                            continue;}}

                FACE_INDEX<TV::m> face=it.face;
                face.index(a)+=sn;
                diag_code+=1<<2*a;
                if(un.type==fluid) coded_entries.Append({uf.global_id,un.global_id,CODE_ID(2*TV::m+a)});
                else if(un.type==wall) rhs_vector(uf.global_id)=de(a)*un.bc_value;}}
        coded_entries.Append({uf.global_id,uf.global_id,CODE_ID(diag_code)});}

    for(FACE_RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()){
        auto& uf=fl.used_faces(it.face);
        auto& uc0=fl.used_cells(it.face.First_Cell_Index());
        auto& uc1=fl.used_cells(it.face.Second_Cell_Index());
        T e=grid.one_over_dX(it.face.axis);

        if(uf.type!=fluid){
            if(uf.type!=wall) continue;
            if(uc0.type==fluid) rhs_vector(uc0.global_id)-=e*uf.bc_value;
            if(uc1.type==fluid) rhs_vector(uc1.global_id)+=e*uf.bc_value;
            continue;}

        if(uc0.type==fluid){
            coded_entries.Append({uf.global_id,uc0.global_id,CODE_ID(TV::m+it.face.axis)});
            coded_entries.Append({uc0.global_id,uf.global_id,CODE_ID(TV::m+it.face.axis)});}
        else if(uc0.type==dirichlet) rhs_vector(uf.global_id)-=e*uc0.bc_value;
        
        if(uc1.type==fluid){
            coded_entries.Append({uf.global_id,uc1.global_id,CODE_ID(it.face.axis)});
            coded_entries.Append({uc1.global_id,uf.global_id,CODE_ID(it.face.axis)});}
        else if(uc1.type==dirichlet) rhs_vector(uf.global_id)+=e*uc1.bc_value;}

    rhs_vector=-rhs_vector;
}

template<class T,class TV>
void Solve_And_Display_Solution(const GRID<TV>& grid,const FLUID_LAYOUT<TV>& fl,
    const SYSTEM_MATRIX_HELPER<T>& MH,const ARRAY<T,DOF_ID>& rhs_vector,
    ARRAY<T,DOF_ID>* sol_out)
{
    SPARSE_MATRIX_FLAT_MXN<T> M;
    MH.Set_Matrix(Value(fl.Total_Dofs()),Value(fl.Total_Dofs()),M);
    
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > KRY_VEC;
    typedef MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<T>,T,KRY_VEC> KRY_MAT;
    KRY_MAT sys(M);
    KRY_VEC rhs,sol;
    rhs.v=reinterpret_cast<const ARRAY<T>&>(rhs_vector);
    sol.v.Resize(Value(fl.Total_Dofs()));

    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    sys.Test_System(sol);
    OCTAVE_OUTPUT<T>("M.txt").Write("M",sys,rhs);
    OCTAVE_OUTPUT<T>("b.txt").Write("b",rhs);

    MINRES<T> mr;
    bool converged=mr.Solve(sys,sol,rhs,av,1e-8,0,100000);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

    OCTAVE_OUTPUT<T>("x.txt").Write("x",sol);
    if(sol_out) reinterpret_cast<ARRAY<T>&>(*sol_out)=sol.v;
    
    ARRAY<T,FACE_INDEX<TV::m> > face_velocity(grid,1);
    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& uf=fl.used_faces(it.Full_Index());
        if(uf.type!=fluid) continue;
        face_velocity(it.Full_Index())=sol.v(Value(uf.global_id));
    }
    Flush_Frame(face_velocity,"solve");
}
template void Compute_Full_Matrix<double,VECTOR<double,2> >(GRID<VECTOR<double,2> > const&,ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID>,int>&,ARRAY<double,CODE_ID>&,ARRAY<double,DOF_ID>&,FLUID_LAYOUT<VECTOR<double,2> > const&,double);
template void Compute_Full_Matrix<double,VECTOR<double,3> >(GRID<VECTOR<double,3> > const&,ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID>,int>&,ARRAY<double,CODE_ID>&,ARRAY<double,DOF_ID>&,FLUID_LAYOUT<VECTOR<double,3> > const&,double);
template void Solve_And_Display_Solution<double,VECTOR<double,2> >(GRID<VECTOR<double,2> > const&,FLUID_LAYOUT<VECTOR<double,2> > const&,SYSTEM_MATRIX_HELPER<double> const&,ARRAY<double,DOF_ID> const&,ARRAY<double,DOF_ID>*);
template void Solve_And_Display_Solution<double,VECTOR<double,3> >(GRID<VECTOR<double,3> > const&,FLUID_LAYOUT<VECTOR<double,3> > const&,SYSTEM_MATRIX_HELPER<double> const&,ARRAY<double,DOF_ID> const&,ARRAY<double,DOF_ID>*);
}
