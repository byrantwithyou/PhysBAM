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
#include "FLUID_LAYOUT.h"
namespace PhysBAM{

template<class T,class TV>
void Compute_Full_Matrix(const GRID<TV>& grid,SYSTEM_MATRIX_HELPER<T>& MH,
    ARRAY<T>& rhs_vector,const FLUID_LAYOUT<TV>& fl,T mu)
{
    typedef VECTOR<int,TV::m> TV_INT;

    rhs_vector.Resize(fl.Total_Dofs(),init_all,0);
    auto face_func=[&fl](const FACE_INDEX<TV::m>& face, T& rhs)
        {
            auto& uf=fl.used_faces(face);
            if(uf.type==fluid) return uf.global_id;
            rhs=uf.bc_value;
            return -1;
        };
    auto face_func_beta=[mu,face_func](const FACE_INDEX<TV::m>& face, T& beta, T& rhs)
        {
            beta=mu;
            return face_func(face,rhs);
        };
    auto cell_func=[&fl](const VECTOR<int,TV::m>& cell, T& rhs)
        {
            auto& uc=fl.used_cells(cell);
            if(uc.type==fluid) return uc.global_id;
            rhs=uc.bc_value;
            return -1;
        };

    Compute_Vector_Laplacian_Matrix(MH,grid,0,face_func_beta,cell_func,&rhs_vector);
    MH.New_Block();
    Compute_Gradient_Matrix(MH,grid,0,face_func,cell_func,&rhs_vector,&rhs_vector);
    MH.Add_Transpose();
    MH.Compact(fl.Total_Dofs());
    rhs_vector=-rhs_vector;
}

template<class T,class TV>
void Solve_And_Display_Solution(const GRID<TV>& grid,const FLUID_LAYOUT<TV>& fl,
    const SYSTEM_MATRIX_HELPER<T>& MH,const ARRAY<T>& rhs_vector,ARRAY<T>* sol_out)
{
    SPARSE_MATRIX_FLAT_MXN<T> M;
    MH.Set_Matrix(fl.Total_Dofs(),fl.Total_Dofs(),M);
    
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > KRY_VEC;
    typedef MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<T>,T,KRY_VEC> KRY_MAT;
    KRY_MAT sys(M);
    KRY_VEC rhs,sol;
    rhs.v=rhs_vector;
    sol.v.Resize(fl.Total_Dofs());

    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    sys.Test_System(sol);
    OCTAVE_OUTPUT<T>("M.txt").Write("M",sys,rhs);
    OCTAVE_OUTPUT<T>("b.txt").Write("b",rhs);

    MINRES<T> mr;
    bool converged=mr.Solve(sys,sol,rhs,av,1e-6,0,1000);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

    OCTAVE_OUTPUT<T>("x.txt").Write("x",sol);
    if(sol_out) *sol_out=sol.v;
    
    ARRAY<T,FACE_INDEX<TV::m> > face_velocity(grid,1);
    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& uf=fl.used_faces(it.Full_Index());
        if(uf.type!=fluid) continue;
        face_velocity(it.Full_Index())=sol.v(uf.global_id);
    }
    Flush_Frame(face_velocity,"solve");
}

template void Compute_Full_Matrix<double,VECTOR<double,2> >(
    GRID<VECTOR<double,2> > const&,SYSTEM_MATRIX_HELPER<double>&,
    ARRAY<double,int>&,FLUID_LAYOUT<VECTOR<double,2> > const&,double);
template void Solve_And_Display_Solution<double,VECTOR<double,2> >(
    GRID<VECTOR<double,2> > const&,FLUID_LAYOUT<VECTOR<double,2> > const&,
    SYSTEM_MATRIX_HELPER<double> const&,ARRAY<double,int> const&,ARRAY<double,int>*);
}
