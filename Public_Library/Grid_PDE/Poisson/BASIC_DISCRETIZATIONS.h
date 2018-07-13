//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIC_DISCRETIZATIONS
//#####################################################################
#ifndef __BASIC_DISCRETIZATIONS__
#define __BASIC_DISCRETIZATIONS__

#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
namespace PhysBAM{

// int face_func(const FACE_INDEX<TV::m>& face, T& rhs)
//   - if regular face, returns dof index of face
//   - if Neumann face, return -1
//     if bc_times_matrix_transpose!=0, sets rhs to a bc value
// int cell_func(const VECTOR<int,TV::m>& cell, T& rhs)
//   - if regular cell, returns dof index of cell
//   - if Dirichlet cell, return -1
//     if bc_times_matrix!=0, sets rhs to a bc value
template<class T,class TV,class FACE_FUNC,class CELL_FUNC>
void Compute_Gradient_Matrix(SYSTEM_MATRIX_HELPER<T>& M,const GRID<TV>& grid,
    int ghost,FACE_FUNC face_func,CELL_FUNC cell_func,
    ARRAY<T>* bc_times_matrix,ARRAY<T>* bc_times_matrix_transpose)
{
    for(FACE_RANGE_ITERATOR<TV::m> it(grid.Domain_Indices(ghost));it.Valid();it.Next()){
        T rhs_p0=0,rhs_p1=0,rhs_u=0,e=grid.one_over_dX(it.face.axis);
        int f=face_func(it.face,rhs_u);
        int c0=cell_func(it.face.First_Cell_Index(),rhs_p0);
        int c1=cell_func(it.face.Second_Cell_Index(),rhs_p1);
        if(f==-1){
            if(bc_times_matrix_transpose && (c0>=0 || c1>=0)){
                T rhs=e*rhs_u;
                if(c0>=0) (*bc_times_matrix_transpose)(c0)-=rhs;
                if(c1>=0) (*bc_times_matrix_transpose)(c1)+=rhs;}
            continue;}

        if(c0>=0) M.data.Append({f,c0,-e});
        else if(bc_times_matrix) (*bc_times_matrix)(f)-=e*rhs_p0;
        
        if(c1>=0) M.data.Append({f,c1,e});
        else if(bc_times_matrix) (*bc_times_matrix)(f)+=e*rhs_p1;}
}

// Computes negative Poisson matrix (SPSD)
// int face_func(const FACE_INDEX<TV::m>& face, T& beta, T& rhs)
//   - if regular face, returns dof index of face
//     beta is a scaling factor
//   - if Neumann face, return -1
//     if bc_times_matrix!=0, sets rhs to a bc value
// int cell_func(const VECTOR<int,TV::m>& cell, T& rhs)
//   - if regular cell, returns dof index of cell
//   - if Dirichlet cell, return -1
//     if bc_times_matrix!=0, sets rhs to a bc value
template<class T,class TV,class FACE_FUNC,class CELL_FUNC>
void Compute_Poisson_Matrix(SYSTEM_MATRIX_HELPER<T>& M,const GRID<TV>& grid,
    int ghost,FACE_FUNC face_func,CELL_FUNC cell_func,
    ARRAY<T>* bc_times_matrix)
{
    TV dX=grid.one_over_dX,dX2=dX*dX;
    for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()){
        T diag=0,bc=0;
        int c=cell_func(it.index,bc);
        if(c<0) continue;

        for(int a=0;a<TV::m;a++){
            for(int s=0,sn=-1;s<2;s++,sn+=2){
                FACE_INDEX<TV::m> face(a,it.index);
                face.index(a)+=s;
                T bc=0,beta=1;
                int f=face_func(face,beta,bc);
                if(f<0){
                    if(bc_times_matrix) (*bc_times_matrix)(c)+=sn*bc*dX(a);
                    continue;}

                VECTOR<int,TV::m> cell=it.index;
                cell(a)+=sn;
                bc=0;
                int b=cell_func(cell,bc);

                T e=dX2(a)*beta;
                diag+=e;
                if(b>=0) M.data.Append({c,b,-e});
                else if(bc_times_matrix) (*bc_times_matrix)(c)=e*bc;}}
        M.data.Append({c,c,diag});}
}

// Computes Laplacian matrices on faces (SPSD)
// Boundary conditions are relative to cells
// int face_func(const FACE_INDEX<TV::m>& face, T& beta, T& rhs)
//   - if regular face, returns dof index of face
//     beta is a scaling factor
//   - if Neumann face, return -1
//     if bc_times_matrix!=0, sets rhs to a bc value
// int cell_func(const VECTOR<int,TV::m>& cell, T& rhs)
//   - if regular cell, returns dof index of cell
//   - if Dirichlet cell, return -1
//     if bc_times_matrix!=0, sets rhs to a bc value
template<class T,class TV,class FACE_FUNC,class CELL_FUNC>
void Compute_Vector_Laplacian_Matrix(SYSTEM_MATRIX_HELPER<T>& M,const GRID<TV>& grid,
    int ghost,FACE_FUNC face_func,CELL_FUNC cell_func,
    ARRAY<T>* bc_times_matrix)
{
    typedef VECTOR<int,TV::m> TV_INT;
    auto p_psi_D=[&cell_func](const TV_INT& c){T rhs;return cell_func(c,rhs)<0;};

    for(int axis=0;axis<TV::m;axis++){
        // Compute as needed rather than store.  face_grid index.
        auto axis_face_func=[axis,p_psi_D](const FACE_INDEX<TV::m>& face, T& beta, T& rhs)
            {
                beta=1;
                rhs=0;
                TV_INT a=face.index,b=a;
                b(axis)--;
                if(face.axis==axis) return p_psi_D(b)?-1:0;
                if(p_psi_D(a) && p_psi_D(b)) return -1;
                a(face.axis)--;
                b(face.axis)--;
                if(p_psi_D(a) && p_psi_D(b)) return -1;
                return 0;
            };
        auto axis_cell_func=[&face_func,axis](const TV_INT& cell, T& rhs)
            {
                T beta;
                return face_func(FACE_INDEX<TV::m>(axis,cell),beta,rhs);
            };
        Compute_Poisson_Matrix(M,grid.Get_Face_Grid(axis),ghost,axis_face_func,axis_cell_func,bc_times_matrix);
    }
}
}
#endif
