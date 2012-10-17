//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHING_CUBES_SYSTEM
//#####################################################################
#ifndef __MARCHING_CUBES_SYSTEM__
#define __MARCHING_CUBES_SYSTEM__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>

namespace PhysBAM{
//#####################################################################
// Class MARCHING_CUBES_VECTOR
//#####################################################################
template<class TV>
class MARCHING_CUBES_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;

public:

    ARRAY<TV> x; 
    ARRAY<TV> n; 

    MARCHING_CUBES_VECTOR(){}
    virtual ~MARCHING_CUBES_VECTOR(){}

    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE{
        x+=debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x;
        n+=debug_cast<const MARCHING_CUBES_VECTOR&>(bv).n;
        return *this;}

    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE{
        x-=debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x;
        n-=debug_cast<const MARCHING_CUBES_VECTOR&>(bv).n;
        return *this;}

    BASE& operator*=(const T a) PHYSBAM_OVERRIDE{
        x*=a; n*=a; return *this;}

    void Copy(const T c1,const BASE& bv1) PHYSBAM_OVERRIDE{
        x=debug_cast<const MARCHING_CUBES_VECTOR&>(bv1).x*c1;
        n=debug_cast<const MARCHING_CUBES_VECTOR&>(bv1).n*c1;}

    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE{
        x=debug_cast<const MARCHING_CUBES_VECTOR&>(bv1).x*c1+debug_cast<const MARCHING_CUBES_VECTOR&>(bv2).x;
        n=debug_cast<const MARCHING_CUBES_VECTOR&>(bv1).n*c1+debug_cast<const MARCHING_CUBES_VECTOR&>(bv2).n;}

    int Raw_Size() const PHYSBAM_OVERRIDE
    {return (x.m+n.m)*TV::m;}

    T& Raw_Get(int i) PHYSBAM_OVERRIDE{
        if(i<x.m*TV::m) return x(i/TV::m)(i%TV::m);
        else{i-=x.m*TV::m;return n(i/TV::m)(i%TV::m);}}

    KRYLOV_VECTOR_BASE<T>* Clone_Default() const PHYSBAM_OVERRIDE{
        MARCHING_CUBES_VECTOR* V=new MARCHING_CUBES_VECTOR;
        V->x.Resize(x.m); V->n.Resize(n.m); return V;}

    void Resize(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE{
        x.Resize(debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x.m);
        n.Resize(debug_cast<const MARCHING_CUBES_VECTOR&>(bv).n.m);}
};
//#####################################################################
// Class MARCHING_CUBES_SYSTEM
//#####################################################################
template<class TV>
class MARCHING_CUBES_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef MARCHING_CUBES_VECTOR<TV> VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;

public:
    struct BLOCK
    {
        ARRAY<int> index;
        MATRIX_MXN<SYMMETRIC_MATRIX<T,TV::m> > matrix_xx;
        ARRAY<MATRIX<T,TV::m> > matrix_xn;
        SYMMETRIC_MATRIX<T,TV::m> matrix_nn;
    };
    
    ARRAY<BLOCK> blocks;
    ARRAY<int> project_flags;

    MARCHING_CUBES_SYSTEM():BASE(false,false){}
    virtual ~MARCHING_CUBES_SYSTEM(){}

//#####################################################################
    void Multiply(const KRYLOV_VECTOR_BASE<T>& bv_input,KRYLOV_VECTOR_BASE<T>& bv_result) const PHYSBAM_OVERRIDE
    {
        const ARRAY<TV>& x_input=debug_cast<const VECTOR_T&>(bv_input).x;
        const ARRAY<TV>& n_input=debug_cast<const VECTOR_T&>(bv_input).n;
        ARRAY<TV>& x_result=debug_cast<VECTOR_T&>(bv_result).x;
        ARRAY<TV>& n_result=debug_cast<VECTOR_T&>(bv_result).n;
        x_result.Fill(TV());
        n_result.Fill(TV());
        for(int b=0;b<blocks.m;b++){
            const BLOCK& block=blocks(b);
            for(int i=0;i<block.index.m;i++){
                const int index_i=block.index(i);
                for(int j=0;j<block.index.m;j++){
                    const int index_j=block.index(j);
                    x_result(index_i)+=block.matrix_xx(i,j)*x_input(index_j);}
                x_result(index_i)+=block.matrix_xn(i)*n_input(b);
                n_result(b)+=block.matrix_xn(i).Transpose_Times(x_input(index_i));}
            n_result(b)+=block.matrix_nn*n_input(b);}
    }

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) const PHYSBAM_OVERRIDE{
        const VECTOR_T& v1=debug_cast<const MARCHING_CUBES_VECTOR<TV>&>(bv1);
        const VECTOR_T& v2=debug_cast<const MARCHING_CUBES_VECTOR<TV>&>(bv2);
        return v1.x.Dot(v2.x)+v1.n.Dot(v2.n);}

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE
    {return sqrt(Inner_Product(bv,bv));}

    void Project(KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE{
        ARRAY<TV>& x=debug_cast<VECTOR_T&>(bv).x;
        for(int i=0;i<project_flags.m;i++)
            for(int j=0;j<TV::m;j++)
                if(!(project_flags(i)&(1<<j)))
                    x(i)(j)=0;}

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE {}
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE {Project(bv);}
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& br,KRYLOV_VECTOR_BASE<T>& bz) const PHYSBAM_OVERRIDE {}
//#####################################################################
    T Set_Matrix_Block_And_Rhs(const ARRAY_VIEW<TV>& particles,const ARRAY<int>& particle_indices,
        const ARRAY<int>& index_map,const ARRAY<int>& reverse_index_map,const TV& n,ARRAY<TV>& positions_rhs_full,TV& rhs_n)
    {
        // NOTE: assumes n is normalized !!!

        blocks.Add_End();
        BLOCK& block=blocks.Last();
        
        const int x_full_count=particle_indices.m;
        
        // compute average and initialize reduced index list
        TV average;
        for(int i=0;i<particle_indices.m;i++){
            const int index=particle_indices(i);
            const int reduced_index=reverse_index_map(index);
            if(reduced_index>=0) block.index.Append(reduced_index);
            average+=particles(index);}
        average/=x_full_count;
        
        const int x_dofs=block.index.m;
        
        // normal differential
        const SYMMETRIC_MATRIX<T,TV::m> nnt=SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(n);
        const SYMMETRIC_MATRIX<T,TV::m> dn=SYMMETRIC_MATRIX<T,TV::m>::Identity_Matrix()-nnt;
        
        T E=0;
        
        // set RHS
        const T average_dot_n=average.Dot(n);
        for(int i=0;i<particle_indices.m;i++){
            const int index=particle_indices(i);
            const TV& x=particles(index);
            const T x_dot_n=x.Dot(n);
            positions_rhs_full(index)+=n*(2*(x_dot_n-average_dot_n));
            rhs_n+=(x-average)*2*x_dot_n;
            E+=sqr(x_dot_n);}
        rhs_n=dn*rhs_n;
        E-=sqr(average_dot_n)*x_full_count;
        
        // set matrix block
        block.matrix_xx.Resize(x_dofs,x_dofs);
        block.matrix_xn.Resize(x_dofs);
        
        for(int i=0;i<x_dofs;i++){
            for(int j=0;j<x_dofs;j++){
                block.matrix_xx(i,j)=nnt;
                block.matrix_xx(i,j)*=2*((i==j)-(T)1/particle_indices.m);}
            const int index=index_map(block.index(i));
            const TV& x=particles(index);
            block.matrix_xn(i)+=dn*(2*(x.Dot(n)-average_dot_n));
            block.matrix_xn(i)+=MATRIX<T,TV::m>::Outer_Product(n,dn*(x-average))*2;}
        
        for(int i=0;i<particle_indices.m;i++){
            const int index=particle_indices(i);
            const TV& x=particles(index);
            const TV xma=x-average;
            SYMMETRIC_MATRIX<T,TV::m> tmp=
                -MATRIX<T,TV::m>::Outer_Product(n,xma).Symmetric_Part()*2
                +(nnt*(T)3-SYMMETRIC_MATRIX<T,TV::m>::Identity_Matrix())*n.Dot(xma);
            block.matrix_nn+=tmp*(x.Dot(n));
            block.matrix_nn+=SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(dn*x);}
        block.matrix_nn-=SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(dn*average)*x_full_count;
        block.matrix_nn*=2;
        
        return E;
    }
    
    T Setup(MARCHING_CUBES_VECTOR<TV>& rhs,MARCHING_CUBES_VECTOR<TV>& sol,const int normals_count,const ARRAY<int>& position_dofs,
        const ARRAY<int>& index_map,const ARRAY<int>& reverse_index_map,const ARRAY<TV_INT>& junction_cells,
        const HASHTABLE<TV_INT,HASHTABLE<VECTOR<int,2>,VECTOR<int,2> > >& cell_to_element,
        const HASHTABLE<VECTOR<int,2>,T_SURFACE*>& surface,const GEOMETRY_PARTICLES<TV>& particles,const ARRAY<TV>& normals)
    {
        T E=0;
        ARRAY<TV> positions_rhs_full(particles.number);
        rhs.n.Resize(normals_count);
        project_flags=position_dofs.Subset(index_map);
        int normal_index=-1;
        for(int jc=0;jc<junction_cells.m;jc++){
            const TV_INT& junction_cell_index=junction_cells(jc);
            const HASHTABLE<VECTOR<int,2>,VECTOR<int,2> >& junction_cell_elements=cell_to_element.Get(junction_cell_index);
            for(typename HASHTABLE<VECTOR<int,2>,VECTOR<int,2> >::CONST_ITERATOR it(junction_cell_elements);it.Valid();it.Next()){
                normal_index++;
                HASHTABLE<int> particle_indices_ht;
                const VECTOR<int,2>& color_pair=it.Key();
                const T_SURFACE& color_pair_surface=*surface.Get(color_pair);
                const RANGE<TV_INT> range(RANGE<TV_INT>::Centered_Box()*2);
                for(RANGE_ITERATOR<TV::m> it2(range);it2.Valid();it2.Next()){
                    const TV_INT& cell_index=junction_cell_index+it2.index;
                    if(cell_to_element.Contains(cell_index)){
                        const HASHTABLE<VECTOR<int,2>,VECTOR<int,2> >& cell_elements=cell_to_element.Get(cell_index);
                        VECTOR<int,2> elements_range;
                        if(cell_elements.Get(color_pair,elements_range))
                            for(int i=elements_range.x;i<elements_range.y;i++){
                                TV_INT element=color_pair_surface.mesh.elements(i);
                                for(int j=0;j<TV::m;j++)
                                    if(!particle_indices_ht.Contains(element(j)))
                                        particle_indices_ht.Insert(element(j));}}}
                ARRAY<int> particle_indices;
                for(typename HASHTABLE<int>::ITERATOR it(particle_indices_ht);it.Valid();it.Next())particle_indices.Append(it.Key());
                E+=Set_Matrix_Block_And_Rhs(particles.X,particle_indices,index_map,reverse_index_map,
                    normals(normal_index),positions_rhs_full,rhs.n(normal_index));}}
        rhs.x=positions_rhs_full.Subset(index_map);
        sol.Resize(rhs);
        return E;
    }
};
}
#endif
