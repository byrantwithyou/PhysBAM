//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLIP_KRYLOV_SYSTEM__
#define __FLIP_KRYLOV_SYSTEM__

#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>

#include "FLIP_KRYLOV_VECTOR.h"

namespace PhysBAM{

template<class TV>
class FLIP_KRYLOV_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef ARRAY<short int,TV_INT> INT_CELL_ARRAY;
    typedef ARRAY<T,TV_INT> T_CELL_ARRAY;
    typedef GRID<TV> T_GRID;
    typedef FACE_INDEX<TV::m> T_FACE;
    typedef ARRAY<T,T_FACE> T_FACE_ARRAY;

    enum CELL_TYPE{
        CELL_NEUMANN=0,
        CELL_DIRICHLET=1,
        CELL_INTERIOR=2};

public:

    const T one_over_dx_squared;
    const T threads;
    const ARRAY<TV_INT>& interior_cells;
    const INT_CELL_ARRAY& cell_type;
    const T_FACE_ARRAY& mass;

    T dt;

    FLIP_KRYLOV_SYSTEM(const T one_over_dx,const int threads_input,const ARRAY<TV_INT>& interior_cells_input,const INT_CELL_ARRAY& cell_type_input,const T_FACE_ARRAY& mass_input)
        :BASE(true,false),one_over_dx_squared(sqr(one_over_dx)),threads(threads_input),interior_cells(interior_cells_input),cell_type(cell_type_input),mass(mass_input)
    {}

    virtual ~FLIP_KRYLOV_SYSTEM(){}

    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const PHYSBAM_OVERRIDE
    {
        const FLIP_KRYLOV_VECTOR<TV>& xx=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(x);
        FLIP_KRYLOV_VECTOR<TV>& r=debug_cast<FLIP_KRYLOV_VECTOR<TV>&>(result);

#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            r.p(index)=T();
            for(int j=0;j<TV::m;j++)
            for(int s=-1;s<=1;s+=2){
                const TV_INT offset_index=index+TV_INT::Axis_Vector(j)*s;
                T mass_scale;
                if(s>0) mass_scale=1/mass(T_FACE(j,offset_index));
                else mass_scale=1/mass(T_FACE(j,index));
                if(cell_type(offset_index)==CELL_INTERIOR)
                    r.p(index)+=(xx.p(offset_index)-xx.p(index))*mass_scale;
                else if (cell_type(offset_index)==CELL_DIRICHLET)
                    r.p(index)+=-xx.p(index)*mass_scale;}
            r.p(index)*=-dt*one_over_dx_squared;}
    }

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const PHYSBAM_OVERRIDE
    {
        const FLIP_KRYLOV_VECTOR<TV>& xx=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(x);
        const FLIP_KRYLOV_VECTOR<TV>& yy=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(y);

        ARRAY<double> result_per_thread(threads);
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            const int tid=omp_get_thread_num();
            result_per_thread(tid)+=(double)xx.p(index)*(double)yy.p(index);}
        double result=0;
        for(int tid=0;tid<threads;tid++)
            result+=result_per_thread(tid);
        return result;
    }

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
    {
        const FLIP_KRYLOV_VECTOR<TV>& xx=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(x);

        ARRAY<T> result_per_thread(threads);
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            const int tid=omp_get_thread_num();
            T& r=result_per_thread(tid);
            r=max(r,(T)fabs(xx.p(index)));}
        T result=0;
        for(int tid=0;tid<threads;tid++)
            result=max(result,result_per_thread(tid));
        return result;
    }

    void Project(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE{}
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE {}
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE {}

protected:

    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const PHYSBAM_OVERRIDE
    {
        const FLIP_KRYLOV_VECTOR<TV>& xx=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(x);
        FLIP_KRYLOV_VECTOR<TV>& r=debug_cast<FLIP_KRYLOV_VECTOR<TV>&>(result);
        
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            T scale=0;
            for(int j=0;j<TV::m;j++)
            for(int s=-1;s<=1;s+=2){
                const TV_INT offset_index=index+TV_INT::Axis_Vector(j)*s;
                T mass_scale;
                if(s>0) mass_scale=1/mass(T_FACE(j,offset_index));
                else mass_scale=1/mass(T_FACE(j,index));
                if(cell_type(offset_index)==CELL_INTERIOR || cell_type(offset_index)==CELL_DIRICHLET)
                    scale+=-mass_scale;}
            scale*=-dt*one_over_dx_squared;
            r.p(index)=xx.p(index)/scale;}
    }
//#####################################################################
};
}
#endif
