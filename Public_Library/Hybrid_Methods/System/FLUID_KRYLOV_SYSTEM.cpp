//#####################################################################
// Copyright 2015, Greg Klar, Andre Pradhana
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/FLUID_KRYLOV_SYSTEM.h>
#include <Hybrid_Methods/System/FLUID_KRYLOV_VECTOR.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_KRYLOV_SYSTEM<TV>::
FLUID_KRYLOV_SYSTEM(MPM_EXAMPLE<TV>& example)
    :KRYLOV_SYSTEM_BASE<T>(false,false),example(example)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_KRYLOV_SYSTEM<TV>::
~FLUID_KRYLOV_SYSTEM()
{
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void FLUID_KRYLOV_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    PHYSBAM_ASSERT(&BV!=&BF);
    const FLUID_KRYLOV_VECTOR<TV>& V=debug_cast<const FLUID_KRYLOV_VECTOR<TV>&>(BV);
    FLUID_KRYLOV_VECTOR<TV>& F=debug_cast<FLUID_KRYLOV_VECTOR<TV>&>(BF);

    for(int t=0;t<example.valid_pressure_cell_indices.m;t++){
        TV_INT index=example.valid_pressure_cell_indices(t);
        T sum((T(0)));
        for(int a=0;a<TV_INT::m;++a){
            const TV_INT axis=TV_INT::Axis_Vector(a);
            if(example.weights->Order()==1){
                if(example.cell_solid(index-axis)&&example.cell_solid(index+axis))
                    PHYSBAM_FATAL_ERROR("Fluid trapped between two solids");
                else if(example.cell_solid(index-axis)&&!example.cell_solid(index+axis))
                    sum-=2*V.p(index+axis);
                else if(!example.cell_solid(index-axis)&&example.cell_solid(index+axis))
                    sum-=2*V.p(index-axis);
                else if(!example.cell_solid(index-axis)&&!example.cell_solid(index+axis))
                    sum-=(V.p(index-axis)+V.p(index+axis));}
            else PHYSBAM_NOT_IMPLEMENTED();}
        F.p(index)=(2*TV::m*V.p(index)+sum)*example.grid.one_over_dX(0)*example.grid.one_over_dX(0);}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double FLUID_KRYLOV_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const FLUID_KRYLOV_VECTOR<TV>& X=debug_cast<const FLUID_KRYLOV_VECTOR<TV>&>(x);
    const FLUID_KRYLOV_VECTOR<TV>& Y=debug_cast<const FLUID_KRYLOV_VECTOR<TV>&>(y);
    T r=0;
#pragma omp parallel for reduction(+:r)
    for(int k=0;k<example.valid_grid_indices.m;k++){
        int i=example.valid_grid_indices(k);
        r+=X.p.array(i)*Y.p.array(i);}
    return r;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_KRYLOV_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    return sqrt(Inner_Product(BR,BR));
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void FLUID_KRYLOV_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void FLUID_KRYLOV_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void FLUID_KRYLOV_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BR) const
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void FLUID_KRYLOV_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
}
template class FLUID_KRYLOV_SYSTEM<VECTOR<float,2> >;
template class FLUID_KRYLOV_SYSTEM<VECTOR<float,3> >;
template class FLUID_KRYLOV_SYSTEM<VECTOR<double,2> >;
template class FLUID_KRYLOV_SYSTEM<VECTOR<double,3> >;
}
