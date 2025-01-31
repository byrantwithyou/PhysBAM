//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/MATRIX.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Forces/LINEAR_POINT_ATTRACTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LINEAR_POINT_ATTRACTION<TV>::
LINEAR_POINT_ATTRACTION(T_MESH& mesh,const TV& pt,T coefficient_input)
    :BASE(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(mesh.particles)),surface(mesh),coefficient(coefficient_input),point(pt),dt(0),apply_explicit_forces(true),apply_implicit_forces(true)
{
    Get_Unique(referenced_particles,mesh.mesh.elements.Flattened());
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LINEAR_POINT_ATTRACTION<TV>::
~LINEAR_POINT_ATTRACTION()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(apply_explicit_forces)
        for(int i=0;i<referenced_particles.m;i++){int p=referenced_particles(i);
            F(p)+=coefficient*(point-surface.particles.X(p));}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    if(apply_implicit_forces)
        for(int i=0;i<referenced_particles.m;i++){int p=referenced_particles(i);
            F(p)-=dt*coefficient*V(p);}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose) const
{
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int LINEAR_POINT_ATTRACTION<TV>::
Velocity_Dependent_Forces_Size() const
{
    int size=0;
    if(apply_implicit_forces) size+=referenced_particles.m*TV::m;
    return size;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    T c=sqrt(dt*coefficient);
    if(apply_implicit_forces)
        for(int i=0;i<referenced_particles.m;i++){int p=referenced_particles(i);
            TV t=c*V(p);
            for(int a=0;a<TV::m;a++) aggregate(i*TV::m+a)=t(a);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    T c=sqrt(dt*coefficient);
    if(apply_implicit_forces)
        for(int i=0;i<referenced_particles.m;i++){int p=referenced_particles(i);
            TV t;
            for(int a=0;a<TV::m;a++) t(a)=aggregate(i*TV::m+a);
            F(p)+=c*t;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_POINT_ATTRACTION<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    surface.mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void LINEAR_POINT_ATTRACTION<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_POINT_ATTRACTION<TV>::
Potential_Energy(const T time) const
{
    T pe=0;
    if(apply_explicit_forces)
        for(int i=0;i<referenced_particles.m;i++){int p=referenced_particles(i);
            pe+=coefficient/2*(surface.particles.X(p)-point).Magnitude_Squared();}
    return pe;
}
namespace PhysBAM{
template class LINEAR_POINT_ATTRACTION<VECTOR<float,2> >;
template class LINEAR_POINT_ATTRACTION<VECTOR<double,2> >;
}
