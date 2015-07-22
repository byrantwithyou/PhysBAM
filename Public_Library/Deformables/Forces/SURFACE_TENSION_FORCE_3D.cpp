//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Matrices/MATRIX.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE_3D.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SURFACE_TENSION_FORCE_3D<TV>::
SURFACE_TENSION_FORCE_3D(TRIANGULATED_SURFACE<T>& surface_input,T surface_tension_coefficient_input)
    :BASE(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(surface_input.particles)),surface(surface_input),surface_tension_coefficient(surface_tension_coefficient_input),dt(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SURFACE_TENSION_FORCE_3D<TV>::
~SURFACE_TENSION_FORCE_3D()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int t=0;t<surface.mesh.elements.m;t++)
        F.Subset(surface.mesh.elements(t))-=dE(t).Append(-dE(t).Sum());
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    typedef DIFF_LAYOUT<T,TV::m,TV::m> LAYOUT;
    E.Resize(surface.mesh.elements.m);
    dE.Resize(surface.mesh.elements.m);
    ddE.Resize(surface.mesh.elements.m);
    for(int i=0;i<surface.mesh.elements.m;i++){
        TV_INT node=surface.mesh.elements(i);
        TV x0=surface.particles.X(node(0));
        TV x1=surface.particles.X(node(1));
        TV x2=surface.particles.X(node(2));
        auto X0=Hess_From_Var<LAYOUT,0>(x0-x2);
        auto X1=Hess_From_Var<LAYOUT,1>(x1-x2);
        auto phi=surface_tension_coefficient/2*X0.Cross(X1).Magnitude();
        E(i)=phi.x;
        Extract<0>(dE(i),phi.dx);
        Extract<0,0>(ddE(i),phi.ddx);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(int t=0;t<surface.mesh.elements.m;t++){
        TV_INT node=surface.mesh.elements(t);
        VECTOR<TV,2> v(V(node(0))-V(node(2)),V(node(1))-V(node(2))),f;
        const MATRIX<MATRIX<T,TV::m>,2>& M=ddE(t);
        for(int i=0;i<2;i++)
            for(int j=0;j<2;j++)
                f(i)+=M(i,j)*v(j);
        F.Subset(node)-=f.Append(-f.Sum());}
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int SURFACE_TENSION_FORCE_3D<TV>::
Velocity_Dependent_Forces_Size() const
{
    return 0;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR SURFACE_TENSION_FORCE_3D<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    surface.mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR SURFACE_TENSION_FORCE_3D<TV>::
Potential_Energy(const T time) const
{
    return E.Sum();
}
namespace PhysBAM{
template class SURFACE_TENSION_FORCE_3D<VECTOR<float,3> >;
template class SURFACE_TENSION_FORCE_3D<VECTOR<double,3> >;
}
