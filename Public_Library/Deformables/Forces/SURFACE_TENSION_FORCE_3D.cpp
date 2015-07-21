//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE_3D.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SURFACE_TENSION_FORCE_3D<TV>::
SURFACE_TENSION_FORCE_3D(TRIANGULATED_SURFACE<T>& surface_input,T surface_tension_coefficient_input)
    :BASE(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(surface_input.particles)),surface(surface_input),surface_tension_coefficient(surface_tension_coefficient_input),dt(0),apply_explicit_forces(true),use_velocity_independent_implicit_forces(false)
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
    if(!apply_explicit_forces) return;
    for(int t=0;t<surface.mesh.elements.m;t++){
        TV_INT node=surface.mesh.elements(t);
        TV x0=surface.particles.X(node(0));
        TV x1=surface.particles.X(node(1));
        TV x2=surface.particles.X(node(2));
        TV x10=x1-x0;
        TV x21=x2-x1;
        TV x02=x0-x2;
        TV scaled_n=(T).5*surface_tension_coefficient*normals(t);
        F(node(0))+=x21.Cross(scaled_n);
        F(node(1))+=x02.Cross(scaled_n);
        F(node(2))+=x10.Cross(scaled_n);}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    areas.Resize(surface.mesh.elements.m);
    normals.Resize(surface.mesh.elements.m);
    scaled_n.Resize(surface.mesh.elements.m);
    for(int i=0;i<surface.mesh.elements.m;i++){
        TV_INT node=surface.mesh.elements(i);
        TV x0=surface.particles.X(node(0));
        TV x1=surface.particles.X(node(1));
        TV x2=surface.particles.X(node(2));
        TV x10=x1-x0;
        TV x21=x2-x1;
        TV x02=x0-x2;
        TV n=TV::Cross_Product(x02,x10);
        areas(i)=n.Normalize()/2;
        normals(i)=n;}
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
    if(!use_velocity_independent_implicit_forces) return;
    for(int t=0;t<surface.mesh.elements.m;t++){
        TV_INT node=surface.mesh.elements(t);
        TV x0=surface.particles.X(node(0));
        TV x1=surface.particles.X(node(1));
        TV x2=surface.particles.X(node(2));
        TV x10=x1-x0;
        TV x21=x2-x1;
        TV x02=x0-x2;
        TV dx0=V(node(0));
        TV dx1=V(node(1));
        TV dx2=V(node(2));
        TV normal=normals(t);
        TV scaled_n=(T).5*surface_tension_coefficient*normals(t);
        MATRIX<T,3> xx(x21,x02,x10);
        MATRIX<T,3> vv(dx0,dx1,dx2);
        TV C=-(T).5/areas(t)*xx.Transpose_Times(xx*(vv.Transpose_Times(scaled_n)));
        F(node(0))+=C(0)*normal+scaled_n.Cross(dx1-dx2);
        F(node(1))+=C(1)*normal+scaled_n.Cross(dx2-dx0);
        F(node(2))+=C(2)*normal+scaled_n.Cross(dx0-dx1);}
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
    T pe=0;
    if(apply_explicit_forces)
        for(int i=0;i<surface.mesh.elements.m;i++)
            pe+=surface_tension_coefficient*areas(i);
    return pe;
}
//#####################################################################
// Function Test_Diff
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
Test_Diff(const T time) 
{
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
    RANDOM_NUMBERS<T> random;
    for(int t=0;t<100;t++){
        T e=(T)1e-5;
        ARRAY<TV> dX(surface.particles.X.m);
        random.Fill_Uniform(dX,-e,e);
        ARRAY<TV> X2a(surface.particles.X+dX);
        ARRAY_VIEW<TV> X1(X2a);
        ARRAY<TV> F0(surface.particles.X.m);
        ARRAY<TV> F1(surface.particles.X.m);
        ARRAY<TV> G0(surface.particles.X.m);
        ARRAY<TV> G1(surface.particles.X.m);
        Update_Position_Based_State(time,true,false);
        T psi0=Potential_Energy(time);
        Add_Velocity_Independent_Forces(F0,time);
        Add_Implicit_Velocity_Independent_Forces(dX,G0,time);
        surface.particles.X.Exchange(X1);
        Update_Position_Based_State(time,true,false);
        T psi1=Potential_Energy(time);
        Add_Velocity_Independent_Forces(F1,time);
        Add_Implicit_Velocity_Independent_Forces(dX,G1,time);
        surface.particles.X.Exchange(X1);
        Update_Position_Based_State(time,true,false);
        
        ARRAY<TV> F0pF1(surface.particles.X.m);
        ARRAY<TV> F1mF0(surface.particles.X.m);
        ARRAY<TV> G0pG1o2(surface.particles.X.m);
        for(int k=0;k<F0.m;k++){
            F0pF1(k)=F0(k)+F1(k);
            G0pG1o2(k)=(G0(k)+G1(k))/2;
            F1mF0(k)=F1(k)-F0(k);}
        T df=(F0pF1.Dot(dX)/2);
        T dpsi=(psi1-psi0);
        T error=(dpsi+df)/e;
        LOG::cout<<"Energy Diff Test: "<<psi0<<" "<<psi1<<" "<<error<<std::endl;        
        T MF=sqrt(F1mF0.Magnitude_Squared());
        T MG=sqrt(G0pG1o2.Magnitude_Squared());
        T ferror=sqrt((F1mF0-G0pG1o2).Magnitude_Squared());
        LOG::cout<<"Force Diff Test: "<<MF<<" "<<MG<<" "<<ferror<<" "<<ferror/max((T)1e-30,MF,MG)<<std::endl;}
}
namespace PhysBAM{
template class SURFACE_TENSION_FORCE_3D<VECTOR<float,3> >;
template class SURFACE_TENSION_FORCE_3D<VECTOR<double,3> >;
}
