//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SURFACE_TENSION_FORCE<TV>::
SURFACE_TENSION_FORCE(SEGMENTED_CURVE_2D<T>& surface_input,T surface_tension_coefficient_input)
    :BASE(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(surface_input.particles)),surface(surface_input),surface_tension_coefficient(surface_tension_coefficient_input),dt(0),apply_explicit_forces(true),apply_implicit_forces(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SURFACE_TENSION_FORCE<TV>::
~SURFACE_TENSION_FORCE()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(apply_explicit_forces)
        for(int i=0;i<surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            TV f=(surface.particles.X(k.x)-surface.particles.X(k.y))*coefficients(i);
            F(k.x)-=f;
            F(k.y)+=f;}
}
//#####################################################################
// Function Tangential_Helper
//#####################################################################
template<class T> static void Tangential_Helper(MATRIX<T,3,2>& tangential,const VECTOR<T,3>& normal)
{
    tangential.Set_Column(1,normal.Unit_Orthogonal_Vector());
    tangential.Set_Column(2,normal.Cross_Product(normal,tangential.Column(1)));
}
template<class T> static void Tangential_Helper(MATRIX<T,2,1>& tangential,const VECTOR<T,2>& normal)
{
    tangential.Set_Column(1,normal.Orthogonal_Vector());
}
template<class T> static void Tangential_Helper(MATRIX<T,1,0>& tangential,const VECTOR<T,1>& normal)
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    coefficients.Resize(surface.mesh.elements.m);
    normal.Resize(surface.mesh.elements.m);
    for(int i=0;i<surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
        normal(i)=(surface.particles.X(k.x)-surface.particles.X(k.y)).Orthogonal_Vector();
        coefficients(i)=surface_tension_coefficient/normal(i).Normalize();}
    if(this->compute_half_forces){
        sqrt_coefficients.Resize(surface.mesh.elements.m);
        for(int i=0;i<surface.mesh.elements.m;i++)
            sqrt_coefficients(i)=sqrt(dt*coefficients(i));}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    // this is different from Add_Velocity_Dependent_Forces because there is no dt in the formular
    for(int i=0;i<surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
        TV f=(TV::Dot_Product((V(k.x)-V(k.y)),normal(i))*coefficients(i))*normal(i);
        F(k.x)-=f;
        F(k.y)+=f;}
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int SURFACE_TENSION_FORCE<TV>::
Velocity_Dependent_Forces_Size() const
{
    int size=0;
    if(apply_implicit_forces) size+=surface.mesh.elements.m;
    return size;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    if(apply_implicit_forces)
        for(int i=0;i<surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            aggregate(i)=TV::Dot_Product((V(k.x)-V(k.y)),normal(i))*sqrt_coefficients(i);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    if(apply_implicit_forces)
        for(int i=0;i<surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            TV f=(aggregate(i)*sqrt_coefficients(i))*normal(i);
            F(k.x)+=f;
            F(k.y)-=f;}
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    if(apply_implicit_forces)
        for(int i=0;i<surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            int off_kx=(k.x-1)*TV::m,off_ky=(k.y-1)*TV::m;
            TV sn=sqrt_coefficients(i)*normal(i);
            for(int j=0;j<TV::m;j++){
                data.Append(TRIPLE<int,int,T>(i,off_kx+j,sn(j)));
                data.Append(TRIPLE<int,int,T>(i,off_ky+j,-sn(j)));}}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR SURFACE_TENSION_FORCE<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    surface.mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR SURFACE_TENSION_FORCE<TV>::
Potential_Energy(const T time) const
{
    T pe=0;
    if(apply_explicit_forces)
        for(int i=0;i<surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            pe+=surface_tension_coefficient*(surface.particles.X(k.x)-surface.particles.X(k.y)).Magnitude();}
    return pe;
}
//#####################################################################
// Function Dump_Curvatures
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
Dump_Curvatures() const
{
    T mn=FLT_MAX,mx=-mn,av=0;
    int n=0;
    ARRAY<TV> F(surface.mesh.elements.m);
    Add_Velocity_Independent_Forces(F,0);
    ARRAY<T> A(coefficients.m);
    for(int i=0;i<surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
        T len=(surface.particles.X(k.x)-surface.particles.X(k.y)).Magnitude()/2;
        A(k.x)+=len;
        A(k.y)+=len;}
    for(int i=0;i<F.m;i++){
        T K=abs(F(i).Magnitude()/A(i)/surface_tension_coefficient);
        TV dx=surface.particles.X(i)-(T).02;
        av+=K;
        n++;
        if(K>mx) mx=K;
        if(K<mn) mn=K;}
    if(n) LOG::cout<<"cstats "<<mn<<"   "<<mx<<std::endl;
    if(coefficients.m) LOG::cout<<"length estimates  "<<coefficients.Max()/coefficients.Min()<<std::endl;
}
//#####################################################################
// Function Test_Diff
//#####################################################################
template<class TV> void SURFACE_TENSION_FORCE<TV>::
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
template class SURFACE_TENSION_FORCE<VECTOR<float,2> >;
template class SURFACE_TENSION_FORCE<VECTOR<double,2> >;
}
