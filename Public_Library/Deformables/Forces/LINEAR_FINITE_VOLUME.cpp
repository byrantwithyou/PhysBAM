//#####################################################################
// Copyright 2004-2007, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/SPARSE_UNION_FIND.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Core/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Forces/LINEAR_FINITE_VOLUME.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> LINEAR_FINITE_VOLUME<TV,d>::
LINEAR_FINITE_VOLUME(T_OBJECT& object,const T youngs_modulus,const T poissons_ratio,const T Rayleigh_coefficient)
    :DEFORMABLES_FORCES<TV>(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(object.particles)),object(object),mesh(object.mesh),use_uniform_density(false),density_list(0)
{
    assert(poissons_ratio!=-1 && poissons_ratio!=.5);
    lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    mu=youngs_modulus/(2*(1+poissons_ratio));
    alpha=Rayleigh_coefficient*lambda;beta=Rayleigh_coefficient*mu;
    Initialize_Material_State(particles.X);
    mesh.Initialize_Incident_Elements();
    if(use_uniform_density){
        T total_mass=0;
        ARRAY<int> mesh_particles;
        Get_Unique(mesh_particles,mesh.elements.Flattened());
        for(int i=0;i<mesh_particles.m;i++) total_mass+=particles.mass(mesh_particles(i));
        density=total_mass/object.Total_Size();
        if(density==0) density=TV::m==1?1:TV::m==2?100:1000;}
    else{
        density_list=new ARRAY<T>(mesh.elements.m);
        for(int i=0;i<mesh.elements.m;i++){
            const VECTOR<int,d+1>& nodes=mesh.elements(i);
            T volume=object.Signed_Size(i);
            for(int j=0;j<nodes.m;j++) (*density_list)(i)+=particles.mass(nodes(j))/(*mesh.incident_elements)(nodes(j)).m/volume;
            if((*density_list)(i)==0) (*density_list)(i)=TV::m==1?1:TV::m==2?100:1000;}}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> LINEAR_FINITE_VOLUME<TV,d>::
~LINEAR_FINITE_VOLUME()
{}
//#####################################################################
// Function Initialize_Material_State
//#####################################################################
namespace{
template<class T,int d> inline VECTOR<T,d> Normal(const MATRIX<T,d>& A)
{
    return VECTOR<T,d>();
}
template<class T> inline VECTOR<T,3> Normal(const MATRIX<T,3,2>& A)
{
    return A.Weighted_Normal().Normalized();
}
template<class T,int d> inline MATRIX<T,d> Pseudoinverse(const MATRIX<T,d>& A)
{
    return A.Inverse();
}
template<class T> inline MATRIX<T,2,3> Pseudoinverse(const MATRIX<T,3,2>& A)
{
    T scale=A.Parallelepiped_Measure();assert(scale);
    return (T)1/scale*A.Cofactor_Matrix().Transposed();
}
}
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Initialize_Material_State(ARRAY_VIEW<const TV> X)
{
    Dm_inverse.Resize(mesh.elements.m,no_init);
    Bm.Resize(mesh.elements.m,no_init);
    if(TV::m>d) normals.Resize(mesh.elements.m,no_init);
    for(int t=0;t<mesh.elements.m;t++){
        MATRIX<T,TV::m,d> Dm=Ds(X,t);
        if(TV::m>d) normals(t)=Normal(Dm);
        Dm_inverse(t)=Pseudoinverse(Dm);
        Bm(t)=-(T)1/factorial(d)*Dm.Cofactor_Matrix();}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    Update_Force_Elements(force_elements,mesh.elements,particle_is_simulated);
    Update_Force_Particles(force_particles,mesh.elements.Flattened(),particle_is_simulated,true);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int t:force_elements){
        MATRIX<T,TV::m,d> G=Stress(t)*Bm(t);
        STRAIN_MEASURE<TV,d>::Distribute_Force(F,mesh.elements(t),G);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(int t:force_elements){
        SYMMETRIC_MATRIX<T,TV::m> cauchy_strain_rate=(Ds(V,t)*Dm_inverse(t)).Symmetric_Part();
        MATRIX<T,TV::m,d> G=(2*beta*cauchy_strain_rate+alpha*cauchy_strain_rate.Trace())*Bm(t);
        STRAIN_MEASURE<TV,d>::Distribute_Force(F,mesh.elements(t),G);}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose) const
{
    for(int t:force_elements){
        MATRIX<T,TV::m,d> dG=Stress_Differential(V,t)*Bm(t);
        STRAIN_MEASURE<TV,d>::Distribute_Force(F,mesh.elements(t),dG);}
}
//#####################################################################
// Function Intialize_CFL
//#####################################################################
namespace{
template<class T,int d> inline T Simplex_Minimum_Altitude(const MATRIX<T,d>& Dm_inverse)
{
    return Dm_inverse.Inverse().Simplex_Minimum_Altitude();
}
template<class T> inline T Simplex_Minimum_Altitude(const MATRIX<T,2,3>& Dm_inverse)
{
    return Dm_inverse.Transposed().R_From_QR_Factorization().Inverse().Simplex_Minimum_Altitude();
}
}
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    // TODO: MPI
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);

    ARRAY<FREQUENCY_DATA> fragment_particle_frequency(frequency.Size());
    for(int t:force_elements){
        T local_density=use_uniform_density?density:(*density_list)(t);
        T altitude_squared=sqr(Simplex_Minimum_Altitude(Dm_inverse(t)));
        T elastic_squared=(lambda+2*mu)/(local_density*altitude_squared)*one_over_cfl_number_squared;
        T damping=(alpha+2*beta)/(local_density*altitude_squared)*one_over_cfl_number;
        const VECTOR<int,d+1>& nodes=mesh.elements(t);
        for(int j=0;j<nodes.m;j++){FREQUENCY_DATA& data=fragment_particle_frequency(nodes[j]);
            data.elastic_squared=max(data.elastic_squared,elastic_squared);data.damping=max(data.damping,damping);}}
    for(int p:force_particles){
        frequency(p).elastic_squared+=fragment_particle_frequency(p).elastic_squared;
        frequency(p).damping+=fragment_particle_frequency(p).damping;}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV,int d> typename TV::SCALAR LINEAR_FINITE_VOLUME<TV,d>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0;
    for(int t:force_elements){
        max_strain_rate=max(max_strain_rate,(Ds(particles.V,t)*Dm_inverse(t)).Max_Abs());}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
namespace PhysBAM{
template class LINEAR_FINITE_VOLUME<VECTOR<float,2>,2>;
template class LINEAR_FINITE_VOLUME<VECTOR<float,3>,2>;
template class LINEAR_FINITE_VOLUME<VECTOR<float,3>,3>;
template class LINEAR_FINITE_VOLUME<VECTOR<double,2>,2>;
template class LINEAR_FINITE_VOLUME<VECTOR<double,3>,2>;
template class LINEAR_FINITE_VOLUME<VECTOR<double,3>,3>;
}
