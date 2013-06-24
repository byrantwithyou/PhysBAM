//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <Tools/Math_Tools/cyclic_shift.h>
#include <Tools/Math_Tools/Robust_Arithmetic.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Forces/AXIAL_BENDING_SPRINGS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> AXIAL_BENDING_SPRINGS<T>::
AXIAL_BENDING_SPRINGS(DEFORMABLE_PARTICLES<TV>& particles_input,TRIANGLE_MESH& triangle_mesh_input)
    :DEFORMABLES_FORCES<TV>(particles_input),triangle_mesh(triangle_mesh_input),verbose(false)
{
    Initialize();
    Set_Stiffness(0);
    restlength.Fill((T)0);
    visual_restlength.Fill((T)0);
    Set_Damping(0);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> AXIAL_BENDING_SPRINGS<T>::
~AXIAL_BENDING_SPRINGS()
{}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    triangle_mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_springs.Update(spring_particles,particle_is_simulated);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Initialize()
{
    bool adjacent_triangles_defined=(triangle_mesh.adjacent_elements!=0);if(!adjacent_triangles_defined) triangle_mesh.Initialize_Adjacent_Elements();

    int number_quadruples=0;
    for(int t=0;t<triangle_mesh.elements.m;t++) for(int a=0;a<(*triangle_mesh.adjacent_elements)(t).m;a++) if((*triangle_mesh.adjacent_elements)(t)(a)>t) number_quadruples++;
    spring_particles.Resize(number_quadruples,false);youngs_modulus.Resize(number_quadruples);restlength.Resize(number_quadruples);
    visual_restlength.Resize(number_quadruples);damping.Resize(number_quadruples);attached_edge_length.Resize(number_quadruples);attached_edge_restlength.Resize(number_quadruples);
    int index=0; // reset number
    for(int t=0;t<triangle_mesh.elements.m;t++){
        int t1,t2,t3;triangle_mesh.elements(t).Get(t1,t2,t3);
        for(int a=0;a<(*triangle_mesh.adjacent_elements)(t).m;a++){
            int s=(*triangle_mesh.adjacent_elements)(t)(a);
            if(s>t){
                int s1,s2,s3;triangle_mesh.elements(s).Get(s1,s2,s3);
                if(t1==s1 || t1==s2 || t1==s3){cyclic_shift(t1,t2,t3);if(t1==s1 || t1==s2 || t1==s3) cyclic_shift(t1,t2,t3);}
                spring_particles(index++).Set(t2,t3,t1,triangle_mesh.Other_Node(t2,t3,s));}}}

    if(!adjacent_triangles_defined){delete triangle_mesh.adjacent_elements;triangle_mesh.adjacent_elements=0;}
}
//#####################################################################
// Function Set_Restlength_From_Particles
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Set_Restlength_From_Particles()
{
    for(int s=0;s<spring_particles.m;s++){const VECTOR<int,4>& nodes=spring_particles(s);
        TV axial_direction;VECTOR<T,2> weights;
        Axial_Vector(nodes,visual_restlength(s),axial_direction,weights,attached_edge_restlength(s));
        restlength(s)=visual_restlength(s);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    Invalidate_CFL();
    for(int s=0;s<spring_particles.m;s++){
        T harmonic_mass=Pseudo_Divide((T)4,particles.one_over_effective_mass.Subset(spring_particles(s)).Sum());
        damping(s)=overdamping_fraction*2*sqrt(youngs_modulus(s)*restlength(s)*harmonic_mass);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction) // 1 is critically damped
{
    Invalidate_CFL();
    for(int s=0;s<spring_particles.m;s++){
        T harmonic_mass=Pseudo_Divide((T)4,particles.one_over_effective_mass.Subset(spring_particles(s)).Sum());
        damping(s)=overdamping_fraction(s)*2*sqrt(youngs_modulus(s)*restlength(s)*harmonic_mass);}
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    Invalidate_CFL();
    ARRAY<T> save_damping(damping);Set_Overdamping_Fraction(overdamping_fraction);
    for(int k=0;k<damping.m;k++) damping(k)=max(damping(k),save_damping(k));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    optimization_current_length.Resize(spring_particles.m,false,false);optimization_direction.Resize(spring_particles.m,false,false);
    optimization_weights.Resize(spring_particles.m,false,false);optimization_coefficient.Resize(spring_particles.m,false,false);

    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,4>& nodes=spring_particles(s);
        VECTOR<T,2> weights;
        Axial_Vector(nodes,optimization_current_length(s),optimization_direction(s),weights,attached_edge_length(s));
        optimization_weights(s).Set(1-weights.x,weights.x,1-weights.y,weights.y);
        optimization_coefficient(s)=damping(s)/restlength(s);}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2,node3,node4;spring_particles(s).Get(node1,node2,node3,node4);
        T w1,w2,w3,w4;optimization_weights(s).Get(w1,w2,w3,w4);
        TV force=youngs_modulus(s)/restlength(s)/*sqr(attached_edge_length(s)/attached_edge_restlength(s)-1)*/*(optimization_current_length(s)-visual_restlength(s))*optimization_direction(s);
        F(node1)-=w1*force;F(node2)-=w2*force;F(node3)+=w3*force;F(node4)+=w4*force;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2,node3,node4;spring_particles(s).Get(node1,node2,node3,node4);
        T w1,w2,w3,w4;optimization_weights(s).Get(w1,w2,w3,w4);
        TV force=(optimization_coefficient(s)*TV::Dot_Product(w1*V(node1)+w2*V(node2)-w3*V(node3)-w4*V(node4),optimization_direction(s)))*optimization_direction(s);
        F(node1)-=w1*force;F(node2)-=w2*force;F(node3)+=w3*force;F(node4)+=w4*force;}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2,node3,node4;spring_particles(s).Get(node1,node2,node3,node4);
        T w1,w2,w3,w4;optimization_weights(s).Get(w1,w2,w3,w4);
        TV dl=w1*V(node1)+w2*V(node2)-w3*V(node3)-w4*V(node4);
        TV force=youngs_modulus(s)/restlength(s)/*sqr(attached_edge_length(s)/attached_edge_restlength(s)-1)*/*dl.Projected_On_Unit_Direction(optimization_direction(s));
        F(node1)-=w1*force;F(node2)-=w2*force;F(node3)+=w3*force;F(node4)+=w4*force;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,4>& nodes=spring_particles(s);
        T one_over_mass_times_restlength=(T).25/restlength(s)*particles.one_over_effective_mass.Subset(nodes).Sum();
        T elastic_hertz_squared=4*youngs_modulus(s)*one_over_mass_times_restlength*one_over_cfl_number_squared;
        T damping_hertz=2*damping(s)*one_over_mass_times_restlength*one_over_cfl_number;
        for(int k=0;k<4;k++){frequency(nodes[k]).elastic_squared+=elastic_hertz_squared;frequency(nodes[k]).damping+=damping_hertz;}}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
CFL_Strain_Rate() const
{
    ARRAY_VIEW<const TV> V(particles.V);
    T max_strain_rate=0;
    if(use_rest_state_for_strain_rate) for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,4>& nodes=spring_particles(s);
        TV dx;VECTOR<T,2> weights;T current_length;T current_edge_length;
        Axial_Vector(nodes,current_length,dx,weights,current_edge_length); // dx is normalized
        T strain_rate=TV::Dot_Product((1-weights.x)*V(nodes[0])+weights.x*V(nodes[1])-(1-weights.y)*V(nodes[2])-weights.y*V(nodes[3]),dx)/restlength(s);
        max_strain_rate=max(max_strain_rate,abs(strain_rate));}
    else for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,4>& nodes=spring_particles(s);
        TV dx;VECTOR<T,2> weights;T current_length;T current_edge_length;
        Axial_Vector(nodes,current_length,dx,weights,current_edge_length); // dx is normalized
        T strain_rate=TV::Dot_Product((1-weights.x)*V(nodes[0])+weights.x*V(nodes[1])-(1-weights.y)*V(nodes[2])-weights.y*V(nodes[3]),dx)/current_length;
        max_strain_rate=max(max_strain_rate,abs(strain_rate));}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
// assumes mass and restlength are already defined
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient)
{
    for(int s=0;s<spring_particles.m;s++){
        T harmonic_mass=Pseudo_Divide((T)4,particles.one_over_effective_mass.Subset(spring_particles(s)).Sum());
        youngs_modulus(s)=scaling_coefficient*harmonic_mass/restlength(s);}
}
//#####################################################################
// Function Axial_Vector
//#####################################################################
// direction points from cross edge to common edge
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Axial_Vector(const VECTOR<int,4>& nodes,T& axial_length,TV& axial_direction,VECTOR<T,2>& weights,T& attached_edge_length) const
{
    axial_direction=SEGMENT_3D<T>(particles.X(nodes[0]),particles.X(nodes[1])).Shortest_Vector_Between_Segments(SEGMENT_3D<T>(particles.X(nodes[2]),particles.X(nodes[3])),weights);
    axial_length=axial_direction.Normalize();
    // TODO: How do we choose the sign of axial_direction?
    if(!axial_length) axial_direction=TV::Cross_Product(particles.X(nodes[1])-particles.X(nodes[0]),particles.X(nodes[3])-particles.X(nodes[2])).Normalized();
    attached_edge_length=(particles.X(nodes[0])-particles.X(nodes[2])).Magnitude()+
        (particles.X(nodes[0])-particles.X(nodes[3])).Magnitude()+
        (particles.X(nodes[1])-particles.X(nodes[2])).Magnitude()+
        (particles.X(nodes[1])-particles.X(nodes[3])).Magnitude();
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Potential_Energy(int s,const T time) const
{
    const VECTOR<int,4>& nodes=spring_particles(s);
    VECTOR<T,2> weights;
    TV axial_direction=SEGMENT_3D<T>(particles.X(nodes[0]),particles.X(nodes[1])).Shortest_Vector_Between_Segments(SEGMENT_3D<T>(particles.X(nodes[2]),particles.X(nodes[3])),weights);
    T axial_length=axial_direction.Normalize();
    return (T).5*youngs_modulus(s)/restlength(s)*sqr(axial_length-visual_restlength(s));
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        potential_energy+=Potential_Energy(s,time);}
    return potential_energy;
}
//#####################################################################
// Function Endpoint_Velocity
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Endpoint_Mass(int s,int b) const
{
    switch(b){
        case 0:
            return particles.mass(spring_particles(s)(0))*optimization_weights(s)(0)+particles.mass(spring_particles(s)(1))*optimization_weights(s)(1);
        case 1:
            return particles.mass(spring_particles(s)(2))*optimization_weights(s)(2)+particles.mass(spring_particles(s)(3))*optimization_weights(s)(3);
        default:
            PHYSBAM_FATAL_ERROR("Invalid endpoint");}
}
//#####################################################################
// Function Endpoint_Velocity
//#####################################################################
template<class T> VECTOR<T,3> AXIAL_BENDING_SPRINGS<T>::
Endpoint_Velocity(int s,int b) const
{
    return Endpoint_Velocity(particles.V,s,b);
}
//#####################################################################
// Function Endpoint_Velocity
//#####################################################################
template<class T> VECTOR<T,3> AXIAL_BENDING_SPRINGS<T>::
Endpoint_Velocity(ARRAY_VIEW<const TV> velocity,int s,int b) const
{
    switch(b){
        case 0:{
            int node1=spring_particles(s)(0),node2=spring_particles(s)(1);
            return velocity(node1)*optimization_weights(s)(0)+velocity(node2)*optimization_weights(s)(1);}
        case 1:{
            int node3=spring_particles(s)(2),node4=spring_particles(s)(3);
            return velocity(node3)*optimization_weights(s)(2)+velocity(node4)*optimization_weights(s)(3);}
        default:
            PHYSBAM_FATAL_ERROR("Invalid endpoint");}
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Endpoint_Kinetic_Energy(int s,int b) const
{
    return Endpoint_Kinetic_Energy(particles.V,s,b);
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Endpoint_Kinetic_Energy(ARRAY_VIEW<const TV> velocity,int s,int b) const
{
    return (T).5*Endpoint_Mass(s,b)*Endpoint_Velocity(velocity,s,b).Magnitude_Squared();
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Endpoint_Kinetic_Energy(int s) const
{
    return Endpoint_Kinetic_Energy(s,1)+Endpoint_Kinetic_Energy(s,2);
}
//#####################################################################
// Function Effective_Impulse_Factor
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Effective_Impulse_Factor(int s) const
{
    T mass1=Endpoint_Mass(s,1);
    T mass2=Endpoint_Mass(s,2);
    T one_over_denom=1/(mass1*mass2);
    return (mass1+mass2)*one_over_denom;
}
//#####################################################################
// Function Create_Axial_Bending_Springs
//#####################################################################
template<class T> AXIAL_BENDING_SPRINGS<T>* PhysBAM::
Create_Axial_Bending_Springs(DEFORMABLE_PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& triangle_mesh,const T clamped_restlength,const T stiffness,const T overdamping_fraction,
    const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate,const bool verbose)
{
    AXIAL_BENDING_SPRINGS<T>* axial=new AXIAL_BENDING_SPRINGS<T>(particles,triangle_mesh);
    axial->verbose=verbose;
    axial->Set_Restlength_From_Particles();
    axial->Clamp_Restlength(clamped_restlength);
    axial->Set_Stiffness(stiffness);
    axial->Set_Overdamping_Fraction(overdamping_fraction);
    axial->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    axial->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    return axial;
}
//#####################################################################
// Function Create_Axial_Bending_Springs
//#####################################################################
template<class T> AXIAL_BENDING_SPRINGS<T>* PhysBAM::
Create_Axial_Bending_Springs(TRIANGULATED_SURFACE<T>& triangulated_surface,
    const T clamped_restlength,const T stiffness,const T overdamping_fraction,const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,
    const bool use_rest_state_for_strain_rate,const bool verbose)
{
    return Create_Axial_Bending_Springs(dynamic_cast<DEFORMABLE_PARTICLES<VECTOR<T,3> >&>(triangulated_surface.particles),triangulated_surface.mesh,clamped_restlength,stiffness,overdamping_fraction,limit_time_step_by_strain_rate,
        max_strain_per_time_step,use_rest_state_for_strain_rate,verbose);
}
//#####################################################################
namespace PhysBAM{
template class AXIAL_BENDING_SPRINGS<float>;
template AXIAL_BENDING_SPRINGS<float>* Create_Axial_Bending_Springs<float>(TRIANGULATED_SURFACE<float>&,float,float,float,bool,float,bool,bool);
template AXIAL_BENDING_SPRINGS<float>* Create_Axial_Bending_Springs<float>(DEFORMABLE_PARTICLES<VECTOR<float,3> >&,TRIANGLE_MESH&,float,float,float,bool,float,bool,bool);
template class AXIAL_BENDING_SPRINGS<double>;
template AXIAL_BENDING_SPRINGS<double>* Create_Axial_Bending_Springs<double>(TRIANGULATED_SURFACE<double>&,double,double,double,bool,double,bool,bool);
template AXIAL_BENDING_SPRINGS<double>* Create_Axial_Bending_Springs<double>(DEFORMABLE_PARTICLES<VECTOR<double,3> >&,TRIANGLE_MESH&,double,double,double,bool,double,bool,bool);
}
