//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_SPRINGS
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Data_Structures/SPARSE_UNION_FIND.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Math_Tools/Robust_Arithmetic.h>
#include <Core/Math_Tools/sign.h>
#include <Core/Matrices/MATRIX_1X1.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <cfloat>
using ::std::sqrt;
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LINEAR_SPRINGS<TV>::
LINEAR_SPRINGS(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh_input)
    :DEFORMABLES_FORCES<TV>(particles),segment_mesh(segment_mesh_input),use_plasticity(false),cache_strain(false),verbose(false)
{
    Set_Stiffness(0);Set_Damping(0);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LINEAR_SPRINGS<TV>::
~LINEAR_SPRINGS()
{
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Stiffness(const T youngs_modulus_input)
{
    constant_youngs_modulus=youngs_modulus_input;
    youngs_modulus.Clean_Memory();
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Stiffness(ARRAY_VIEW<const T> youngs_modulus_input)
{
    constant_youngs_modulus=0;
    youngs_modulus=youngs_modulus_input;
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Restlength
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Restlength(ARRAY_VIEW<const T> restlength_input)
{
    restlength=restlength_input;
    visual_restlength=restlength;
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Damping
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Damping(const T damping_input)
{
    constant_damping=damping_input;
    damping.Clean_Memory();
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Damping
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Damping(ARRAY_VIEW<const T> damping_input)
{
    constant_damping=0;
    damping=damping_input;
    Invalidate_CFL();
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Set_Restlength_From_Particles
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Restlength_From_Particles()
{
    Set_Restlength_From_Material_Coordinates(particles.X);
}
//#####################################################################
// Function Set_Restlength_From_Material_Coordinates
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<TV> material_coordinates)
{
    restlength.Resize(segment_mesh.elements.m,no_init);
    Invalidate_CFL();
    for(int i=0;i<segment_mesh.elements.m;i++)
        restlength(i)=(material_coordinates(segment_mesh.elements(i)(0))-material_coordinates(segment_mesh.elements(i)(1))).Magnitude();
    visual_restlength=restlength;
}
//#####################################################################
// Function Clamp_Restlength
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Clamp_Restlength(const T clamped_restlength)
{
    Invalidate_CFL();
    for(int i=0;i<restlength.m;i++) restlength(i)=max(visual_restlength(i),clamped_restlength);
}
//#####################################################################
// Function Enable_Plasticity
//#####################################################################
template<class TV> template<class T_FIELD> void LINEAR_SPRINGS<TV>::
Enable_Plasticity(const T_FIELD& plastic_yield_strain_input,const T_FIELD& plastic_hardening_input,const T plasticity_clamp_ratio_input)
{
    use_plasticity=true;plasticity_clamp_ratio=plasticity_clamp_ratio_input;
    plastic_yield_strain.Resize(segment_mesh.elements.m,false,false);plastic_yield_strain.Fill(plastic_yield_strain_input);
    plastic_hardening.Resize(segment_mesh.elements.m,false,false);plastic_hardening.Fill(plastic_hardening_input);
    plastic_visual_restlength=visual_restlength;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    segment_mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    Update_Force_Elements(force_segments,segment_mesh.elements,particle_is_simulated);
    if(cache_strain) strains_of_segment.Resize(segment_mesh.elements.m,no_init);
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient)
{
    constant_youngs_modulus=0;
    youngs_modulus.Resize(segment_mesh.elements.m,no_init);
    Invalidate_CFL();
    for(int i=0;i<segment_mesh.elements.m;i++){
        int end1,end2;
        segment_mesh.elements(i).Get(end1,end2);
        T reduced_mass=Pseudo_Inverse(particles.one_over_effective_mass(end1)+particles.one_over_effective_mass(end2));
        youngs_modulus(i)=scaling_coefficient*reduced_mass/restlength(i);}
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Stiffness_Based_On_Reduced_Mass(ARRAY_VIEW<const T> scaling_coefficient)
{
    PHYSBAM_ASSERT(scaling_coefficient.Size()==segment_mesh.elements.m);
    constant_youngs_modulus=0;
    youngs_modulus.Resize(segment_mesh.elements.m,no_init);
    Invalidate_CFL();
    for(int i=0;i<segment_mesh.elements.m;i++){
        int end1,end2;
        segment_mesh.elements(i).Get(end1,end2);
        T reduced_mass=Pseudo_Inverse(particles.one_over_effective_mass(end1)+particles.one_over_effective_mass(end2));
        youngs_modulus(i)=scaling_coefficient(i)*reduced_mass/restlength(i);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;
    damping.Resize(segment_mesh.elements.m,no_init);
    Invalidate_CFL();
    for(int i=0;i<segment_mesh.elements.m;i++){
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(segment_mesh.elements(i)(0))+particles.one_over_effective_mass(segment_mesh.elements(i)(1)));
        T ym;
        if(!youngs_modulus.m) ym=constant_youngs_modulus;
        else ym=youngs_modulus(i);
        damping(i)=overdamping_fraction*2*sqrt(ym*restlength(i)*harmonic_mass);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;
    damping.Resize(segment_mesh.elements.m,no_init);
    Invalidate_CFL();
    for(int i=0;i<segment_mesh.elements.m;i++){
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(segment_mesh.elements(i)(0))+particles.one_over_effective_mass(segment_mesh.elements(i)(1)));
        T ym;
        if(!youngs_modulus.m) ym=constant_youngs_modulus;
        else ym=youngs_modulus(i);
        damping(i)=overdamping_fraction(i)*2*sqrt(ym*restlength(i)*harmonic_mass);}
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;damping.Resize(segment_mesh.elements.m);Invalidate_CFL();
    ARRAY<T> save_damping(damping);Set_Overdamping_Fraction(overdamping_fraction);
    for(int k=0;k<damping.m;k++) damping(k)=max(damping(k),save_damping(k));
}
//#####################################################################
// Function Clamp_Restlength_With_Fraction_Of_Springs
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Clamp_Restlength_With_Fraction_Of_Springs(const T fraction)
{
    ARRAY<T> length(restlength);
    length.Sort();
    T minimum_restlength=length(min((int)(fraction*length.m)+1,length.m));LOG::cout<<"Enlarging the restlength of all linear springs below "<<minimum_restlength<<std::endl;
    Clamp_Restlength(minimum_restlength);
}
//#####################################################################
// Function Enlarge_Restlength_By
//#####################################################################
template <class TV> void LINEAR_SPRINGS<TV>:: 
Enlarge_Restlength_By(const T enlargement_ratio)
{
    restlength=restlength*(1+enlargement_ratio);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    states.Resize(segment_mesh.elements.m,no_init);
    current_lengths.Resize(segment_mesh.elements.m,no_init);
    
    ARRAY_VIEW<const TV> X(particles.X);
    if(!damping.m) for(int s:force_segments){
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        STATE& state=states(s);
        state.nodes=VECTOR<int,2>(node1,node2);
        state.dX=X(node2)-X(node1);
        state.direction=state.dX;
        current_lengths(s)=state.direction.Normalize();
        state.coefficient=constant_damping/restlength(s);}
    else for(int s:force_segments){
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        STATE& state=states(s);
        state.nodes=VECTOR<int,2>(node1,node2);
        state.dX=X(node2)-X(node1);
        state.direction=state.dX;
        current_lengths(s)=state.direction.Normalize();
        state.coefficient=damping(s)/restlength(s);}
    if(use_plasticity) for(int s:force_segments){
        T strain=(current_lengths(s)-visual_restlength(s))/restlength(s);T strain_sign=sign(strain),strain_magnitude=abs(strain);
        if(strain_magnitude>plastic_yield_strain(s)){
            visual_restlength(s)=clamp(current_lengths(s)-strain_sign*plastic_yield_strain(s)*restlength(s),plastic_visual_restlength(s)/plasticity_clamp_ratio,
                plastic_visual_restlength(s)*plasticity_clamp_ratio);
            //visual_restlength(s)=clamp(optimization_current_length(s)-strain_sign*plastic_yield_strain(s)*restlength(s),restlength(s)/plasticity_clamp_ratio,restlength(s)*plasticity_clamp_ratio);
            plastic_yield_strain(s)+=plastic_hardening(s)*(strain_magnitude-plastic_yield_strain(s));}}
    if(compute_half_forces) for(int i=0;i<states.m;i++) states(i).sqrt_coefficient=sqrt(states(i).coefficient);
    if(!segment_mesh.incident_elements) segment_mesh.Initialize_Incident_Elements();
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(!youngs_modulus.m) for(int s:force_segments){
        const STATE& state=states(s);
        TV force=constant_youngs_modulus/restlength(s)*(current_lengths(s)-visual_restlength(s))*state.direction;
        F(state.nodes[0])+=force;F(state.nodes[1])-=force;}
    else for(int s:force_segments){
        const STATE& state=states(s);
        TV force=youngs_modulus(s)/restlength(s)*(current_lengths(s)-visual_restlength(s))*state.direction;
        F(state.nodes[0])+=force;F(state.nodes[1])-=force;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    //int springs_processed=0;
    //LOG::Time("spring add velocity dependent");
    for(int s:force_segments){
        const STATE& state=states(s);
        //springs_processed++;
        TV force=(state.coefficient*TV::Dot_Product(V(state.nodes[0])-V(state.nodes[1]),state.direction))*state.direction;
        F(state.nodes[0])-=force;F(state.nodes[1])+=force;}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose) const
{
    if(!youngs_modulus.m) for(int s:force_segments){
        const STATE& state=states(s);
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV dl=V(node2)-V(node1);
        T l=current_lengths(s),l3=cube(l);
        TV dforce=constant_youngs_modulus*(l-restlength(s))/(restlength(s)*l)*dl+constant_youngs_modulus/l3*TV::Dot_Product(state.dX,dl)*state.dX;
        F(node1)+=dforce;F(node2)-=dforce;}
    else for(int s:force_segments){
        const STATE& state=states(s);
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV dl=V(node2)-V(node1);
        T l=current_lengths(s),l3=cube(l);
        TV dforce=youngs_modulus(s)*(l-restlength(s))/(restlength(s)*l)*dl+youngs_modulus(s)/l3*TV::Dot_Product(state.dX,dl)*state.dX;
        F(node1)+=dforce;F(node2)-=dforce;}
}
//#####################################################################
// Function Velocity_Dependent_Size
//#####################################################################
template<class TV> int LINEAR_SPRINGS<TV>::
Velocity_Dependent_Forces_Size() const
{
    return force_segments.m;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    int aggregate_id=0;
    for(int s:force_segments){
        const STATE& state=states(s);
        aggregate(aggregate_id++)+=state.sqrt_coefficient*TV::Dot_Product(V(state.nodes[0])-V(state.nodes[1]),state.direction);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    int aggregate_id=0;
    for(int s:force_segments){
        const STATE& state=states(s);
        TV force=state.sqrt_coefficient*aggregate(aggregate_id++)*state.direction;
        F(state.nodes[0])+=force;F(state.nodes[1])-=force;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
    for(int s:force_segments){
        const VECTOR<int,2>& nodes=segment_mesh.elements(s);
        T ym=youngs_modulus.m?youngs_modulus(s):constant_youngs_modulus;
        T d=damping.m?damping(s):constant_damping;
        for(int k=0;k<2;k++){
            frequency(nodes[k]).elastic_squared+=particles.one_over_effective_mass(nodes[k])/restlength(s)*4*ym*one_over_cfl_number_squared;
            frequency(nodes[k]).damping+=particles.one_over_effective_mass(nodes[k])/restlength(s)*2*d*one_over_cfl_number;}}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0,strain_rate;TV dx;
    if(use_rest_state_for_strain_rate) for(int s:force_segments){
        int i,j;segment_mesh.elements(s).Get(i,j);
        dx=particles.X(j)-particles.X(i);T magnitude=dx.Magnitude();if(magnitude!=0) dx.Normalize();
        strain_rate=TV::Dot_Product(particles.V(j)-particles.V(i),dx)/restlength(s);
        max_strain_rate=max(max_strain_rate,abs(strain_rate));
        if(cache_strain){strains_of_segment(s)=VECTOR<T,2>(abs(strain_rate),abs((magnitude-visual_restlength(s))/restlength(s)));}}
    else for(int s:force_segments){
        int i,j;segment_mesh.elements(s).Get(i,j);
        dx=particles.X(j)-particles.X(i);
        strain_rate=TV::Dot_Product(particles.V(j)-particles.V(i),dx)/TV::Dot_Product(dx,dx);
        max_strain_rate=max(max_strain_rate,abs(strain_rate));
        if(cache_strain){strains_of_segment(s)=VECTOR<T,2>(abs(strain_rate),abs((dx.Magnitude()-visual_restlength(s))/restlength(s)));}}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Average_Restlength
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Average_Restlength() const
{
    return restlength.Average();
}
//#####################################################################
// Function Print_Restlength_Statistics
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Print_Restlength_Statistics() const
{
    LOG::SCOPE scope("linear spring statistics","linear spring statistics");
    LOG::Stat("count",restlength.m);
    ARRAY<T> length(restlength);
    length.Sort();
    ARRAY<T> visual_length(visual_restlength);
    visual_length.Sort();
    if(length.m){
        LOG::Stat("smallest restlength",length(0));LOG::Stat("smallest visual restlength",visual_length(0));
        LOG::Stat("one percent restlength",length((int)(.01*length.m)+1));LOG::Stat("one percent visual restlength",visual_length((int)(.01*length.m)+1));
        LOG::Stat("ten percent restlength",length((int)(.1*length.m)+1));LOG::Stat("ten percent visual restlength",visual_length((int)(.1*length.m)+1));
        LOG::Stat("median restlength",length((int)(.5*length.m)+1));LOG::Stat("median visual restlength",visual_length((int)(.5*length.m)+1));}
}
//#####################################################################
// Function Print_Deformation_Statistics
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Print_Deformation_Statistics() const
{
    LOG::SCOPE scope("linear spring deformation","linear spring deformation");
    ARRAY<T> deformation(segment_mesh.elements.m,no_init);
    for(int s=0;s<segment_mesh.elements.m;s++){
        int i,j;segment_mesh.elements(s).Get(i,j);
        T length=(particles.X(i)-particles.X(j)).Magnitude(),rl=visual_restlength(s);
        deformation(s)=rl?abs(length-rl)/rl:length==0?0:FLT_MAX;}
    deformation.Sort();
    LOG::Stat("maximum deformation",deformation.Last());
    LOG::Stat("one percent deformation",deformation((int)(.99*deformation.m)+1));
    LOG::Stat("ten percent deformation",deformation((int)(.9*deformation.m)+1));
    LOG::Stat("twenty percent deformation",deformation((int)(.8*deformation.m)+1));
    LOG::Stat("median deformation",deformation((int)(.5*deformation.m)+1));
}
//#####################################################################
// Function Maximum_Compression_Or_Expansion_Fraction
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Maximum_Compression_Or_Expansion_Fraction(int* index) const
{
    T max_compression=0;int max_index=-1;
    for(int s=0;s<segment_mesh.elements.m;s++){
        int i,j;segment_mesh.elements(s).Get(i,j);
        T length=(particles.X(i)-particles.X(j)).Magnitude();
        T rl=visual_restlength(s);
        T compression=(rl)?(abs(length-rl)/rl):((length==0)?0:FLT_MAX);
        if(compression>max_compression){max_compression=compression;max_index=s;}}
    if(index) (*index)=max_index;
    return max_compression;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Potential_Energy(const int s,const T time) const
{
    T current_length=(particles.X(segment_mesh.elements(s)(0))-particles.X(segment_mesh.elements(s)(1))).Magnitude();
    T spring_youngs_modulus=youngs_modulus.m?youngs_modulus(s):constant_youngs_modulus;
    return (T).5*spring_youngs_modulus/restlength(s)*sqr(current_length-visual_restlength(s));
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(int s:force_segments){
        potential_energy+=Potential_Energy(s,time);}
    return potential_energy;
}
//#####################################################################
// Function Add_Force_Data
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{   
    ARRAY_VIEW<const TV> X=particles.X;
    FORCE_DATA<TV> force_data;
    if(force_name.empty()) force_data.name="LINEAR_SPRINGS";
    else force_data.name=force_name;
    for(int s:force_segments){
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        force_data.first_action_point=X(node1);
        force_data.second_action_point=X(node2);
        T current_length=(X(node2)-X(node1)).Magnitude();
        force_data.state=visual_restlength(s)?(current_length-visual_restlength(s))/visual_restlength(s):(T)0;
        force_data_list.Append(force_data);}
}
//#####################################################################
// Function Endpoint_Velocity
//#####################################################################
template<class TV> TV LINEAR_SPRINGS<TV>::
Endpoint_Velocity(int s,int b) const
{
    return particles.V(segment_mesh.elements(s)(b));
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Endpoint_Kinetic_Energy(int s,int b) const
{
    return (T).5*particles.mass(segment_mesh.elements(s)(b))*particles.V(segment_mesh.elements(s)(b)).Magnitude_Squared();
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Endpoint_Kinetic_Energy(int s) const
{
    return Endpoint_Kinetic_Energy(s,1)+Endpoint_Kinetic_Energy(s,2);
}
//#####################################################################
// Function Effective_Impulse_Factor
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Effective_Impulse_Factor(int s) const
{
    T one_over_denom=particles.one_over_mass(segment_mesh.elements(s)(0))*particles.one_over_mass(segment_mesh.elements(s)(1));
    return (particles.mass(segment_mesh.elements(s)(0))+particles.mass(segment_mesh.elements(s)(0)))*one_over_denom;
}
//#####################################################################
// Function Create_Edge_Springs
//#####################################################################
template<class TV> LINEAR_SPRINGS<TV>*
Create_Edge_Springs(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,
    const typename TV::SCALAR stiffness,const typename TV::SCALAR overdamping_fraction,
    const bool limit_time_step_by_strain_rate,const typename TV::SCALAR max_strain_per_time_step,
    const bool use_rest_state_for_strain_rate,const typename TV::SCALAR restlength_enlargement_fraction,
    const bool verbose,const typename TV::SCALAR rest_length_material_length_enlargement_ratio)
{
    LINEAR_SPRINGS<TV>* ls=new LINEAR_SPRINGS<TV>(particles,segment_mesh);
    ls->Set_Restlength_From_Particles();
    if(restlength_enlargement_fraction) ls->Clamp_Restlength_With_Fraction_Of_Springs(restlength_enlargement_fraction);
    if(rest_length_material_length_enlargement_ratio!=0) ls->Enlarge_Restlength_By(rest_length_material_length_enlargement_ratio);
    ls->Set_Stiffness(stiffness);
    ls->Set_Overdamping_Fraction(overdamping_fraction);
    ls->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    ls->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    if(verbose) ls->Print_Restlength_Statistics();
    ls->verbose=verbose;
    return ls;
}
//#####################################################################
// Function Create_Edge_Springs
//#####################################################################
template<class T_OBJECT> LINEAR_SPRINGS<typename T_OBJECT::VECTOR_T>*
Create_Edge_Springs(T_OBJECT& object,
    const typename T_OBJECT::SCALAR stiffness,const typename T_OBJECT::SCALAR overdamping_fraction,
    const bool limit_time_step_by_strain_rate,const typename T_OBJECT::SCALAR max_strain_per_time_step,
    const bool use_rest_state_for_strain_rate,const typename T_OBJECT::SCALAR restlength_enlargement_fraction,
    const bool verbose,const typename T_OBJECT::SCALAR rest_length_material_length_enlargement_ratio)
{
    return Create_Edge_Springs(dynamic_cast<DEFORMABLE_PARTICLES<typename T_OBJECT::VECTOR_T>&>(object.particles),object.Get_Segment_Mesh(),stiffness,overdamping_fraction,limit_time_step_by_strain_rate,max_strain_per_time_step,
        use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose,rest_length_material_length_enlargement_ratio);
}
template class LINEAR_SPRINGS<VECTOR<double,2> >;
template class LINEAR_SPRINGS<VECTOR<double,3> >;
template class LINEAR_SPRINGS<VECTOR<float,2> >;
template class LINEAR_SPRINGS<VECTOR<float,3> >;
template LINEAR_SPRINGS<TETRAHEDRALIZED_VOLUME<double>::VECTOR_T>*
Create_Edge_Springs<TETRAHEDRALIZED_VOLUME<double> >(TETRAHEDRALIZED_VOLUME<double>&,
    TETRAHEDRALIZED_VOLUME<double>::SCALAR,TETRAHEDRALIZED_VOLUME<double>::SCALAR,bool,
    TETRAHEDRALIZED_VOLUME<double>::SCALAR,bool,TETRAHEDRALIZED_VOLUME<double>::SCALAR,
    bool,TETRAHEDRALIZED_VOLUME<double>::SCALAR);
template LINEAR_SPRINGS<TETRAHEDRALIZED_VOLUME<float>::VECTOR_T>*
Create_Edge_Springs<TETRAHEDRALIZED_VOLUME<float> >(TETRAHEDRALIZED_VOLUME<float>&,
    TETRAHEDRALIZED_VOLUME<float>::SCALAR,TETRAHEDRALIZED_VOLUME<float>::SCALAR,bool,
    TETRAHEDRALIZED_VOLUME<float>::SCALAR,bool,TETRAHEDRALIZED_VOLUME<float>::SCALAR,
    bool,TETRAHEDRALIZED_VOLUME<float>::SCALAR);
template LINEAR_SPRINGS<TRIANGULATED_AREA<double>::VECTOR_T>*
Create_Edge_Springs<TRIANGULATED_AREA<double> >(TRIANGULATED_AREA<double>&,
    TRIANGULATED_AREA<double>::SCALAR,TRIANGULATED_AREA<double>::SCALAR,bool,
    TRIANGULATED_AREA<double>::SCALAR,bool,TRIANGULATED_AREA<double>::SCALAR,
    bool,TRIANGULATED_AREA<double>::SCALAR);
template LINEAR_SPRINGS<TRIANGULATED_SURFACE<float>::VECTOR_T>*
Create_Edge_Springs<TRIANGULATED_SURFACE<float> >(TRIANGULATED_SURFACE<float>&,
    TRIANGULATED_SURFACE<float>::SCALAR,TRIANGULATED_SURFACE<float>::SCALAR,bool,
    TRIANGULATED_SURFACE<float>::SCALAR,bool,TRIANGULATED_SURFACE<float>::SCALAR,
    bool,TRIANGULATED_SURFACE<float>::SCALAR);
template LINEAR_SPRINGS<SEGMENTED_CURVE<VECTOR<double,3> >::VECTOR_T>*
Create_Edge_Springs<SEGMENTED_CURVE<VECTOR<double,3> > >(SEGMENTED_CURVE<VECTOR<double,3> >&,
    SEGMENTED_CURVE<VECTOR<double,3> >::SCALAR,SEGMENTED_CURVE<VECTOR<double,3> >::SCALAR,
    bool,SEGMENTED_CURVE<VECTOR<double,3> >::SCALAR,bool,SEGMENTED_CURVE<VECTOR<double,3> >::SCALAR,
    bool,SEGMENTED_CURVE<VECTOR<double,3> >::SCALAR);
template LINEAR_SPRINGS<SEGMENTED_CURVE<VECTOR<float,2> >::VECTOR_T>*
Create_Edge_Springs<SEGMENTED_CURVE<VECTOR<float,2> > >(SEGMENTED_CURVE<VECTOR<float,2> >&,
    SEGMENTED_CURVE<VECTOR<float,2> >::SCALAR,SEGMENTED_CURVE<VECTOR<float,2> >::SCALAR,
    bool,SEGMENTED_CURVE<VECTOR<float,2> >::SCALAR,bool,SEGMENTED_CURVE<VECTOR<float,2> >::SCALAR,
    bool,SEGMENTED_CURVE<VECTOR<float,2> >::SCALAR);
template LINEAR_SPRINGS<SEGMENTED_CURVE<VECTOR<float,3> >::VECTOR_T>*
Create_Edge_Springs<SEGMENTED_CURVE<VECTOR<float,3> > >(SEGMENTED_CURVE<VECTOR<float,3> >&,
    SEGMENTED_CURVE<VECTOR<float,3> >::SCALAR,SEGMENTED_CURVE<VECTOR<float,3> >::SCALAR,
    bool,SEGMENTED_CURVE<VECTOR<float,3> >::SCALAR,bool,SEGMENTED_CURVE<VECTOR<float,3> >::SCALAR,
    bool,SEGMENTED_CURVE<VECTOR<float,3> >::SCALAR);
template LINEAR_SPRINGS<SEGMENTED_CURVE_2D<double>::VECTOR_T>*
Create_Edge_Springs<SEGMENTED_CURVE_2D<double> >(SEGMENTED_CURVE_2D<double>&,
    SEGMENTED_CURVE_2D<double>::SCALAR,SEGMENTED_CURVE_2D<double>::SCALAR,bool,
    SEGMENTED_CURVE_2D<double>::SCALAR,bool,SEGMENTED_CURVE_2D<double>::SCALAR,
    bool,SEGMENTED_CURVE_2D<double>::SCALAR);
template LINEAR_SPRINGS<TRIANGULATED_SURFACE<double>::VECTOR_T>*
Create_Edge_Springs<TRIANGULATED_SURFACE<double> >(TRIANGULATED_SURFACE<double>&,
    TRIANGULATED_SURFACE<double>::SCALAR,TRIANGULATED_SURFACE<double>::SCALAR,bool,
    TRIANGULATED_SURFACE<double>::SCALAR,bool,TRIANGULATED_SURFACE<double>::SCALAR,
    bool,TRIANGULATED_SURFACE<double>::SCALAR);
template LINEAR_SPRINGS<VECTOR<float,3> >* Create_Edge_Springs<VECTOR<float,3> >(
    DEFORMABLE_PARTICLES<VECTOR<float,3> >&,SEGMENT_MESH&,VECTOR<float,3>::SCALAR,
    VECTOR<float,3>::SCALAR,bool,VECTOR<float,3>::SCALAR,bool,VECTOR<float,3>::SCALAR,
    bool,VECTOR<float,3>::SCALAR);
template LINEAR_SPRINGS<VECTOR<double,3> >* Create_Edge_Springs<VECTOR<double,3> >(
    DEFORMABLE_PARTICLES<VECTOR<double,3> >&,SEGMENT_MESH&,VECTOR<double,3>::SCALAR,
    VECTOR<double,3>::SCALAR,bool,VECTOR<double,3>::SCALAR,bool,VECTOR<double,3>::SCALAR,
    bool,VECTOR<double,3>::SCALAR);
//#####################################################################
}
