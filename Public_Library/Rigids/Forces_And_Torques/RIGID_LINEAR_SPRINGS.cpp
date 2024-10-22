//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_LINEAR_SPRINGS
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/SPARSE_UNION_FIND.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Math_Tools/Robust_Arithmetic.h>
#include <Core/Math_Tools/sign.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Forces_And_Torques/RIGID_LINEAR_SPRINGS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <cfloat>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_LINEAR_SPRINGS<TV>::
RIGID_LINEAR_SPRINGS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :BASE(rigid_body_collection_input)
{
    Invalidate_CFL();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_LINEAR_SPRINGS<TV>::
~RIGID_LINEAR_SPRINGS()
{
}
//#####################################################################
// Function Set_Restlength_From_Material_Coordinates
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Restlengths()
{
    restlength.Resize(attachment_radius.m);
    Invalidate_CFL();
    for(int i=0;i<segment_mesh.elements.m;i++) restlength(i)=Spring_Length(i);
    visual_restlength=restlength;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    segment_mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated)
{
    Update_Force_Elements(force_segments,segment_mesh.elements,particle_is_simulated);
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(int b,const T overdamping_fraction) // 1 is critically damped
{
    TV X0=Attachment_Location(b,0),X1=Attachment_Location(b,1);
    TV direction=(X1-X0).Normalized();
    RIGID_BODY<TV>& b1=Body(b,0);
    RIGID_BODY<TV>& b2=Body(b,1);

    T harmonic_mass=TV::Dot_Product((b1.Impulse_Factor(X0)+b2.Impulse_Factor(X1))*direction,direction);
    Set_Damping(b,overdamping_fraction*2*sqrt(youngs_modulus(b)*restlength(b)/harmonic_mass));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Update_Position_Based_State(const T time)
{
    states.Resize(segment_mesh.elements.m);
    current_lengths.Resize(segment_mesh.elements.m);
    
    for(int s:force_segments){
        TV X0=Attachment_Location(s,0),X1=Attachment_Location(s,1);
        STATE& state=states(s);
        state.nodes=segment_mesh.elements(s);
        state.direction=X1-X0;
        current_lengths(s)=state.direction.Normalize();
        state.coefficient=damping(s)/restlength(s);
        state.r(0)=X0-Body(s,0).Frame().t;
        state.r(1)=X1-Body(s,1).Frame().t;}
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Force(ARRAY_VIEW<TWIST<TV> > rigid_F,const STATE& state,const TV& force) const
{
    if(!rigid_body_collection.Rigid_Body(state.nodes(0)).Has_Infinite_Inertia()){
        rigid_F(state.nodes(0)).linear+=force;
        rigid_F(state.nodes(0)).angular+=TV::Cross_Product(state.r(0),force);}
    if(!rigid_body_collection.Rigid_Body(state.nodes(1)).Has_Infinite_Inertia()){
        rigid_F(state.nodes(1)).linear-=force;
        rigid_F(state.nodes(1)).angular-=TV::Cross_Product(state.r(1),force);}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    for(int s:force_segments){
        const STATE& state=states(s);
        TV force=youngs_modulus(s)/restlength(s)*(current_lengths(s)-visual_restlength(s))*state.direction;
        Add_Force(rigid_F,state,force);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    //LOG::Time("spring add velocity dependent");
    for(int s:force_segments){
        const STATE& state=states(s);
        VECTOR<TV,2> V;
        for(int i=0;i<2;i++){
            const TWIST<TV>& twist=rigid_V(segment_mesh.elements(s)(i));
            V(i)=twist.linear+TV::Cross_Product(twist.angular,state.r(i));}
        TV force=(state.coefficient*TV::Dot_Product(V(1)-V(0),state.direction))*state.direction;
        Add_Force(rigid_F,state,force);}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time,bool transpose) const
{
    for(int s:force_segments){
        const STATE& state=states(s);
        VECTOR<TV,2> V;
        for(int i=0;i<2;i++){
            const TWIST<TV>& twist=rigid_V(segment_mesh.elements(s)(i));
            V(i)=twist.linear+TV::Cross_Product(twist.angular,state.r(i));}
        TV dl=V(1)-V(0),dl_projected=dl.Projected_On_Unit_Direction(state.direction);
        TV force=youngs_modulus(s)/restlength(s)*dl_projected;
        Add_Force(rigid_F,state,force);}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
//     T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
//     for(int s:force_segments){
//         const VECTOR<int,2>& nodes=segment_mesh.elements(s);
//         T ym=youngs_modulus(s);
//         T d=damping(s);
//         for(int k=0;k<2;k++){
//             frequency(nodes[k]).elastic_squared+=particles.one_over_effective_mass(nodes[k])/restlength(s)*4*ym*one_over_cfl_number_squared;
//             frequency(nodes[k]).damping+=particles.one_over_effective_mass(nodes[k])/restlength(s)*2*d*one_over_cfl_number;}}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
CFL_Strain_Rate() const
{
//     T max_strain_rate=0,strain_rate;TV dx;
//     if(use_rest_state_for_strain_rate) for(int s:force_segments){
//         int i,j;segment_mesh.elements(s).Get(i,j);
//         dx=particles.X(j)-particles.X(i);T magnitude=dx.Magnitude();if(magnitude!=0) dx.Normalize();
//         strain_rate=TV::Dot_Product(particles.V(j)-particles.V(i),dx)/restlength(s);
//         max_strain_rate=max(max_strain_rate,abs(strain_rate));}
//     else for(int s:force_segments){
//         int i,j;segment_mesh.elements(s).Get(i,j);
//         dx=particles.X(j)-particles.X(i);
//         strain_rate=TV::Dot_Product(particles.V(j)-particles.V(i),dx)/TV::Dot_Product(dx,dx);
//         max_strain_rate=max(max_strain_rate,abs(strain_rate));}
//     return Robust_Divide(max_strain_per_time_step,max_strain_rate);
    return FLT_MAX;
}
//#####################################################################
// Function Average_Restlength
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Average_Restlength() const
{
    return restlength.Average();
}
//#####################################################################
// Function Print_Restlength_Statistics
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
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
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Print_Deformation_Statistics() const
{
    LOG::SCOPE scope("linear spring deformation","linear spring deformation");
    ARRAY<T> deformation(segment_mesh.elements.m,no_init);
    for(int s=0;s<segment_mesh.elements.m;s++){
        int i,j;segment_mesh.elements(s).Get(i,j);
        T length=Spring_Length(i),rl=visual_restlength(s);
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
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Maximum_Compression_Or_Expansion_Fraction(int* index) const
{
    T max_compression=0;int max_index=1;
    for(int s=0;s<segment_mesh.elements.m;s++){
        int i,j;segment_mesh.elements(s).Get(i,j);
        T length=Spring_Length(s);
        T rl=visual_restlength(s);
        T compression=(rl)?(abs(length-rl)/rl):((length==0)?0:FLT_MAX);
        if(compression>max_compression){max_compression=compression;max_index=s;}}
    if(index) (*index)=max_index;
    return max_compression;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Potential_Energy(int s,const T time) const
{
    return (T).5*youngs_modulus(s)/restlength(s)*sqr(Spring_Length(s)-visual_restlength(s));
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(int s:force_segments){
        potential_energy+=Potential_Energy(s,time);}
    return potential_energy;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Compute_Total_Energy(const T time) const
{
    T total_energy=0;
    for(int s:force_segments){
        total_energy+=Potential_Energy(s,time);}
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
        total_energy+=rigid_body_collection.Rigid_Body(i).Kinetic_Energy();
    }
    return total_energy;
}
//#####################################################################
// Function Attachment_Location
//#####################################################################
template<class TV> TV RIGID_LINEAR_SPRINGS<TV>::
Attachment_Location(int s,int b) const
{
    return Body(s,b).World_Space_Point(attachment_radius(s)(b));
}
template<class TV> TV RIGID_LINEAR_SPRINGS<TV>::
Endpoint_Velocity(int s,int b) const
{
    return Body(s,b).Pointwise_Object_Velocity(Attachment_Location(s,b));
}
//#####################################################################
// Function Spring_Length
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Spring_Length(int s) const
{
    return (Attachment_Location(s,0)-Attachment_Location(s,1)).Magnitude();
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Stiffness(int b,T stiffness)
{
    if(youngs_modulus.m<=b) youngs_modulus.Resize(segment_mesh.elements.m);
    youngs_modulus(b)=stiffness;
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Damping(int b,T damp)
{
    if(damping.m<=b) damping.Resize(segment_mesh.elements.m);
    damping(b)=damp;
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Restlength
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Restlength(int b,T length,T visual)
{
    if(restlength.m<=b) restlength.Resize(segment_mesh.elements.m);
    if(visual_restlength.m<=b) visual_restlength.Resize(segment_mesh.elements.m);
    restlength(b)=length>=0?length:Spring_Length(b);
    visual_restlength(b)=visual>=0?visual:restlength(b);
    Invalidate_CFL();
}
//#####################################################################
// Function Add_Spring
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Spring(int body0,int body1,const TV& r1,const TV& r2)
{
    attachment_radius.Append(VECTOR<TV,2>(r1,r2));
    segment_mesh.elements.Append(VECTOR<int,2>(body0,body1));
    Invalidate_CFL();
}
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Effective_Impulse_Factor(int s,int b) const
{
    return TV::Dot_Product(Body(s,b).Impulse_Factor(Attachment_Location(s,b))*states(s).direction,states(s).direction);
}
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Effective_Impulse_Factor(int s) const
{
    return Effective_Impulse_Factor(s,0)+Effective_Impulse_Factor(s,1);
}
template<class TV> const RIGID_BODY<TV>& RIGID_LINEAR_SPRINGS<TV>::
Body(int s,int b) const
{
    return rigid_body_collection.Rigid_Body(segment_mesh.elements(s)(b));
}
template<class TV> RIGID_BODY<TV>& RIGID_LINEAR_SPRINGS<TV>::
Body(int s,int b)
{
    return rigid_body_collection.Rigid_Body(segment_mesh.elements(s)(b));
}
//#####################################################################
namespace PhysBAM{
template class RIGID_LINEAR_SPRINGS<VECTOR<float,1> >;
template class RIGID_LINEAR_SPRINGS<VECTOR<float,2> >;
template class RIGID_LINEAR_SPRINGS<VECTOR<float,3> >;
template class RIGID_LINEAR_SPRINGS<VECTOR<double,1> >;
template class RIGID_LINEAR_SPRINGS<VECTOR<double,2> >;
template class RIGID_LINEAR_SPRINGS<VECTOR<double,3> >;
}
