//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES(DEFORMABLE_PARTICLES<TV>& particles,
    DEFORMABLE_PARTICLES<TV>& undeformed_particles,T_OBJECT& collision_body,
    T_SURFACE& undeformed_triangulated_surface,IMPLICIT_OBJECT<TV>& implicit_surface,
    T stiffness,T separation_parameter)
    :COLLISION_FORCE<TV>(particles),undeformed_particles(undeformed_particles),collision_body(collision_body),
    undeformed_triangulated_surface(undeformed_triangulated_surface),
    triangulated_surface(collision_body.Get_Boundary_Object()),implicit_surface(implicit_surface),
    stiffness(stiffness),separation_parameter(separation_parameter),pe(0)
{
    triangulated_surface.Initialize_Hierarchy();
    undeformed_triangulated_surface.Initialize_Hierarchy();
    collision_body.Initialize_Hierarchy();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
~DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES()
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Closest_Surface_Triangle
//#####################################################################
template<class TV> int DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Estimate_Closest_Undeformed_Surface_Triangle(const TV& X,int p) const
{
    ARRAY<int> nearby_surface_triangles;
    TV surface_location=implicit_surface.Closest_Point_On_Boundary(X);
    for(T thickness=separation_parameter;;thickness*=(T)2){
        undeformed_triangulated_surface.hierarchy->Intersection_List(surface_location,nearby_surface_triangles,thickness);
        for(int i=0;i<nearby_surface_triangles.m;i++){
            int ct=nearby_surface_triangles(i);
            if(!triangulated_surface.mesh.elements(ct).Contains(p))
                return ct;}}
}
//#####################################################################
// Function Update_Penetrating Particles
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Penetrating_Particles(int p)
{
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX T_SIMPLEX;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX T_SURFACE_SIMPLEX;
    const auto& mesh=collision_body.mesh;
    TV x=particles.X(p);
    ARRAY<int> intersection_list;
    ARRAY<int> particles_to_ignore;
    int num_collisions=0;
    T closest_distance_squared=FLT_MAX;
    TV_INT surface_nodes;
    collision_body.hierarchy->Intersection_List(x,intersection_list);
    for(int n=0;n<intersection_list.m;n++){
        bool ignore=false;
        int t=intersection_list(n);
        const VECTOR<int,TV::m+1>& nodes=mesh.elements(t);
        if(nodes.Contains(p))
            continue;
        for(int q=0;q<particles_to_ignore.m;q++)
            if(nodes.Contains(particles_to_ignore(q))){
                ignore=true;
                break;}
        if(ignore) continue;
        const VECTOR<T,TV::m+1>& w=T_SIMPLEX::Barycentric_Coordinates(x,particles.X.Subset(nodes));
        if(w.Min()<-1e-20) continue;
        num_collisions++;
        if(num_collisions==1){
            ARRAY<int> nearby_surface_triangles;
            int ct=closest_surface_triangle(p);
            if(ct==-1) ct=Estimate_Closest_Undeformed_Surface_Triangle(undeformed_particles.X.Subset(nodes).Weighted_Sum(w),p);
            TV projected_point=undeformed_triangulated_surface.Get_Element(ct).Closest_Point(x);
            closest_distance_squared=(x-projected_point).Magnitude_Squared();
            T closest_distance_upper_bound=sqrt(closest_distance_squared);
            triangulated_surface.hierarchy->Intersection_List(x,nearby_surface_triangles,1.05*closest_distance_upper_bound);
            PHYSBAM_ASSERT(nearby_surface_triangles.m!=0);
            for(int k=0;k<nearby_surface_triangles.m;k++){
                int tri=nearby_surface_triangles(k);
                if(triangulated_surface.mesh.elements(tri).Contains(p)) continue;
                TV new_point=undeformed_triangulated_surface.Get_Element(tri).Closest_Point(x);
                T new_distance_squared=(x-new_point).Magnitude_Squared();
                if(new_distance_squared<closest_distance_squared){
                    closest_distance_squared=new_distance_squared;
                    ct=tri;}}
            closest_surface_triangle(p)=ct;
            surface_nodes=triangulated_surface.mesh.elements(ct);
            const ARRAY<TV_INT>& extra=extra_surface_triangles(p);
            for(int k=0;k<extra.m;k++){
                T_SURFACE_SIMPLEX t(particles.X.Subset(extra(k)));
                TV new_point=t.Closest_Point(x);
                T new_distance_squared=(x-new_point).Magnitude_Squared();
                if(new_distance_squared<closest_distance_squared){
                    closest_distance_squared=new_distance_squared;
                    surface_nodes=extra(k);}}}
        particles_to_ignore.Append_Elements(nodes.Remove_Index(w.Arg_Min()));
    }
    for(int i=0;i<num_collisions;i++)
        penetrating_particles.Append(surface_nodes.Insert(p,0));
}
//#####################################################################
// Function Update_Surface_Triangles
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Surface_Triangles()
{
    penetrating_particles.Remove_All();
    collision_body.hierarchy->Update_Boxes(separation_parameter);
    triangulated_surface.hierarchy->Update_Boxes(separation_parameter);
    if(colliding_particles.m)
        for(int pp=0;pp<colliding_particles.m;pp++){
            int p=colliding_particles(pp);
            Update_Penetrating_Particles(p);}
    else
        for(int p=0;p<particles.X.m;p++)
            Update_Penetrating_Particles(p);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    closest_surface_triangle.Resize(particles.X.m,true,true,-1);
    if(extra_surface_triangles.m!=particles.X.m){
        extra_surface_triangles.Resize(particles.X.m);
        ARRAY<bool> surface_node(particles.X.m);
        surface_node.Subset(triangulated_surface.mesh.elements.Flattened()).Fill(true);
        for(int i=0;i<collision_body.mesh.elements.m;i++){
            VECTOR<int,TV::m+1> nodes=collision_body.mesh.elements(i);
            for(int j=0;j<TV::m+1;j++)
                if(surface_node(nodes(j)))
                    extra_surface_triangles(nodes(j)).Append(nodes.Remove_Index(j));}}

    pe=0;
    grad_pe.Remove_All();
    H_pe.Remove_All();
    stored_weights.Remove_All();
    Update_Surface_Triangles();
    for(int pp=0;pp<penetrating_particles.m;pp++)
        Update_Position_Based_State_Particle(pp);
}
//#####################################################################
// Function Penalty
//#####################################################################
template<class T,class TV> void
Penalty(VECTOR<int,4> nodes,const INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,TV::m+1>& >&X,T& e,VECTOR<TV,TV::m+1>& de,VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1>& he,ARRAY<VECTOR<T,TV::m+1> >& stored_weights,T stiffness)
{
    auto w=Hess_From_Var<3,0>(X(nodes(0))-X(nodes(3)));
    auto v=Hess_From_Var<3,1>(X(nodes(1))-X(nodes(3)));
    auto u=Hess_From_Var<3,2>(X(nodes(2))-X(nodes(3)));
    auto uu=u.Dot(u);
    auto vv=v.Dot(v);
    auto uv=u.Dot(v);
    auto uw=u.Dot(w);
    auto vw=v.Dot(w);
    auto d=uu*vv-sqr(uv);
    auto g=vv*uw-uv*vw;
    auto h=uu*vw-uv*uw;
    auto a=g/d;
    auto b=h/d;
    auto z=a*u+b*v-w;
    auto phi_sq=z.Magnitude_Squared();
    auto ee=stiffness*phi_sq*sqrt(phi_sq+1e-15);
    e=ee.x;
    for(int i=0;i<3;i++){
        TV t=ee.dx(i);
        de(nodes(i))=t;
        de(nodes(3))-=t;}
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            MATRIX<T,TV::m> t=ee.ddx(i,j);
            he(nodes(i))(nodes(j))=t;
            he(nodes(i))(nodes(3))-=t;
            he(nodes(3))(nodes(3))+=t;}
        he(nodes(3))(nodes(i))=he(nodes(i))(nodes(3)).Transposed();}
    VECTOR<T,TV::m+1> weight;
    weight(nodes(0))=-1;
    weight(nodes(1))=b.x;
    weight(nodes(2))=a.x;
    weight(nodes(3))=1-a.x-b.x;
    stored_weights.Append(weight);
}
//#####################################################################
// Function Penalty
//#####################################################################
template<class T,class TV> void
Penalty(VECTOR<int,3> nodes,const INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,TV::m+1>& >&X,T& e,VECTOR<TV,TV::m+1>& de,VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1>& he,ARRAY<VECTOR<T,TV::m+1> >& stored_weights,T stiffness)
{
    auto w=Hess_From_Var<2,0>(X(nodes(0))-X(nodes(2)));
    auto v=Hess_From_Var<2,1>(X(nodes(1))-X(nodes(2)));
    auto vv=v.Dot(v);
    auto vw=v.Dot(w);
    auto a=min(max(vw/vv,0),1);
    auto z=a*v-w;
    auto phi_sq=z.Magnitude_Squared();
    auto ee=stiffness*phi_sq*sqrt(phi_sq+1e-15);
    e=ee.x;
    for(int i=0;i<2;i++){
        TV t=ee.dx(i);
        de(nodes(i))=t;
        de(nodes(2))-=t;}
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            MATRIX<T,TV::m> t=ee.ddx(i,j);
            he(nodes(i))(nodes(j))=t;
            he(nodes(i))(nodes(2))-=t;
            he(nodes(2))(nodes(2))+=t;}
        he(nodes(2))(nodes(i))=he(nodes(i))(nodes(2)).Transposed();}
    VECTOR<T,TV::m+1> weight;
    weight(nodes(0))=-1;
    weight(nodes(1))=a.x;
    weight(nodes(2))=1-a.x;
    stored_weights.Append(weight);
}
//#####################################################################
// Function Penalty
//#####################################################################
template<class T,class TV> void
Penalty(VECTOR<int,2> nodes,const INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,TV::m+1>& >&X,T& e,VECTOR<TV,TV::m+1>& de,VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1>& he,ARRAY<VECTOR<T,TV::m+1> >& stored_weights,T stiffness)
{
    auto w=Hess_From_Var<1,0>(X(nodes(0))-X(nodes(1)));
    auto phi_sq=w.Magnitude_Squared();
    auto ee=stiffness*phi_sq*sqrt(phi_sq+1e-15);
    e=ee.x;
    TV t=ee.dx(0);
    de(nodes(0))=t;
    de(nodes(1))=-t;
    MATRIX<T,TV::m> m=ee.ddx(0,0);
    he(nodes(0))(nodes(0))=m;
    he(nodes(0))(nodes(1))=-m;
    he(nodes(1))(nodes(1))=m;
    he(nodes(1))(nodes(0))=-m.Transposed();
    VECTOR<T,TV::m+1> weight;
    weight(nodes(0))=-1;
    weight(nodes(1))=1;
    stored_weights.Append(weight);
}
//#####################################################################
// Function Update_Position_Based_State_Particle
//#####################################################################
template<class T,class TV> void
Update_Position_Based_State_Particle_Helper(const VECTOR<int,4>& nodes,ARRAY_VIEW<TV> particles_X,
    T& e,VECTOR<TV,TV::m+1>& de,VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1>& he,ARRAY<VECTOR<T,TV::m+1> >& stored_weights,T stiffness)
{
    TV location=particles_X(nodes(0));
    TRIANGLE_3D<T> t(particles_X.Subset(nodes.Remove_Index(0)));
    TV weights=t.Barycentric_Coordinates(location);
    const INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,TV::m+1>& > X=particles_X.Subset(nodes);
    if(weights.x<0){
        T a12=SEGMENT_3D<T>::Interpolation_Fraction(location,t.X.y,t.X.z); // Check edge X.y--X.z
        if(a12<0){
            if(weights.z<0) Penalty(VECTOR<int,3>(0,1,2),X,e,de,he,stored_weights,stiffness); // Closest point is on edge X.x--X.y
            else Penalty(VECTOR<int,2>(0,2),X,e,de,he,stored_weights,stiffness);} // Closest point is t.X.y
        else if(a12>1){
            if(weights.y<0) Penalty(VECTOR<int,3>(0,1,3),X,e,de,he,stored_weights,stiffness); // Closest point is on edge t.X.x--t.X.z
            else Penalty(VECTOR<int,2>(0,3),X,e,de,he,stored_weights,stiffness);} // Closest point is t.X.z
        else Penalty(VECTOR<int,3>(0,2,3),X,e,de,he,stored_weights,stiffness);} // Closest point is on edge t.X.y--t.X.z
    else if(weights.y<0){
        T a02=SEGMENT_3D<T>::Interpolation_Fraction(location,t.X.x,t.X.z); // Check edge t.X.x--t.X.z
        if(a02<0){
            if(weights.z<0) Penalty(VECTOR<int,3>(0,1,2),X,e,de,he,stored_weights,stiffness); // Closest point is on edge t.X.x--t.X.y
            else Penalty(VECTOR<int,2>(0,1),X,e,de,he,stored_weights,stiffness);} // Closest point is t.X.x
        else if(a02>1) Penalty(VECTOR<int,2>(0,3),X,e,de,he,stored_weights,stiffness); // Closest point is t.X.z
        else Penalty(VECTOR<int,3>(0,1,3),X,e,de,he,stored_weights,stiffness);} // Closest point is on edge t.X.x--t.X.z
    else if(weights.z<0) Penalty(VECTOR<int,3>(0,1,2),X,e,de,he,stored_weights,stiffness); // Closest point is on edge t.X.x--t.X.y
    else Penalty(VECTOR<int,4>(0,1,2,3),X,e,de,he,stored_weights,stiffness);
}
//#####################################################################
// Function Update_Position_Based_State_Particle
//#####################################################################
template<class T,class TV> void
Update_Position_Based_State_Particle_Helper(const VECTOR<int,3>& nodes,ARRAY_VIEW<TV> particles_X,
    T& e,VECTOR<TV,TV::m+1>& de,VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1>& he,ARRAY<VECTOR<T,TV::m+1> >& stored_weights,T stiffness)
{
    TV location=particles_X(nodes(0));
    SEGMENT_2D<T> t(particles_X.Subset(nodes.Remove_Index(0)));
    TV weights=t.Barycentric_Coordinates(location);
    const INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,TV::m+1>& > X=particles_X.Subset(nodes);
    if(weights.x<0) Penalty(VECTOR<int,2>(0,1),X,e,de,he,stored_weights,stiffness); // Closest point is t.X.x
    else if(weights.y<0) Penalty(VECTOR<int,2>(0,2),X,e,de,he,stored_weights,stiffness); // Closest point is t.X.y
    else Penalty(VECTOR<int,3>(0,1,2),X,e,de,he,stored_weights,stiffness); // Closest point is on edge X.x--X.y
}
//#####################################################################
// Function Update_Position_Based_State_Particle
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Position_Based_State_Particle(int pp)
{
    T e;
    VECTOR<TV,TV::m+1> de;
    VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1> he;
    Update_Position_Based_State_Particle_Helper(penetrating_particles(pp),particles.X,e,de,he,stored_weights,stiffness);
    pe+=e;
    grad_pe.Append(de);
    H_pe.Append(he);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int pp=0;pp<penetrating_particles.m;pp++)
        F.Subset(penetrating_particles(pp))-=grad_pe(pp);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(int pp=0;pp<penetrating_particles.m;pp++){
        const VECTOR<int,TV::m+1>& n=penetrating_particles(pp);
        for(int i=0;i<n.m;i++){
            int p=n(i);
            for(int j=0;j<n.m;j++)
                F(p)-=H_pe(pp)(i)(j)*V(n(j));}}
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Potential_Energy(const T time) const
{
    return pe;
}
//#####################################################################
// Function Apply_Friction
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Apply_Friction(ARRAY_VIEW<TV> V,const T time) const
{
    for(int pp=0;pp<penetrating_particles.m;pp++){
        VECTOR<int,TV::m+1> nodes=penetrating_particles(pp);
        TV normal=grad_pe(pp)(0);
        T normal_force=normal.Normalize();
        const VECTOR<T,TV::m+1>& weights=stored_weights(pp);
        TV v_hat=V.Subset(nodes).Weighted_Sum(weights);
        T mass_hat=(particles.one_over_mass.Subset(nodes)*weights).Dot(weights);
        TV force_dir=-1/mass_hat*v_hat.Projected_Orthogonal_To_Unit_Direction(normal);
        T force_mag=force_dir.Normalize();
        TV force=min(force_mag,normal_force*coefficient_of_friction)*force_dir;
        for(int i=0;i<TV::m+1;i++) V(nodes(i))+=weights(i)*particles.one_over_mass(nodes(i))*force;}
}
template class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<float,2> >;
template class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<float,3> >;
template class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<double,2> >;
template class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<double,3> >;
}
