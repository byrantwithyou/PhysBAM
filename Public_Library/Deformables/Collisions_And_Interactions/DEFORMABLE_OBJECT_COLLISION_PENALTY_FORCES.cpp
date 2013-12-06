//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{
template<class T> static TRIANGULATED_SURFACE<T>& Triangulated_Surface_Helper(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    if(!tetrahedralized_volume.triangulated_surface) tetrahedralized_volume.Initialize_Triangulated_Surface();
    return *tetrahedralized_volume.triangulated_surface;
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
    DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES(DEFORMABLE_PARTICLES<TV>& particles,DEFORMABLE_PARTICLES<TV>& undeformed_particles,TETRAHEDRALIZED_VOLUME<T>& collision_body,TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface,IMPLICIT_OBJECT<TV>& implicit_surface,
        T stiffness,T separation_parameter)
    :DEFORMABLES_FORCES<TV>(particles),undeformed_particles(undeformed_particles),collision_body(collision_body),
    undeformed_triangulated_surface(undeformed_triangulated_surface),triangulated_surface(Triangulated_Surface_Helper(collision_body)),
    implicit_surface(implicit_surface),closest_surface_triangle(particles.X.m,true,-1),extra_surface_triangles(particles.X.m),stiffness(stiffness),separation_parameter(separation_parameter),pe(0)
{
    triangulated_surface.Update_Triangle_List();
    triangulated_surface.Initialize_Hierarchy();
    undeformed_triangulated_surface.Update_Triangle_List();
    undeformed_triangulated_surface.Initialize_Hierarchy();
    triangulated_surface.avoid_normal_interpolation_across_sharp_edges=false;
    collision_body.Initialize_Hierarchy();
    ARRAY<bool> surface_node(particles.X.m);
    surface_node.Subset(triangulated_surface.mesh.elements.Flattened()).Fill(true);
    for(int i=0;i<collision_body.mesh.elements.m;i++){
        VECTOR<int,4> nodes=collision_body.mesh.elements(i);
        for(int j=0;j<4;j++)
            if(surface_node(nodes(j)))
                extra_surface_triangles(nodes(j)).Append(nodes.Remove_Index(j));}
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
    const TETRAHEDRON_MESH& tetrahedron_mesh=collision_body.mesh;
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
        const VECTOR<int,4>& nodes=tetrahedron_mesh.elements(t);
        if(nodes.Contains(p))
            continue;
        for(int q=0;q<particles_to_ignore.m;q++)
            if(nodes.Contains(particles_to_ignore(q))){
                ignore=true;
                break;}
        if(ignore) continue;
        const VECTOR<T,TV::m+1>& w=TETRAHEDRON<T>::Barycentric_Coordinates(x,particles.X.Subset(nodes));
        if(w.Min()<-1e-20) continue;
        num_collisions++;
        TV weights;
        if(num_collisions==1){
            ARRAY<int> nearby_surface_triangles;
            int ct=closest_surface_triangle(p);
            if(ct==-1) ct=Estimate_Closest_Undeformed_Surface_Triangle(undeformed_particles.X.Subset(nodes).Weighted_Sum(w),p);
            TV projected_point=(*triangulated_surface.triangle_list)(ct).Closest_Point(x,weights);
            closest_distance_squared=(x-projected_point).Magnitude_Squared();
            T closest_distance_upper_bound=sqrt(closest_distance_squared);
            triangulated_surface.hierarchy->Intersection_List(x,nearby_surface_triangles,1.05*closest_distance_upper_bound);
            PHYSBAM_ASSERT(nearby_surface_triangles.m!=0);
            for(int k=0;k<nearby_surface_triangles.m;k++){
                int tri=nearby_surface_triangles(k);
                if(triangulated_surface.mesh.elements(tri).Contains(p)) continue;
                TV new_point=(*triangulated_surface.triangle_list)(tri).Closest_Point(x,weights);
                T new_distance_squared=(x-new_point).Magnitude_Squared();
                if(new_distance_squared<closest_distance_squared){
                    closest_distance_squared=new_distance_squared;
                    ct=tri;}}
            closest_surface_triangle(p)=ct;
            surface_nodes=triangulated_surface.mesh.elements(ct);
            const ARRAY<TV_INT>& extra=extra_surface_triangles(p);
            for(int k=0;k<extra.m;k++){
                TRIANGLE_3D<T> t(particles.X.Subset(extra(k)));
                TV new_point=t.Closest_Point(x,weights);
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
    triangulated_surface.Update_Triangle_List();
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
Update_Position_Based_State(const T time,const bool is_position_update)
{
    pe=0;
    grad_pe.Remove_All();
    H_pe.Remove_All();
    Update_Surface_Triangles();
    for(int pp=0;pp<penetrating_particles.m;pp++)
        Update_Position_Based_State_Particle(pp);
}
//#####################################################################
// Function Penalty
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Penalty(VECTOR<int,4> nodes, const INDIRECT_ARRAY<ARRAY_VIEW<TV, int>, VECTOR<int,4>& >&X, T& e, VECTOR<TV,4>& de, VECTOR<VECTOR<MATRIX<T,TV::m>,4>,4>& he)
{
    auto w=From_Var<3,0>(X(nodes(0))-X(nodes(3)));
    auto v=From_Var<3,1>(X(nodes(1))-X(nodes(3)));
    auto u=From_Var<3,2>(X(nodes(2))-X(nodes(3)));
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
    auto ee=stiffness*phi_sq;
    e=ee.x;
    for(int i=0;i<3;i++){
        VECTOR<T,3> t=ee.dx(i);
        de(nodes(i))=t;
        de(nodes(3))-=t;}
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            MATRIX<T,3> t=ee.ddx(i,j);
            he(nodes(i))(nodes(j))=t;
            he(nodes(i))(nodes(3))-=t;
            he(nodes(3))(nodes(3))+=t;}
        he(nodes(3))(nodes(i))=he(nodes(i))(nodes(3)).Transposed();}
}
//#####################################################################
// Function Penalty
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Penalty(VECTOR<int,3> nodes, const INDIRECT_ARRAY<ARRAY_VIEW<TV, int>, VECTOR<int,4>& >&X, T& e, VECTOR<TV,4>& de, VECTOR<VECTOR<MATRIX<T,TV::m>,4>,4>& he)
{
    auto w=From_Var<2,0>(X(nodes(0))-X(nodes(2)));
    auto v=From_Var<2,1>(X(nodes(1))-X(nodes(2)));
    auto vv=v.Dot(v);
    auto vw=v.Dot(w);
    auto a=min(max(vw/vv,0),1);
    auto z=a*v-w;
    auto phi_sq=z.Magnitude_Squared();
    auto ee=stiffness*phi_sq;
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
}
//#####################################################################
// Function Penalty
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Penalty(VECTOR<int,2> nodes, const INDIRECT_ARRAY<ARRAY_VIEW<TV, int>, VECTOR<int,4>& >&X, T& e, VECTOR<TV,4>& de, VECTOR<VECTOR<MATRIX<T,TV::m>,4>,4>& he)
{
    auto w=From_Var<1,0>(X(nodes(0))-X(nodes(1)));
    auto phi_sq=w.Magnitude_Squared();
    auto ee=stiffness*phi_sq;
    e=ee.x;
    TV t=ee.dx(0);
    de(nodes(0))=t;
    de(nodes(1))=-t;
    MATRIX<T,TV::m> m=ee.ddx(0,0);
    he(nodes(0))(nodes(0))=m;
    he(nodes(0))(nodes(1))=-m;
    he(nodes(1))(nodes(1))=m;
    he(nodes(1))(nodes(0))=-m.Transposed();
}
//#####################################################################
// Function Update_Position_Based_State_Particle
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Position_Based_State_Particle(int pp)
{
    const VECTOR<int,TV::m+1>& nodes=penetrating_particles(pp);
    T e;
    VECTOR<TV,TV::m+1> de;
    VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1> he;
    TV location=particles.X(nodes(0));
    TRIANGLE_3D<T> t(particles.X.Subset(nodes.Remove_Index(0)));
    TV weights=t.Barycentric_Coordinates(location);
    const INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,4>& > X=particles.X.Subset(nodes);
    if(weights.x<0){
        T a23=SEGMENT_3D<T>::Interpolation_Fraction(location,t.X.y,t.X.z); // Check edge X.y--X.z
        if(a23<0){
            if(weights.z<0) Penalty(VECTOR<int,3>(0,1,2),X,e,de,he); // Closest point is on edge X.x--X.y
            else Penalty(VECTOR<int,2>(0,2),X,e,de,he);} // Closest point is t.X.y
        else if(a23>1){
            if(weights.y<0) Penalty(VECTOR<int,3>(0,1,3),X,e,de,he); // Closest point is on edge t.X.x--t.X.z
            else Penalty(VECTOR<int,2>(0,3),X,e,de,he);} // Closest point is t.X.z
        else Penalty(VECTOR<int,3>(0,2,3),X,e,de,he);} // Closest point is on edge t.X.y--t.X.z
    else if(weights.y<0){
        T a13=SEGMENT_3D<T>::Interpolation_Fraction(location,t.X.x,t.X.z); // Check edge t.X.x--t.X.z
        if(a13<0){
            if(weights.z<0) Penalty(VECTOR<int,3>(0,1,2),X,e,de,he); // Closest point is on edge t.X.x--t.X.y
            else Penalty(VECTOR<int,2>(0,1),X,e,de,he);} // Closest point is t.X.x
        else if(a13>1) Penalty(VECTOR<int,2>(0,3),X,e,de,he); // Closest point is t.X.z
        else Penalty(VECTOR<int,3>(0,1,3),X,e,de,he);} // Closest point is on edge t.X.x--t.X.z
    else if(weights.z<0) Penalty(VECTOR<int,3>(0,1,2),X,e,de,he); // Closest point is on edge t.X.x--t.X.y
    else Penalty(VECTOR<int,4>(0,1,2,3),X,e,de,he);
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
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const
{
    for(int pp=0;pp<penetrating_particles.m;pp++){
        const VECTOR<int,TV::m+1>& n=penetrating_particles(pp);
        for(int i=0;i<n.m;i++){
            int p=n(i);
            for(int j=0;j<n.m;j++)
                F(p)-=H_pe(pp)(i)(j)*V(n(j))*scale;}}
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
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
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
template class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<float,3> >;
template class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<double,3> >;
}
