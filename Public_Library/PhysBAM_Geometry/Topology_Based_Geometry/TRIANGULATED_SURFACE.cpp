//#####################################################################
// Copyright 2002-2008, Chris Allocco, Robert Bridson, Kevin Der, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Neil Molino, Igor Neverov, Andrew Selle, Eftychios Sifakis, Jonathan Su, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGULATED_SURFACE
//##################################################################### 
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SPHERE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_3D_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>::
TRIANGULATED_SURFACE()
    :MESH_OBJECT<TV,TRIANGLE_MESH>(*new TRIANGLE_MESH,*new GEOMETRY_PARTICLES<TV>),
    triangle_list(0),segment_lengths(0),hierarchy(0),particle_hierarchy(0),
    avoid_normal_interpolation_across_sharp_edges(false),normal_variance_threshold((T).1),
    vertex_normals(0),face_vertex_normals(0)
{
    Use_Face_Normals();
    this->need_destroy_mesh=this->need_destroy_particles=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>::
TRIANGULATED_SURFACE(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,TRIANGLE_MESH>(triangle_mesh_input,particles_input),
    triangle_list(0),segment_lengths(0),hierarchy(0),particle_hierarchy(0),
    avoid_normal_interpolation_across_sharp_edges(false),normal_variance_threshold((T).1),
    vertex_normals(0),face_vertex_normals(0)
{
    Use_Face_Normals();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>::
~TRIANGULATED_SURFACE()
{
    delete triangle_list;delete hierarchy;delete segment_lengths;delete particle_hierarchy;
    delete vertex_normals;delete face_vertex_normals;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Clean_Memory()
{
    MESH_OBJECT<TV,TRIANGLE_MESH>::Clean_Memory();
    delete triangle_list;triangle_list=0;delete hierarchy;hierarchy=0;
    delete segment_lengths;segment_lengths=0;delete particle_hierarchy;particle_hierarchy=0;delete vertex_normals;vertex_normals=0;
    delete face_vertex_normals;face_vertex_normals=0;
    Use_Face_Normals();
}
//#####################################################################
// Function Refresh_Auxiliary_Structures_Helper
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Refresh_Auxiliary_Structures_Helper()
{
    if(triangle_list) Update_Triangle_List();
    if(hierarchy) Initialize_Hierarchy();
    if(vertex_normals || face_vertex_normals) Update_Vertex_Normals();
    if(segment_lengths) Initialize_Segment_Lengths();
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Hierarchy(const bool update_boxes,const int triangles_per_group) // creates and updates the boxes as well
{
    delete hierarchy;
    if(triangle_list) hierarchy=new TRIANGLE_HIERARCHY<T>(mesh,particles,*triangle_list,update_boxes,triangles_per_group);
    else hierarchy=new TRIANGLE_HIERARCHY<T>(mesh,particles,update_boxes,triangles_per_group);
}
//#####################################################################
// Function Initialize_Particle_Hierarchy
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Particle_Hierarchy(const INDIRECT_ARRAY<ARRAY_VIEW<TV> >& particle_subset_input,const bool update_boxes,const int particles_per_group) // creates and updates the boxes as well
{
    typedef VECTOR<T,3> TV;
    delete particle_hierarchy;
    particle_hierarchy=new PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >(particle_subset_input,update_boxes,particles_per_group);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Rescale(const T scaling_x,const T scaling_y,const T scaling_z)
{
    typedef VECTOR<T,3> TV;
    if(scaling_x*scaling_y*scaling_z<=0) PHYSBAM_FATAL_ERROR();
    for(int k=0;k<particles.Size();k++) particles.X(k)*=TV(scaling_x,scaling_y,scaling_z);
    if(triangle_list) Update_Triangle_List();if(hierarchy) hierarchy->Update_Boxes();if(bounding_box) Update_Bounding_Box();
}
//#####################################################################
// Function Update_Triangle_List
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Update_Triangle_List()
{
    Update_Triangle_List(particles.X);
}
//#####################################################################
// Function Update_Triangle_List
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Update_Triangle_List(ARRAY_VIEW<const TV> X)
{
    if(!triangle_list) triangle_list=new ARRAY<TRIANGLE_3D<T> >;
    triangle_list->Resize(mesh.elements.m);
    for(int t=0;t<mesh.elements.m;t++)
        (*triangle_list)(t)=TRIANGLE_3D<T>(X.Subset(mesh.elements(t)));
}
//#####################################################################
// Function Initialize_Torus_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Torus_Mesh_And_Particles(const int m,const int n,const T major_radius,const T minor_radius)
{
    typedef VECTOR<T,3> TV;
    T di=T(2*pi)/m,dj=T(2*pi)/n;
    for(int j=0;j<n;j++){
        T phi=-dj*j,radius=major_radius+minor_radius*cos(phi),z=minor_radius*sin(phi);
        for(int i=0;i<m;i++){
            int p=particles.Add_Element();T theta=di*(i-(T).5*(j&1));
            particles.X(p)=TV(radius*cos(theta),radius*sin(theta),z);}}
    mesh.Initialize_Torus_Mesh(m,n);
}
//#####################################################################
// Function Initialize_Cylinder_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Cylinder_Mesh_And_Particles(const int m,const int n,const T length,const T radius,const bool create_caps)
{
    typedef VECTOR<T,3> TV;
    particles.Delete_All_Elements();T dtheta=(T)two_pi/n;T dlength=length/(m-1);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        int p=particles.Add_Element();T theta=j*dtheta;
        particles.X(p)=TV(dlength*i,radius*sin(theta),radius*cos(theta));}
    if(create_caps){int p_1=particles.Add_Element();int p_2=particles.Add_Element();particles.X(p_1)=TV(0,0,0);particles.X(p_2)=TV(length,0,0);}
    mesh.Initialize_Cylinder_Mesh(m,n,create_caps);
}
//#####################################################################
// Function Initialize_Segment_Lengths
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Initialize_Segment_Lengths()
{
    bool segment_mesh_defined=mesh.segment_mesh!=0;if(!segment_mesh_defined) mesh.Initialize_Segment_Mesh();
    delete segment_lengths;segment_lengths=new ARRAY<T>(mesh.segment_mesh->elements.m);
    for(int t=0;t<mesh.segment_mesh->elements.m;t++) 
        (*segment_lengths)(t)=(particles.X(mesh.segment_mesh->elements(t)(0))-particles.X(mesh.segment_mesh->elements(t)(1))).Magnitude();
    if(!segment_mesh_defined){delete mesh.segment_mesh;mesh.segment_mesh=0;}
}
//#####################################################################
// Function Update_Vertex_Normals
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Update_Vertex_Normals()
{  
    typedef VECTOR<T,3> TV;
    bool incident_elements_defined=mesh.incident_elements!=0;if(!incident_elements_defined) mesh.Initialize_Incident_Elements();
    bool triangle_list_defined=triangle_list!=0;if(!triangle_list_defined) Update_Triangle_List();

    if(avoid_normal_interpolation_across_sharp_edges){
        delete vertex_normals;vertex_normals=0;
        if(!face_vertex_normals) face_vertex_normals=new ARRAY<VECTOR<TV,3> >(mesh.elements.m);else face_vertex_normals->Resize(mesh.elements.m);
        TV face_normal,normal,zero_vector(0,0,0);
        for(int t=0;t<mesh.elements.m;t++){
            face_normal=Face_Normal(t);
            for(int i=0;i<3;i++){
                normal=zero_vector;int node=mesh.elements(t)(i); 
                for(int k=0;k<(*mesh.incident_elements)(node).m;k++){
                    TRIANGLE_3D<T>& triangle=(*triangle_list)((*mesh.incident_elements)(node)(k));
                    TV local_normal=triangle.Normal();
                    if((face_normal-local_normal).Magnitude_Squared() < normal_variance_threshold) normal+=triangle.Area()*local_normal;}
                if(normal != zero_vector) normal.Normalize();
                (*face_vertex_normals)(t)(i)=normal;}}}
    else{
        delete face_vertex_normals;face_vertex_normals=0;
        if(!vertex_normals) vertex_normals=new ARRAY<TV>(mesh.number_nodes);
        else vertex_normals->Resize(mesh.number_nodes);
        for(int k=0;k<vertex_normals->m;k++){
            (*vertex_normals)(k)=TV(); // initialize to zero
            for(int kk=0;kk<(*mesh.incident_elements)(k).m;kk++) 
                (*vertex_normals)(k)+=(*triangle_list)((*mesh.incident_elements)(k)(kk)).Area()*(*triangle_list)((*mesh.incident_elements)(k)(kk)).Normal();
            if((*vertex_normals)(k) != TV()) (*vertex_normals)(k).Normalize();}}

    if(!triangle_list_defined){delete triangle_list;triangle_list=0;}
    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> TRIANGULATED_SURFACE<T>::
Normal(const TV& location,const int aggregate) const 
{
    typedef VECTOR<T,3> TV;
    assert(aggregate >= 1);

    if(use_vertex_normals){
        TV normal1,normal2,normal3;
        int node1,node2,node3;mesh.elements(aggregate).Get(node1,node2,node3);
        if(avoid_normal_interpolation_across_sharp_edges) (*face_vertex_normals)(aggregate).Get(normal1,normal2,normal3);
        else{normal1=(*vertex_normals)(node1);normal2=(*vertex_normals)(node2);normal3=(*vertex_normals)(node3);}
        TV weights=TRIANGLE_3D<T>::Barycentric_Coordinates(location,particles.X(node1),particles.X(node2),particles.X(node3));
        return (weights.x*normal1+weights.y*normal2+weights.z*normal3).Normalized();}
    else return Face_Normal(aggregate);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Inside(const TV& location,const T thickness_over_two) const
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(bounding_box && hierarchy);
    if(bounding_box->Outside(location,thickness_over_two)) return false;
    if(hierarchy->box_hierarchy(hierarchy->root).Outside(location,thickness_over_two)) return false;
    RAY<TV> ray(location,TV(0,0,1),true);
    return Inside_Using_Ray_Test(ray,thickness_over_two);
}
//#####################################################################
// Function Inside_Relative_To_Triangle
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Inside_Relative_To_Triangle(const TV& location,const int triangle_index_for_ray_test,const T thickness_over_two) const
{
    PHYSBAM_ASSERT(triangle_list);
    RAY<TV> ray(SEGMENT_3D<T>(location,(*triangle_list)(triangle_index_for_ray_test).Center()));
    ray.t_max+=2*thickness_over_two;
    return Inside_Using_Ray_Test(ray,thickness_over_two); 
}
//#####################################################################
// Function Inside_Using_Ray_Test
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Inside_Using_Ray_Test(RAY<TV>& ray,const T thickness_over_two) const
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(mesh.adjacent_elements && triangle_list);
    bool inside=false;
    if(INTERSECTION::Intersects(ray,*this,thickness_over_two) && ray.t_max > 0){ // otherwise missed the object or on the boundary
        TV point=ray.Point(ray.t_max); // point is inside if and only if location is inside
        T thickness=2*thickness_over_two;
        TRIANGLE_3D<T>& triangle=(*triangle_list)(ray.aggregate_id);
        int region_id,region=triangle.Region(point,region_id,thickness);
        if(region==0){ // vertex
            if(Signed_Solid_Angle_Of_Triangle_Web(point,mesh.elements(ray.aggregate_id)(region_id)) > 0) inside=true;} 
        else if(region==1) { // edge
            int node1,node2,node3,neighbor=-1;mesh.elements(ray.aggregate_id).Get(node1,node2,node3);
            if(region_id==0) neighbor=mesh.Adjacent_Triangle(ray.aggregate_id,node1,node2);
            else if(region_id==1) neighbor=mesh.Adjacent_Triangle(ray.aggregate_id,node2,node3);
            else neighbor=mesh.Adjacent_Triangle(ray.aggregate_id,node3,node1);
            if(neighbor==-1){if(triangle.Lazy_Inside_Plane(ray.endpoint)) inside=true;}
            else{
                TRIANGLE_3D<T>& triangle2=(*triangle_list)(neighbor);
                int convex=0;
                if(region_id==0){if(TV::Dot_Product(triangle2.Normal(),triangle.X.z-triangle2.X.x) >= 0) convex=1;}
                else if(region_id==1){if(TV::Dot_Product(triangle2.Normal(),triangle.X.x-triangle2.X.x) >= 0) convex=1;}
                else if(region_id==2){if(TV::Dot_Product(triangle2.Normal(),triangle.X.y-triangle2.X.x) >= 0) convex=1;}
                if(convex){if(triangle.Lazy_Inside_Plane(ray.endpoint) && triangle2.Lazy_Inside_Plane(ray.endpoint)) inside=true;} // inside both - can use location or point
                else{if(triangle.Lazy_Inside_Plane(ray.endpoint) || triangle2.Lazy_Inside_Plane(ray.endpoint)) inside=true;}}} // inside either - can use location or point
        else{if(triangle.Lazy_Inside_Plane(ray.endpoint)) inside=true;}} // region=2 - face - can use location or point
    return inside;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Outside(const TV& location,const T thickness_over_two) const
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(bounding_box && hierarchy && mesh.adjacent_elements && triangle_list);
    if(bounding_box->Outside(location,thickness_over_two)) return true;
    if(hierarchy->box_hierarchy(hierarchy->root).Outside(location,thickness_over_two)) return true;
       
    bool outside=false;
    RAY<TV> ray(location,TV(0,0,1));
    if(!INTERSECTION::Intersects(ray,*this,thickness_over_two)) outside=true; // missed the object, outside
    else if(ray.t_max > 0){ // not in boundary region
        TV point=ray.Point(ray.t_max); // point is inside if and only if location is inside
        T thickness=2*thickness_over_two;
        TRIANGLE_3D<T>& triangle=(*triangle_list)(ray.aggregate_id);
        int region_id,region=triangle.Region(point,region_id,thickness);
        if(region==0){ // vertex
            if(Signed_Solid_Angle_Of_Triangle_Web(point,mesh.elements(ray.aggregate_id)(region_id)) < 0) outside=true;} 
        else if(region==1){ // edge
            int node1,node2,node3,neighbor=-1;mesh.elements(ray.aggregate_id).Get(node1,node2,node3);
            if(region_id==0) neighbor=mesh.Adjacent_Triangle(ray.aggregate_id,node1,node2);
            else if(region_id==1) neighbor=mesh.Adjacent_Triangle(ray.aggregate_id,node2,node3);
            else neighbor=mesh.Adjacent_Triangle(ray.aggregate_id,node3,node1);
            if(neighbor==-1){if(triangle.Lazy_Outside_Plane(location)) outside=true;}
            else{
                TRIANGLE_3D<T>& triangle2=(*triangle_list)(neighbor);
                bool convex=false;
                if(region_id==0){if(TV::Dot_Product(triangle2.Normal(),triangle.X.z-triangle2.X.x)>=0) convex=true;}
                else if(region_id==1){if(TV::Dot_Product(triangle2.Normal(),triangle.X.x-triangle2.X.x)>=0) convex=true;}
                else if(region_id==2){if(TV::Dot_Product(triangle2.Normal(),triangle.X.y-triangle2.X.x)>=0) convex=true;}
                if(convex){if(triangle.Lazy_Outside_Plane(location) || triangle2.Lazy_Outside_Plane(location)) outside=true;} // outside either - can use location or point
                else{if(triangle.Lazy_Outside_Plane(location) && triangle2.Lazy_Outside_Plane(location)) outside=true;}}} // outside both - can use location or point
        else{if(triangle.Lazy_Outside_Plane(location)) outside=true;}} // region=2 - face - can use location or point
    return outside;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Boundary(const TV& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Inside_Any_Triangle
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Inside_Any_Triangle(const TV& location,int& triangle_id,const T thickness_over_two) const
{
    assert(hierarchy);assert(triangle_list);
    ARRAY<int> nearby_triangles;hierarchy->Intersection_List(location,nearby_triangles,thickness_over_two);
    for(int k=0;k<nearby_triangles.m;k++){
        TRIANGLE_3D<T>& triangle=(*triangle_list)(nearby_triangles(k));
        if(triangle.Point_Inside_Triangle(location,thickness_over_two)){triangle_id=nearby_triangles(k);return true;}}
    return false;
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,3> TRIANGULATED_SURFACE<T>::
Surface(const TV& location,const T max_depth,const T thickness_over_2,int* closest_triangle,T* distance) const
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(triangle_list);
    TV point;
    
    if(max_depth){
        PHYSBAM_ASSERT(hierarchy);
        RANGE<TV> box(location-max_depth,location+max_depth);
        ARRAY<int> nearby_triangles;hierarchy->Intersection_List(box,nearby_triangles);
        if(!nearby_triangles.m){ // grab any point assuming far from the interface
            RAY<TV> ray(location,particles.X(mesh.elements(0)(0))-location);if(closest_triangle) *closest_triangle=0;
            if(INTERSECTION::Intersects(ray,*this,thickness_over_2)){if(distance) *distance=ray.t_max;return ray.Point(ray.t_max);}
            else{if(distance) *distance=(particles.X(mesh.elements(0)(0))-location).Magnitude();return particles.X(mesh.elements(0)(0));}} 
        else{
            TV weights;
            {TRIANGLE_3D<T>& triangle=(*triangle_list)(nearby_triangles(0));point=triangle.Closest_Point(location,weights);}
            T distance_temp=(location-point).Magnitude_Squared();if(closest_triangle) *closest_triangle=nearby_triangles(0);
            for(int k=1;k<nearby_triangles.m;k++){
                TRIANGLE_3D<T>& triangle=(*triangle_list)(nearby_triangles(k));
                TV new_point=triangle.Closest_Point(location,weights);
                T new_distance=(location-new_point).Magnitude_Squared();
                if(new_distance < distance_temp){distance_temp=new_distance;point=new_point;if(closest_triangle) *closest_triangle=nearby_triangles(k);}}
            if(distance) *distance=sqrt(distance_temp);return point;}}

    // slow method
    TV weights;
    {TRIANGLE_3D<T>& triangle=(*triangle_list)(0);point=triangle.Closest_Point(location,weights);}
    T distance_temp=(location-point).Magnitude_Squared();if(closest_triangle) *closest_triangle=0;
    for(int k=1;k<mesh.elements.m;k++){
        TRIANGLE_3D<T>& triangle=(*triangle_list)(k);
        TV new_point=triangle.Closest_Point(location,weights);
        T new_distance=(location-new_point).Magnitude_Squared();
        if(new_distance < distance_temp){distance_temp=new_distance;point=new_point;if(closest_triangle) *closest_triangle=k;}}
    if(distance) *distance=sqrt(distance_temp);return point;
}
//#####################################################################
// Function Oriented_Surface
//#####################################################################
template<class T> VECTOR<T,3> TRIANGULATED_SURFACE<T>::
Oriented_Surface(const TV& location,const TV& normal,const T max_depth,const T thickness_over_2,int* closest_triangle,T* distance) const
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(triangle_list);
    TV point;
    
    if(max_depth){
        PHYSBAM_ASSERT(hierarchy);
        RANGE<TV> box(location-max_depth,location+max_depth);
        ARRAY<int> nearby_triangles;hierarchy->Intersection_List(box,nearby_triangles);
        if(!nearby_triangles.m){ // grab any point ignoring normal assuming far from the interface
            RAY<TV> ray(location,particles.X(mesh.elements(0)(0))-location);if(closest_triangle) *closest_triangle=0;
            if(INTERSECTION::Intersects(ray,*this,thickness_over_2)){if(distance) *distance=ray.t_max;return ray.Point(ray.t_max);}
            else{if(distance) *distance=(particles.X(mesh.elements(0)(0))-location).Magnitude();return particles.X(mesh.elements(0)(0));}} 
        else{
            TV weights;
            {TRIANGLE_3D<T>& triangle=(*triangle_list)(nearby_triangles(0));point=triangle.Closest_Point(location,weights);}
            T distance_temp=(location-point).Magnitude_Squared();if(closest_triangle) *closest_triangle=nearby_triangles(0);
            for(int k=1;k<nearby_triangles.m;k++){
                TRIANGLE_3D<T>& triangle=(*triangle_list)(nearby_triangles(k));
                TV new_point=triangle.Closest_Point(location,weights);
                T new_distance=(location-new_point).Magnitude_Squared();
                if(new_distance<distance_temp && TV::Dot_Product(normal,triangle.Raw_Normal())>0){distance_temp=new_distance;point=new_point;if(closest_triangle) *closest_triangle=nearby_triangles(k);}}
            if(distance) *distance=sqrt(distance_temp);return point;}}

    // slow method
    TV weights;
    {TRIANGLE_3D<T>& triangle=(*triangle_list)(0);point=triangle.Closest_Point(location,weights);}
    T distance_temp=(location-point).Magnitude_Squared();if(closest_triangle) *closest_triangle=0;
    for(int k=1;k<mesh.elements.m;k++){
        TRIANGLE_3D<T>& triangle=(*triangle_list)(k);
        TV new_point=triangle.Closest_Point(location,weights);
        T new_distance=(location-new_point).Magnitude_Squared();
        if(new_distance<distance_temp&&TV::Dot_Product(normal,triangle.Raw_Normal())>0){distance_temp=new_distance;point=new_point;if(closest_triangle) *closest_triangle=k;}}
    if(distance) *distance=sqrt(distance_temp);return point;
}
//#####################################################################
// Function Signed_Solid_Angle_Of_Triangle_Web
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Signed_Solid_Angle_Of_Triangle_Web(const TV& location,int web_root_node) const
{
    assert(mesh.incident_elements);
    T signed_solid_angle=0;
    for(int i=0;i<(*mesh.incident_elements)(web_root_node).m;i++){
        int t=(*mesh.incident_elements)(web_root_node)(i);
        int node1,node2,node3;mesh.elements(t).Get(node1,node2,node3);
        TV x1=particles.X(node1),x2=particles.X(node2),x3=particles.X(node3);
        SPHERE<TV> sphere(location,1);
        // extend two of the triangle edges so their farther endpoints are on the unit sphere
        if(web_root_node == node1){RAY<TV> ray1(x1,x2-x1);INTERSECTION::Intersects(ray1,sphere,(T)0);x2=ray1.Point(ray1.t_max);RAY<TV> ray2(x1,x3-x1);INTERSECTION::Intersects(ray2,sphere,(T)0);x3=ray2.Point(ray2.t_max);} 
        else if(web_root_node == node2){RAY<TV> ray1(x2,x1-x2);INTERSECTION::Intersects(ray1,sphere,(T)0);x1=ray1.Point(ray1.t_max);RAY<TV> ray2(x2,x3-x2);INTERSECTION::Intersects(ray2,sphere,(T)0);x3=ray2.Point(ray2.t_max);} 
        else{RAY<TV> ray1(x3,x1-x3);INTERSECTION::Intersects(ray1,sphere,(T)0);x1=ray1.Point(ray1.t_max);RAY<TV> ray2(x3,x2-x3);INTERSECTION::Intersects(ray2,sphere,(T)0);x2=ray2.Point(ray2.t_max);}
        signed_solid_angle+=TRIANGLE_3D<T>(x1,x2,x3).Signed_Solid_Angle(location);} 
    return signed_solid_angle;
}
//#####################################################################
// Function Check_For_Self_Intersection
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Check_For_Self_Intersection(const T thickness_over_2,const bool update_bounding_boxes,ARRAY<VECTOR<int,2> >* intersecting_segment_triangle_pairs)
{  
    bool segment_mesh_defined=mesh.segment_mesh!=0;if(!segment_mesh_defined) mesh.Initialize_Segment_Mesh();
    bool intersection=Segment_Triangle_Intersection(*mesh.segment_mesh,particles.X,thickness_over_2,update_bounding_boxes,intersecting_segment_triangle_pairs);
    if(!segment_mesh_defined){delete mesh.segment_mesh;mesh.segment_mesh=0;}
    return intersection;
}
//#####################################################################
// Function Find_First_Segment_Triangle_Intersection
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Find_First_Segment_Triangle_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,const int max_coarsening_attempts,
    const bool update_bounding_boxes)
{
    if(Segment_Triangle_Intersection(test_segment_mesh,X,thickness_over_2,update_bounding_boxes)){
        LOG::cout<<"SELF INTERSECTIONS !"<<std::endl;
        return true;}
    else{
        for(int loops=0;loops<max_coarsening_attempts;loops++){
            T distance=(1<<loops)*thickness_over_2;
            if(Segment_Triangle_Intersection(test_segment_mesh,X,distance,false)){
                LOG::cout<<"collision at a proximity < "<<distance<<std::endl;
                return true;}
            else LOG::cout<<"ok at a proximity = "<<distance<<std::endl;
    }}
    return false;
}
//#####################################################################
// Function Segment_Triangle_Intersection
//#####################################################################
template<class T> bool TRIANGULATED_SURFACE<T>::
Segment_Triangle_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,const bool update_bounding_boxes,
    ARRAY<VECTOR<int,2> >* intersecting_segment_triangle_pairs)
{  
    bool intersection=false;
    ARRAY<ARRAY<int> > triangles_near_edges(test_segment_mesh.elements.m);Get_Triangles_Near_Edges(triangles_near_edges,test_segment_mesh,X,thickness_over_2,update_bounding_boxes);
    for(int e=0;e<test_segment_mesh.elements.m;e++){
        SEGMENT_3D<T> segment(X(test_segment_mesh.elements(e)(0)),X(test_segment_mesh.elements(e)(1)));
        for(int k=0;k<triangles_near_edges(e).m;k++){int t=triangles_near_edges(e)(k);
            TRIANGLE_3D<T> triangle(X(mesh.elements(t)(0)),X(mesh.elements(t)(1)),X(mesh.elements(t)(2)));
            if(INTERSECTION::Intersects(segment,triangle,thickness_over_2)){intersection=true;
                if(intersecting_segment_triangle_pairs) intersecting_segment_triangle_pairs->Append(VECTOR<int,2>(e,t));else return true;}}}
    return intersection;
}
//#####################################################################
// Function Get_Triangles_Near_Edges
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Get_Triangles_Near_Edges(ARRAY<ARRAY<int> >& triangles_near_edges,const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,
    const bool update_bounding_boxes)
{  
    bool hierarchy_defined=hierarchy!=0;if(!hierarchy_defined) Initialize_Hierarchy(false);
    if(!hierarchy_defined || update_bounding_boxes) hierarchy->Update_Boxes(X);

    for(int k=0;k<test_segment_mesh.elements.m;k++){
        int node1,node2;test_segment_mesh.elements(k).Get(node1,node2);
        RANGE<VECTOR<T,3> > box(X(node1));box.Enlarge_To_Include_Point(X(node2));
        hierarchy->Intersection_List(box,triangles_near_edges(k),thickness_over_2);
        for(int kk=0;kk<triangles_near_edges(k).m;kk++){int t=triangles_near_edges(k)(kk); // remove elements that contain a node of the edge being tested
            for(int i=0;i<3;i++) if(mesh.elements(t)(i) == node1 || mesh.elements(t)(i) == node2){triangles_near_edges(k).Remove_Index_Lazy(kk);kk--;break;}}}
    
    if(!hierarchy_defined){delete hierarchy;hierarchy=0;}
}
//#####################################################################
// Function Centroid_Of_Neighbors
//#####################################################################
template<class T> VECTOR<T,3> TRIANGULATED_SURFACE<T>::
Centroid_Of_Neighbors(const int node) const
{
    assert(mesh.neighbor_nodes);
    TV target;
    int number_of_neighbors=(*mesh.neighbor_nodes)(node).m;
    for(int j=0;j<number_of_neighbors;j++) target+=particles.X((*mesh.neighbor_nodes)(node)(j));
    if(number_of_neighbors != 0) target/=(T)number_of_neighbors;
    else target=particles.X(node); // if no neighbors, return the current node location
    return target;
}
//#####################################################################
// Function Calculate_Signed_Distance
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Calculate_Signed_Distance(const TV& location,T thickness) const
{
    TV closest_point=Surface(location); // this uses the slow (but accurate) method
    T distance=(closest_point-location).Magnitude();
    if(Inside(location,thickness)) distance*=-1;
    return distance;
}
//#####################################################################
// Function Linearly_Subdivide
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Linearly_Subdivide()
{
    TRIANGLE_SUBDIVISION subdivision(mesh);
    TRIANGLE_MESH refined_mesh;
    subdivision.Refine_Mesh(refined_mesh);
    ARRAY<VECTOR<T,3> > X_save(particles.X),V_save(particles.V);
    particles.Add_Elements(refined_mesh.number_nodes-particles.Size());
    if(X_save.Size()) subdivision.Apply_Linear_Subdivision(X_save,particles.X);
    if(V_save.Size()) subdivision.Apply_Linear_Subdivision(V_save,particles.V);
    mesh.Initialize_Mesh(refined_mesh);
}
//#####################################################################
// Function Loop_Subdivide
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Loop_Subdivide()
{
    TRIANGLE_SUBDIVISION subdivision(mesh);
    TRIANGLE_MESH refined_mesh;
    subdivision.Refine_Mesh(refined_mesh);
    ARRAY<VECTOR<T,3> > X_save(particles.X),V_save(particles.V);
    particles.Add_Elements(refined_mesh.number_nodes-particles.Size());
    if(X_save.Size()) subdivision.Apply_Loop_Subdivision(X_save,particles.X);
    if(V_save.Size()) subdivision.Apply_Loop_Subdivision(V_save,particles.V);
    mesh.Initialize_Mesh(refined_mesh);
}
//#####################################################################
// Function Root_Three_Subdivide
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Root_Three_Subdivide()
{
    TRIANGLE_SUBDIVISION subdivision(mesh);
    TRIANGLE_MESH refined_mesh;
    subdivision.Refine_Mesh_Dual(refined_mesh);
    ARRAY<VECTOR<T,3> > X_save(particles.X),V_save(particles.V);
    particles.Add_Elements(refined_mesh.number_nodes-particles.Size());
    if(X_save.Size()) subdivision.Apply_Root_Three_Subdivision(X_save,particles.X);
    if(V_save.Size()) subdivision.Apply_Root_Three_Subdivision(V_save,particles.V);
    mesh.Initialize_Mesh(refined_mesh);
}
//#####################################################################
// Funcion Total_Area
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Total_Area() const
{
    T area=0;
    for(int t=0;t<mesh.elements.m;t++){
        int node0=mesh.elements(t)(0),node1=mesh.elements(t)(1),node2=mesh.elements(t)(2);
        area+=TRIANGLE_3D<T>::Area(particles.X(node0),particles.X(node1),particles.X(node2));}
    return area;
}
//#####################################################################
// Function Minimum_Angle
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Minimum_Angle(int* index) const
{
    TV u,v,w;T min_cosine=1;
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        u=(particles.X(i)-particles.X(j)).Normalized();v=(particles.X(j)-particles.X(k)).Normalized();w=(particles.X(k)-particles.X(i)).Normalized();
        T local_min=min(TV::Dot_Product(u,v),TV::Dot_Product(v,w),TV::Dot_Product(w,u));
        if(local_min < min_cosine){min_cosine=local_min;if(index) *index=t;}}
    return acos(max(-min_cosine,(T)-1));
}
//#####################################################################
// Function Maximum_Angle
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Maximum_Angle(int* index) const
{
    TV u,v,w;T max_cosine=-1;
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        u=(particles.X(i)-particles.X(j)).Normalized();v=(particles.X(j)-particles.X(k)).Normalized();w=(particles.X(k)-particles.X(i)).Normalized();
        T local_max=max(TV::Dot_Product(u,v),TV::Dot_Product(v,w),TV::Dot_Product(w,u));
        if(local_max > max_cosine){max_cosine=local_max;if(index) *index=t;}}
    return acos(min(-max_cosine,(T)1));
}
//#####################################################################
// Function Average_Minimum_Angle
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Average_Minimum_Angle() const
{
    TV u,v,w;T total_min_angle=0;
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        u=(particles.X(i)-particles.X(j)).Normalized();v=(particles.X(j)-particles.X(k)).Normalized();w=(particles.X(k)-particles.X(i)).Normalized();
        total_min_angle+=acos(max(-min(TV::Dot_Product(u,v),TV::Dot_Product(v,w),TV::Dot_Product(w,u)),(T)-1));}
    return total_min_angle/=mesh.elements.m;
}
//#####################################################################
// Function Average_Maximum_Angle
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Average_Maximum_Angle() const
{
    TV u,v,w;T total_max_angle=0;
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        u=(particles.X(i)-particles.X(j)).Normalized();v=(particles.X(j)-particles.X(k)).Normalized();w=(particles.X(k)-particles.X(i)).Normalized();
        total_max_angle+=acos(min(-max(TV::Dot_Product(u,v),TV::Dot_Product(v,w),TV::Dot_Product(w,u)),(T)1));}
    return total_max_angle/=mesh.elements.m;
}
//#####################################################################
// Function Minimum_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Minimum_Edge_Length(int* index) const
{
    T min_edge_length_squared=FLT_MAX;
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        T local_min_squared=min((particles.X(i)-particles.X(j)).Magnitude_Squared(),(particles.X(j)-particles.X(k)).Magnitude_Squared(),(particles.X(k)-particles.X(i)).Magnitude_Squared());
        if(local_min_squared < min_edge_length_squared){min_edge_length_squared=local_min_squared;if(index) *index=t;}}
    return sqrt(min_edge_length_squared);
}
//#####################################################################
// Function Maximum_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Maximum_Edge_Length(int* index) const
{
    T max_edge_length_squared=0;
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        T local_max_squared=max((particles.X(i)-particles.X(j)).Magnitude_Squared(),(particles.X(j)-particles.X(k)).Magnitude_Squared(),(particles.X(k)-particles.X(i)).Magnitude_Squared());
        if(local_max_squared > max_edge_length_squared){max_edge_length_squared=local_max_squared;if(index) *index=t;}}
    return sqrt(max_edge_length_squared);
}
//#####################################################################
// Function Average_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Average_Edge_Length() const
{
    bool segment_mesh_defined=mesh.segment_mesh!=0;if(!segment_mesh_defined) mesh.Initialize_Segment_Mesh();
    T average_edge_length=0;
    for(int s=0;s<mesh.segment_mesh->elements.m;s++){
        int i,j;mesh.segment_mesh->elements(s).Get(i,j);
        average_edge_length+=(particles.X(i)-particles.X(j)).Magnitude();}
    if(mesh.segment_mesh->elements.m) average_edge_length/=mesh.segment_mesh->elements.m;
    if(!segment_mesh_defined){delete mesh.segment_mesh;mesh.segment_mesh=0;}
    return average_edge_length;
}
//#####################################################################
// Function Maximum_Aspect_Ratio
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Maximum_Aspect_Ratio(int* index) const
{
    T max_aspect_ratio=0;
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        T local_max=TRIANGLE_3D<T>::Aspect_Ratio(particles.X(i),particles.X(j),particles.X(k));
        if(local_max > max_aspect_ratio){max_aspect_ratio=local_max;if(index) *index=t;}}
    return max_aspect_ratio;
}
//#####################################################################
// Function Average_Aspect_Ratio
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Average_Aspect_Ratio() const
{
    T total_aspect_ratio=0;
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        total_aspect_ratio+=TRIANGLE_3D<T>::Aspect_Ratio(particles.X(i),particles.X(j),particles.X(k));}
    return total_aspect_ratio/mesh.elements.m;
}
//#####################################################################
// Funcion Minimum_Area
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Minimum_Area(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=0;t<mesh.elements.m;t++){
        int node0,node1,node2;mesh.elements(t).Get(node0,node1,node2);
        T temp=TRIANGLE_3D<T>::Area(particles.X(node0),particles.X(node1),particles.X(node2));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Altitude
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Minimum_Altitude(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=0;t<mesh.elements.m;t++){
        int node0=mesh.elements(t)(0),node1=mesh.elements(t)(1),node2=mesh.elements(t)(2);
        T temp=TRIANGLE_3D<T>::Minimum_Altitude(particles.X(node0),particles.X(node1),particles.X(node2));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Function Maximum_Magnitude_Phi
//#####################################################################
template<class T> T TRIANGULATED_SURFACE<T>::
Maximum_Magnitude_Phi(const IMPLICIT_OBJECT<TV>& implicit_surface,int* index)
{
    T phi=0,max_phi=0;int k=0;
    for(int i=0;i<particles.Size();i++){phi=abs(implicit_surface(particles.X(i)));if(phi > max_phi){max_phi=phi;k=i;}}
    if(index)*index=k;return max_phi;
}
//#####################################################################
// Function Make_Orientations_Consistent_With_Implicit_Surface
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Make_Orientations_Consistent_With_Implicit_Surface(const IMPLICIT_OBJECT<TV>& implicit_surface)
{
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        TV centroid=(T)one_third*(particles.X(i)+particles.X(j)+particles.X(k));
        TV normal=PLANE<T>::Normal(particles.X(i),particles.X(j),particles.X(k));
        if(TV::Dot_Product(normal,implicit_surface.Normal(centroid))<0)mesh.elements(t).Set(i,k,j);}
}
//#####################################################################
// Function Close_Surface
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Close_Surface(const bool merge_coincident_vertices,const T merge_coincident_vertices_threshold,const bool fill_holes,const bool verbose)
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(merge_coincident_vertices_threshold<1);
    bool incident_elements_defined=(mesh.incident_elements!=0);if(!incident_elements_defined) mesh.Initialize_Incident_Elements();
    bool boundary_mesh_defined=(mesh.boundary_mesh!=0);if(!boundary_mesh_defined) mesh.Initialize_Boundary_Mesh();
    bool node_on_boundary_defined=(mesh.node_on_boundary!=0);if(!node_on_boundary_defined) mesh.Initialize_Node_On_Boundary();

    // for each node on boundary, get minimum length of boundary segments adjacent to it - also keep track of longest boundary segment
    ARRAY<T> minimum_incident_boundary_segment_length(mesh.boundary_mesh->number_nodes,false);
    minimum_incident_boundary_segment_length.Fill(FLT_MAX);
    T maximum_boundary_segment_length=0;
    for(int i=0;i<mesh.boundary_mesh->elements.m;i++){
        int node1,node2;mesh.boundary_mesh->elements(i).Get(node1,node2);
        T length=(particles.X(node1)-particles.X(node2)).Magnitude();
        minimum_incident_boundary_segment_length(node1)=min(minimum_incident_boundary_segment_length(node1),length);
        minimum_incident_boundary_segment_length(node2)=min(minimum_incident_boundary_segment_length(node2),length);
        maximum_boundary_segment_length=max(maximum_boundary_segment_length,length);}

    if(merge_coincident_vertices && maximum_boundary_segment_length>0){
        if(verbose) LOG::cout<<"MERGING VERTICES (threshold="<<merge_coincident_vertices_threshold<<")"<<std::endl;
        int number_merged=0;
        PARTICLE_3D_SPATIAL_PARTITION<T> particle_spatial_partition(particles,2*maximum_boundary_segment_length);
        particle_spatial_partition.Reinitialize();particle_spatial_partition.Reset_Pair_Finder();
        int index1;ARRAY<int> nearby_particle_indices;nearby_particle_indices.Preallocate(100);
        while(particle_spatial_partition.Get_Next_Particles_Potentially_Within_Interaction_Radius(index1,nearby_particle_indices)) if((*mesh.node_on_boundary)(index1)){
            TV position1=particles.X(index1);
            int closest_index2=0;T closest_distance_squared=FLT_MAX;
            for(int k=0;k<nearby_particle_indices.m;k++){
                int index2=nearby_particle_indices(k);if(!(*mesh.node_on_boundary)(index2)) continue;
                T distance_squared=(particles.X(index2)-position1).Magnitude_Squared();
                if(distance_squared < closest_distance_squared){closest_index2=index2;closest_distance_squared=distance_squared;}}
            if(closest_index2){ // merge index1 to closest_index2
                T real_threshold=merge_coincident_vertices_threshold*min(minimum_incident_boundary_segment_length(index1),minimum_incident_boundary_segment_length(closest_index2));
                if(closest_distance_squared<=sqr(real_threshold)){
                    number_merged++;
                    if(verbose) LOG::cout<<index1<<"->"<<closest_index2<<" "<<std::flush;
                    for(int j=0;j<(*mesh.incident_elements)(index1).m;j++){
                        int t=(*mesh.incident_elements)(index1)(j);
                        if(mesh.elements(t)(0)==index1)mesh.elements(t)(0)=closest_index2;
                        else if(mesh.elements(t)(1)==index1)mesh.elements(t)(1)=closest_index2;
                        else if(mesh.elements(t)(2)==index1)mesh.elements(t)(2)=closest_index2;}
                    (*mesh.incident_elements)(closest_index2).Append_Elements((*mesh.incident_elements)(index1));}}} // dynamically update the incident triangles list
        if(verbose) LOG::cout<<std::endl<<number_merged<<" vertices merged"<<std::endl<<std::endl;
        Refresh_Auxiliary_Structures();}

    if(fill_holes){
        if(verbose) LOG::cout<<"HOLE FILLING"<<std::endl;
        bool connected_segments_defined=(mesh.boundary_mesh->connected_segments!=0);if(!connected_segments_defined) mesh.boundary_mesh->Initialize_Connected_Segments();
        ARRAY<ARRAY<VECTOR<int,2> > >& connected_segments=*mesh.boundary_mesh->connected_segments;
        for(int i=0;i<connected_segments.m;i++){
            TV centroid;
            for(int j=0;j<connected_segments(i).m;j++){int node1,node2;connected_segments(i)(j).Get(node1,node2);centroid+=particles.X(node1)+particles.X(node2);}
            centroid/=(T)(2*connected_segments(i).m); // assuming we count each node exactly twice, this gives us the average
            int new_particle_index=particles.Add_Element();mesh.number_nodes++;particles.X(new_particle_index)=centroid;
            if(particles.store_velocity){
                TV velocity;
                for(int j=0;j<connected_segments(i).m;j++){int node1,node2;connected_segments(i)(j).Get(node1,node2);velocity+=particles.V(node1)+particles.V(node2);}
                particles.V(new_particle_index)=velocity/(T)(2*connected_segments(i).m);}
            if(verbose){LOG::cout<<"Adding particle "<<new_particle_index<<": "<<particles.X(new_particle_index)<<std::endl;LOG::cout<<"Adding triangles: "<<std::flush;}
            for(int j=0;j<connected_segments(i).m;j++){ // assumes segment orientation is consistent with triangle orientation!
                int node1,node2;connected_segments(i)(j).Get(node1,node2);
                if(verbose) LOG::cout<<"("<<node2<<","<<node1<<","<<new_particle_index<<") "<<std::flush;
                mesh.elements.Append(VECTOR<int,3>(node2,node1,new_particle_index));}
            if(verbose) LOG::cout<<std::endl<<std::endl;
        }
        if(!connected_segments_defined){delete mesh.boundary_mesh->connected_segments;mesh.boundary_mesh->connected_segments=0;}}

    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
    if(!boundary_mesh_defined){delete mesh.boundary_mesh;mesh.boundary_mesh=0;}
    if(!node_on_boundary_defined){delete mesh.node_on_boundary;mesh.node_on_boundary=0;}

    Discard_Valence_Zero_Particles_And_Renumber();  // refreshes auxiliary structures too 
}
//#####################################################################
// Function Remove_Degenerate_Triangles
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Remove_Degenerate_Triangles(const T area_threshold)
{
    T area_threshold_squared=sqr(area_threshold);
    for(int t=mesh.elements.m-1;t>=0;t--){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        if(TRIANGLE_3D<T>::Area_Squared(particles.X(i),particles.X(j),particles.X(k))<area_threshold_squared) mesh.elements.Remove_Index_Lazy(t);}
    Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Create_Compact_Copy
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* TRIANGULATED_SURFACE<T>::
Create_Compact_Copy() const
{
    ARRAY<int> new_to_old;
    HASHTABLE<int,int> old_to_new;
    TRIANGLE_MESH* triangle_mesh=new TRIANGLE_MESH;
    triangle_mesh->elements.Resize(mesh.elements.m);
    for(int i=0;i<mesh.elements.m;i++){const VECTOR<int,3>& element=mesh.elements(i);
        for(int j=0;j<3;j++){
            int& a=old_to_new.Get_Or_Insert(element(j));
            if(a<0) a=new_to_old.Append(element(j));
            triangle_mesh->elements(i)(j)=a;}}

    GEOMETRY_PARTICLES<TV>* deformable_geometry_particle=new GEOMETRY_PARTICLES<TV>;
    deformable_geometry_particle->Add_Elements(new_to_old.Size());
    deformable_geometry_particle->X.Prefix(deformable_geometry_particle->Size())=particles.X.Subset(new_to_old);
    TRIANGULATED_SURFACE* triangulated_surface=new TRIANGULATED_SURFACE(*triangle_mesh,*deformable_geometry_particle);
    triangulated_surface->Update_Number_Nodes();
    return triangulated_surface;
}
//#####################################################################
// Function Print_Statistics
//#####################################################################
template<class T> void TRIANGULATED_SURFACE<T>::
Print_Statistics(std::ostream& output,const T thickness_over_2)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    int index;Update_Bounding_Box();
    if(!mesh.incident_elements) mesh.Initialize_Incident_Elements();

    output<<"triangles = "<<mesh.elements.m<<std::endl;
    output<<"particles = "<<particles.Size()<<std::endl;
    {int particles_touched=0;for(int p=0;p<particles.Size();p++) if((*mesh.incident_elements)(p).m) particles_touched++;
    output<<"particles touched = "<<particles_touched<<std::endl;}
    output<<"bounding box = "<<*bounding_box<<std::endl;
    if(particles.store_velocity){
        int index=particles.V.Arg_Maximum_Magnitude();
        output<<"max_speed = "<<particles.V(index).Magnitude()<<" ("<<index<<")"<<std::endl;}
    output<<"total area = "<<Total_Area()<<std::endl;
    output<<"min area = "<<Minimum_Area(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"min altitude = "<<Minimum_Altitude(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"min edge length = "<<Minimum_Edge_Length(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"max edge length = "<<Maximum_Edge_Length(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"min angle = "<<Minimum_Angle(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"max angle = "<<Maximum_Angle(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"ave min angle = "<<Average_Minimum_Angle()<<std::endl;
    output<<"ave max angle = "<<Average_Maximum_Angle()<<std::endl;
    output<<"max aspect ratio = "<<Maximum_Aspect_Ratio(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"ave aspect ratio = "<<Average_Aspect_Ratio()<<std::endl;
    if(Check_For_Self_Intersection(thickness_over_2)) output<<"found self intersections"<<std::endl;else output<<"no self intersections"<<std::endl;
}
//#####################################################################
template class TRIANGULATED_SURFACE<float>;
template class TRIANGULATED_SURFACE<double>;
