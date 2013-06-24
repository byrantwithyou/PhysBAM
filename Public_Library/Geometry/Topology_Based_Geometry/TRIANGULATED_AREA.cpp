//#####################################################################
// Copyright 2003-2008, Christopher Allocco, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Neil Molino, Duc Nguyen, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Log/LOG.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> TRIANGULATED_AREA<T>::
TRIANGULATED_AREA()
    :MESH_OBJECT<TV,TRIANGLE_MESH>(*new TRIANGLE_MESH,*new GEOMETRY_PARTICLES<TV>),segmented_curve(0),hierarchy(0),triangle_area_fractions(0),triangle_areas(0),
    nodal_areas(0)
{
    this->need_destroy_mesh=this->need_destroy_particles=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> TRIANGULATED_AREA<T>::
TRIANGULATED_AREA(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,TRIANGLE_MESH>(triangle_mesh_input,particles_input),segmented_curve(0),hierarchy(0),triangle_area_fractions(0),triangle_areas(0),
    nodal_areas(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TRIANGULATED_AREA<T>::
~TRIANGULATED_AREA()
{
    delete segmented_curve;delete hierarchy;delete triangle_area_fractions;delete triangle_areas;delete nodal_areas;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Clean_Memory()
{
    MESH_OBJECT<TV,TRIANGLE_MESH>::Clean_Memory();
    delete segmented_curve;segmented_curve=0;delete hierarchy;hierarchy=0;
    delete triangle_area_fractions;triangle_area_fractions=0;delete triangle_areas;triangle_areas=0;delete nodal_areas;nodal_areas=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures_Helper
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Refresh_Auxiliary_Structures_Helper()
{
    if(hierarchy) Initialize_Hierarchy();
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Hierarchy(const bool update_boxes)
{
    delete hierarchy;hierarchy=new TRIANGLE_HIERARCHY_2D<T>(mesh,particles,update_boxes);
}
//#####################################################################
// Funcion Initialize_Square_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Square_Mesh_And_Particles(const GRID<TV>& grid,const bool reverse_triangles)
{
    int m=grid.counts.x,n=grid.counts.y,particle=0;
    particles.Delete_All_Elements();mesh.Initialize_Square_Mesh(m,n,reverse_triangles);particles.Add_Elements(m*n);
    for(int j=0;j<n;j++) for(int i=0;i<m;i++) particles.X(particle++)=grid.X(TV_INT(i,j));
}
//#####################################################################
// Funcion Initialize_Circle_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Circle_Mesh_And_Particles(const T outer_radius,const T inner_radius,const int num_radial,const int num_tangential)
{
    int particle=0;particles.Delete_All_Elements();
    mesh.Initialize_Circle_Mesh(num_radial,num_tangential);particles.Add_Elements(num_radial*num_tangential);
    for(int j=0;j<num_radial;j++) for(int i=0;i<num_tangential;i++){
        T r=T(j)/T(num_tangential-1)*(outer_radius-inner_radius)+inner_radius,theta=T(i)/T(num_tangential)*(T)2*(T)pi; 
        particles.X(particle++)=VECTOR<T,2>(r*cos(theta),r*sin(theta));}
}
//#####################################################################
// Funcion Initialize_Herringbone_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Herring_Bone_Mesh_And_Particles(const GRID<TV>& grid)
{
    int m=grid.counts.x,n=grid.counts.y,particle=0;
    particles.Delete_All_Elements();mesh.Initialize_Herring_Bone_Mesh(m,n);particles.Add_Elements(m*n);
    for(int j=0;j<n;j++) for(int i=0;i<m;i++) particles.X(particle++)=grid.X(TV_INT(i,j));
}
//#####################################################################
// Funcion Initialize_Equilateral_Mesh_And_Particles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Equilateral_Mesh_And_Particles(const GRID<TV>& grid)
{
    int m=grid.counts.x,n=grid.counts.y,particle=0;
    particles.Delete_All_Elements();mesh.Initialize_Equilateral_Mesh(m,n);particles.Add_Elements(m*n);
    for(int j=0;j<n;j++)
        for(int i=0;i<m;i++){
            TV X=grid.X(TV_INT(i,j));
            if(j%2) X.x+=(T).5*grid.dX.x;
            particles.X(particle++)=X;}
}
//#####################################################################
// Function Inside
// note: return the first triangle that it is inside of (including boundary), otherwise returns 0
//#####################################################################
template<class T> int TRIANGULATED_AREA<T>::
Inside(const TV& location,const T thickness_over_two) const
{
    PHYSBAM_ASSERT(bounding_box && hierarchy);
    if(bounding_box->Outside(location,thickness_over_two)) return -1;
    if(hierarchy->box_hierarchy(hierarchy->root).Outside(location,thickness_over_two)) return -1;
    ARRAY<int> triangles_to_check;hierarchy->Intersection_List(location,triangles_to_check,thickness_over_two);
    for(int l=0;l<triangles_to_check.m;l++){
        int t=triangles_to_check(l);int i,j,k;mesh.elements(t).Get(i,j,k);
        if(!TRIANGLE_2D<T>::Outside(location,particles.X(i),particles.X(j),particles.X(k),thickness_over_two)) return t;}
    return -1;
}
//#####################################################################
// Function Inside_Any_Simplex
//#####################################################################
template<class T> bool TRIANGULATED_AREA<T>::
Inside_Any_Simplex(const VECTOR<T,2>& location,int& triangle_id,const T thickness_over_two) const
{
    assert(hierarchy);
    ARRAY<int> nearby_triangles;hierarchy->Intersection_List(location,nearby_triangles,thickness_over_two);
    for(int k=0;k<nearby_triangles.m;k++){
        const TRIANGLE_2D<T>& triangle=Get_Element(nearby_triangles(k));
        if(triangle.Inside(location,thickness_over_two)){triangle_id=nearby_triangles(k);return true;}}
    return false;
}
//#####################################################################
// Function Check_Signed_Area_And_Make_Consistent
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Check_Signed_Area_And_Make_Consistent(const int triangle,const bool verbose)
{
    VECTOR<int,3>& nodes=mesh.elements(triangle);
    if(TRIANGLE_2D<T>::Signed_Area(particles.X(nodes[0]),particles.X(nodes[1]),particles.X(nodes[2])) < 0){
        if(verbose) LOG::cout<<"triangle number "<<triangle<<" is oriented improperly."<<std::endl;
        exchange(nodes[1],nodes[2]);}
}
//#####################################################################
// Function Check_Signed_Areas_And_Make_Consistent
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Check_Signed_Areas_And_Make_Consistent(const bool verbose)
{
    for(int t=0;t<mesh.elements.m;t++) Check_Signed_Area_And_Make_Consistent(t,verbose);
}
//#####################################################################
// Function Initialize_Segmented_Curve
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Segmented_Curve()
{
    delete segmented_curve;if(!mesh.boundary_mesh) mesh.Initialize_Boundary_Mesh();
    segmented_curve=new SEGMENTED_CURVE_2D<T>(*mesh.boundary_mesh,particles);
}
//#####################################################################
// Function Centroid
//#####################################################################
template<class T> VECTOR<T,2> TRIANGULATED_AREA<T>::
Centroid(const int triangle) const
{
    int i,j,k;mesh.elements(triangle).Get(i,j,k);
    return (T)one_third*(particles.X(i)+particles.X(j)+particles.X(k));
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Rescale(const T scaling_factor)
{
    Rescale(scaling_factor,scaling_factor);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Rescale(const T scaling_x,const T scaling_y)
{
    for(int k=0;k<particles.Size();k++) particles.X(k)*=TV(scaling_x,scaling_y);
    Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Area(const int triangle) const
{
    int i,j,k;mesh.elements(triangle).Get(i,j,k);
    return TRIANGLE_2D<T>::Area(particles.X(i),particles.X(j),particles.X(k));
}
//#####################################################################
// Function Signed_Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Signed_Area(const int triangle) const
{
    int i,j,k;mesh.elements(triangle).Get(i,j,k);
    return TRIANGLE_2D<T>::Signed_Area(particles.X(i),particles.X(j),particles.X(k));
}
//#####################################################################
// Funcion Minimum_Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Minimum_Area(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=0;t<mesh.elements.m;t++){
        int node1=mesh.elements(t)(0),node2=mesh.elements(t)(1),node3=mesh.elements(t)(2);
        T temp=TRIANGLE_2D<T>::Area(particles.X(node1),particles.X(node2),particles.X(node3));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Signed_Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Minimum_Signed_Area(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=0;t<mesh.elements.m;t++){
        int node1=mesh.elements(t)(0),node2=mesh.elements(t)(1),node3=mesh.elements(t)(2);
        T temp=TRIANGLE_2D<T>::Signed_Area(particles.X(node1),particles.X(node2),particles.X(node3));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Funcion Total_Area
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Total_Area() const
{
    T area=0;
    for(int t=0;t<mesh.elements.m;t++){
        int node1=mesh.elements(t)(0),node2=mesh.elements(t)(1),node3=mesh.elements(t)(2);
        area+=TRIANGLE_2D<T>::Area(particles.X(node1),particles.X(node2),particles.X(node3));}
    return area;
}
//#####################################################################
// Funcion Minimum_Altitude
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Minimum_Altitude(int* index) const
{
    int k=0;T minimum=FLT_MAX;
    for(int t=0;t<mesh.elements.m;t++){
        int node1=mesh.elements(t)(0),node2=mesh.elements(t)(1),node3=mesh.elements(t)(2);
        T temp=TRIANGLE_2D<T>::Minimum_Altitude(particles.X(node1),particles.X(node2),particles.X(node3));
        if(temp < minimum){minimum=temp;k=t;}}
    if(index) *index=k;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Minimum_Edge_Length(int* index) const
{
    int t_save=0;T minimum_squared=FLT_MAX;
    for(int t=0;t<mesh.elements.m;t++){int i,j,k;mesh.elements(t).Get(i,j,k);
        T temp=TRIANGLE_2D<T>::Minimum_Edge_Length_Squared(particles.X(i),particles.X(j),particles.X(k));if(temp < minimum_squared){minimum_squared=temp;t_save=t;}}
    if(index) *index=t_save;
    return sqrt(minimum_squared);
}
//#####################################################################
// Funcion Maximum_Edge_Length
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Maximum_Edge_Length(int* index) const
{
    int t_save=0;T maximum_squared=FLT_MAX;
    for(int t=0;t<mesh.elements.m;t++){int i,j,k;mesh.elements(t).Get(i,j,k);
        T temp=TRIANGLE_2D<T>::Maximum_Edge_Length_Squared(particles.X(i),particles.X(j),particles.X(k));if(temp < maximum_squared){maximum_squared=temp;t_save=t;}}
    if(index) *index=t_save;
    return sqrt(maximum_squared);
}
//#####################################################################
// Function Inverted_Triangles
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Inverted_Triangles(ARRAY<int>& inverted_triangles) const
{
    inverted_triangles.Resize(0);
    for(int t=0;t<mesh.elements.m;t++){
        int node1=mesh.elements(t)(0),node2=mesh.elements(t)(1),node3=mesh.elements(t)(2);
        if(TRIANGLE_2D<T>::Signed_Area(particles.X(node1),particles.X(node2),particles.X(node3)) < 0) inverted_triangles.Append(t);}
}
//#####################################################################
// Function Volume_Incident_On_A_Particle
//#####################################################################
template<class T> T TRIANGULATED_AREA<T>::
Area_Incident_On_A_Particle(const int particle_index)
{
    int incident_elements_defined=mesh.incident_elements!=0;if(!incident_elements_defined) mesh.Initialize_Incident_Elements();
    T total_incident_area=0;
    for(int t=0;t<(*mesh.incident_elements)(particle_index).m;t++){
        int i,j,k;mesh.elements((*mesh.incident_elements)(particle_index)(t)).Get(i,j,k);
        total_incident_area+=TRIANGLE_2D<T>::Area(particles.X(i),particles.X(j),particles.X(k));}
    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
    return total_incident_area;
}
//#####################################################################
// Function Split_Node
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements
template<class T> int TRIANGULATED_AREA<T>::
Split_Node(const int particle_index,const TV& normal)
{
    int incident_elements_defined=mesh.incident_elements!=0;if(!incident_elements_defined) mesh.Initialize_Incident_Elements();
    VECTOR<T,2> x0=particles.X(particle_index);ARRAY<int> tris_incident_on_old_particle,tris_incident_on_new_particle;
    int t;for(t=0;t<(*mesh.incident_elements)(particle_index).m;t++){
        int this_incident_tri=(*mesh.incident_elements)(particle_index)(t);int i,j,k;mesh.elements(this_incident_tri).Get(i,j,k);
        VECTOR<T,2> x1=particles.X(i),x2=particles.X(j),x3=particles.X(k),centroid=(T)one_third*(x1+x2+x3);
        if((centroid.x-x0.x)*normal.x+(centroid.y-x0.y)*normal.y < 0) tris_incident_on_new_particle.Append(this_incident_tri); 
        else tris_incident_on_old_particle.Append(this_incident_tri);}
    int new_particle=0;
    if(tris_incident_on_old_particle.m != 0 && tris_incident_on_new_particle.m != 0){ 
        // new particle - assumes we're storing position, and velocity - user must fix mass outside this function call
        new_particle=particles.Add_Element();mesh.number_nodes=particles.Size();
        particles.X(new_particle)=particles.X(particle_index);particles.V(new_particle)=particles.V(particle_index);
        for(t=0;t<(*mesh.incident_elements)(particle_index).m;t++){
            int this_incident_tri=(*mesh.incident_elements)(particle_index)(t);int i,j,k;mesh.elements(this_incident_tri).Get(i,j,k);
            VECTOR<T,2> x1=particles.X(i),x2=particles.X(j),x3=particles.X(k),centroid=(T)one_third*(x1+x2+x3);
            if((centroid.x-x0.x)*normal.x+(centroid.y-x0.y)*normal.y < 0){ // relabel with duplicate node
                if(i == particle_index) i=new_particle;if(j == particle_index) j=new_particle;if(k == particle_index) k=new_particle;
                mesh.elements(this_incident_tri).Set(i,j,k);}}        
        if(incident_elements_defined){ //repair incident triangles if necessary
            (*mesh.incident_elements)(particle_index).Clean_Memory();
            (*mesh.incident_elements)(particle_index).Append_Elements(tris_incident_on_old_particle);
            (*mesh.incident_elements).Append(tris_incident_on_new_particle);}}
    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
    return new_particle;
}
//#####################################################################
// Function Discard_Triangles_Outside_Implicit_Curve
//#####################################################################
// uses Whitney-like criterion to discard only those triangles that are for sure outside levelset (assuming accurate signed distance)
template<class T> void TRIANGULATED_AREA<T>::
Discard_Triangles_Outside_Implicit_Curve(IMPLICIT_OBJECT<TV>& implicit_curve)
{
    VECTOR<T,2> xi,xj,xk;int t=1;
    while(t <= mesh.elements.m){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        xi=particles.X(i);xj=particles.X(j);xk=particles.X(k);
        T max_length=sqrt(max((xi-xj).Magnitude_Squared(),(xj-xk).Magnitude_Squared(),(xk-xi).Magnitude_Squared()));
        T min_phi=min(implicit_curve(xi),implicit_curve(xj),implicit_curve(xk));
        if(min_phi > max_length) mesh.elements.Remove_Index_Lazy(t);else t++;}
}
//#####################################################################
// Function Initialize_Triangle_Area_Fractions_From_Voronoi_Regions
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Initialize_Triangle_Area_Fractions_From_Voronoi_Regions()
{
    if(!triangle_area_fractions) triangle_area_fractions=new ARRAY<VECTOR<T,2> >;
    triangle_area_fractions->Resize(mesh.elements.m);
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        VECTOR<T,3> fractions=VECTOR<T,3>(1,1,1)-TRIANGLE_2D<T>::Circumcenter_Barycentric_Coordinates(particles.X(i),particles.X(j),particles.X(k));
        fractions.x=max(T(0),fractions.x);fractions.y=max(T(0),fractions.y);fractions.z=max(T(0),fractions.z);
        fractions/=(fractions.x+fractions.y+fractions.z);
        (*triangle_area_fractions)(t).Set(fractions.x,fractions.y);}
}
//#####################################################################
// Function Compute_Triangle_Areas
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Compute_Triangle_Areas()
{
    if(!triangle_areas) triangle_areas=new ARRAY<T>;
    triangle_areas->Resize(mesh.elements.m);
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        (*triangle_areas)(t)=TRIANGLE_2D<T>::Signed_Area(particles.X(i),particles.X(j),particles.X(k));}
}
//#####################################################################
// Function Compute_Nodal_Areas
//#####################################################################
template<class T> void TRIANGULATED_AREA<T>::
Compute_Nodal_Areas(bool save_triangle_areas)
{
    if(!nodal_areas) nodal_areas=new ARRAY<T>;
    nodal_areas->Resize(particles.Size());
    if(save_triangle_areas){
        if(!triangle_areas) triangle_areas=new ARRAY<T>;
        triangle_areas->Resize(mesh.elements.m);}
    nodal_areas->Fill(0);
    for(int t=0;t<mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        T area=TRIANGLE_2D<T>::Signed_Area(particles.X(i),particles.X(j),particles.X(k));
        if(save_triangle_areas) (*triangle_areas)(t)=area;
        if(!triangle_area_fractions){area*=(T)one_third;(*nodal_areas)(i)+=area;(*nodal_areas)(j)+=area;(*nodal_areas)(k)+=area;}
        else{T wi,wj;(*triangle_area_fractions)(t).Get(wi,wj);(*nodal_areas)(i)+=wi*area;(*nodal_areas)(j)+=wj*area;(*nodal_areas)(k)+=(1-wi-wj)*area;}}
}
//#####################################################################
// Function Triangle_In_Direction_Uninverted
//#####################################################################
template<class T> int TRIANGULATED_AREA<T>::
Triangle_In_Direction_Uninverted(const int node,const TV& direction) const
{
    assert(mesh.topologically_sorted_neighbor_nodes && mesh.topologically_sorted_incident_elements);
    ARRAY<int> &neighbor_nodes=(*mesh.topologically_sorted_neighbor_nodes)(node),&incident_elements=(*mesh.topologically_sorted_incident_elements)(node);
    TV xnode=particles.X(node);int t=1;
    while(t<neighbor_nodes.m){ // find right edge in right half-space of direction
        if(TV::Cross_Product(particles.X(neighbor_nodes(t))-xnode,direction).x<0){t++;continue;}
        while(t<neighbor_nodes.m){ // find left edge in left half-space of direction
            if(TV::Cross_Product(direction,particles.X(neighbor_nodes(t+1))-xnode).x<0){t++;continue;}
            return t;}
        return t<incident_elements.m?incident_elements(t)-1:-1;}
    return t<incident_elements.m?incident_elements(t)-1:-1;
}
//#####################################################################
// Function Triangle_Walk_Uninverted
//#####################################################################
template<class T> int TRIANGULATED_AREA<T>::
Triangle_Walk_Uninverted(const int start_node,const TV& dX) const
{
    int t=Triangle_In_Direction_Uninverted(start_node,dX);if(t<0) return -1;
    TV start=particles.X(start_node),goal=start+dX;
    int e1,e2;mesh.Other_Two_Nodes(start_node,t,e1,e2);
    if(TV::Cross_Product(particles.X(e2)-particles.X(e1),goal-particles.X(e1)).x>=0) return t;
    for(;;){
        t=mesh.Adjacent_Triangle(t,e1,e2);if(t<0) return -1;
        int e3=mesh.Other_Node(e1,e2,t);
        TV w=TRIANGLE_2D<T>::First_Two_Barycentric_Coordinates(goal,particles.X(e2),particles.X(e1),particles.X(e3));
        if(w.x>=0){if(w.y>=0) return t;else e1=e3;}
        else if(w.y>=0 || TV::Cross_Product(dX,particles.X(e3)-start).x>0) e2=e3;
        else e1=e3;}
}
//#####################################################################
// Function Triangle_Walk_Uninverted
//#####################################################################
template<class T> bool TRIANGULATED_AREA<T>::
Fix_Pair_For_Delaunay(const int triangle1,const int triangle2)
{
    VECTOR<int,3> &t1=mesh.elements(triangle1),&t2=mesh.elements(triangle2);
    int a1i=0;for(int j=0;j<3;j++) if(!t2.Contains(t1(j))) a1i=j;
    int a2i=0;for(int j=0;j<3;j++) if(!t1.Contains(t2(j))) a2i=j;
    int a1=t1(a1i),a2=t2(a2i),b1=t1(a1i>2?a1i-2:a1i+1),b2=t1(a1i>1?a1i-1:a1i+2);

    assert(TRIANGLE_2D<T>(particles.X.Subset(t1)).Check_Orientation() && TRIANGLE_2D<T>(particles.X.Subset(t2)).Check_Orientation());
    assert(b2==t2(a2i>2?a2i-2:a2i+1) && b1==t2(a2i>1?a2i-1:a2i+2));

    if(TRIANGLE_2D<T>::Check_Delaunay_Criterion(particles.X(a1),particles.X(b1),particles.X(b2),particles.X(a2))) return false;

    t1=VECTOR<int,3>(a1,b1,a2);t2=VECTOR<int,3>(a2,b2,a1);
    assert(TRIANGLE_2D<T>(particles.X.Subset(t1)).Check_Orientation() && TRIANGLE_2D<T>(particles.X.Subset(t2)).Check_Orientation());
    return true;
}
//#####################################################################
namespace PhysBAM{
template class TRIANGULATED_AREA<float>;
template class TRIANGULATED_AREA<double>;
}
