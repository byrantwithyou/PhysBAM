//#####################################################################
// Copyright 2003-2009, Ron Fedkiw, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Duc Nguyen, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENTED_CURVE
//##################################################################### 
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SEGMENTED_CURVE<TV>::
SEGMENTED_CURVE()
    :MESH_OBJECT<TV,SEGMENT_MESH>(*new SEGMENT_MESH,*new GEOMETRY_PARTICLES<TV>),hierarchy(0),segment_list(0),point_simplices_1d(0)
{
    this->need_destroy_mesh=this->need_destroy_particles=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SEGMENTED_CURVE<TV>::
SEGMENTED_CURVE(SEGMENT_MESH& segment_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,SEGMENT_MESH>(segment_mesh_input,particles_input),hierarchy(0),segment_list(0),point_simplices_1d(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SEGMENTED_CURVE<TV>::
~SEGMENTED_CURVE()
{
    delete hierarchy;delete segment_list;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::
Clean_Memory()
{
    MESH_OBJECT<TV,SEGMENT_MESH>::Clean_Memory();
    delete hierarchy;hierarchy=0;delete segment_list;segment_list=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures_Helper
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::
Refresh_Auxiliary_Structures_Helper()
{
    if(segment_list) Update_Segment_List();
    if(hierarchy) Initialize_Hierarchy();
}
//#####################################################################
// Function Update_Segment_List
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::
Update_Segment_List() // updates the segments assuming the particle positions are already updated
{
    if(!segment_list) segment_list=new ARRAY<typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT>(mesh.elements.m);
    for(int s=0;s<mesh.elements.m;s++){
        int i,j;mesh.elements(s).Get(i,j);
        (*segment_list)(s)=typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT(particles.X(i),particles.X(j));}
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::  
Initialize_Hierarchy(const bool update_boxes) // creates and updates the boxes as well
{
    delete hierarchy;
    if(segment_list) hierarchy=new SEGMENT_HIERARCHY<TV>(mesh,particles,*segment_list,update_boxes);
    else hierarchy=new SEGMENT_HIERARCHY<TV>(mesh,particles,update_boxes);
}
//#####################################################################
// Function Initialize_Straight_Mesh_And_Particles
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::  
Initialize_Straight_Mesh_And_Particles(const GRID<VECTOR<T,1> >& grid)
{
    Clean_Memory();
    mesh.Initialize_Straight_Mesh(grid.counts.x);
    particles.Preallocate(grid.counts.x);
    for(int i=0;i<grid.counts.x;i++) particles.X(particles.Add_Element()).x=grid.X(VECTOR<int,1>(i)).x;
}
//#####################################################################
// Function Initialize_Circle_Mesh_And_Particles
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::  
Initialize_Circle_Mesh_And_Particles(const int m,const T radius)
{
    Clean_Memory();
    mesh.Initialize_Straight_Mesh(m,true);
    particles.Add_Elements(m);
    for(int p=0;p<m;p++){
        COMPLEX<T> X=COMPLEX<T>::Polar(radius,(T)2*(T)pi/m*p);
        particles.X(p)=TV(VECTOR<T,2>(X.re,X.im));}
}
template<> void SEGMENTED_CURVE<VECTOR<double,1> >::  
Initialize_Circle_Mesh_And_Particles(const int m,const T radius)
{PHYSBAM_NOT_IMPLEMENTED();}
template<> void SEGMENTED_CURVE<VECTOR<float,1> >::  
Initialize_Circle_Mesh_And_Particles(const int m,const T radius)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Closest_Point_On_Curve
//#####################################################################
template<class TV> TV SEGMENTED_CURVE<TV>::  
Closest_Point_On_Curve(const TV& location,T thickness_over_two,int* closest_segment,T* distance) const
{
    T min_distance_squared=FLT_MAX;TV point;
    for(int i=0;i<mesh.elements.m;i++){
        T_SEGMENT segment=segment_list?(*segment_list)(i):Get_Element(i);
        TV new_point=segment.Closest_Point_On_Segment(location);
        T distance_squared=(new_point-location).Magnitude_Squared();
        if(distance_squared < min_distance_squared){
            min_distance_squared=distance_squared;point=new_point;
            if(closest_segment) *closest_segment=i;}}
    if(distance) *distance=sqrt(min_distance_squared);
    return point;
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void SEGMENTED_CURVE<TV>::
Rescale(const TV& scaling)
{
    for(int k=0;k<particles.Size();k++) particles.X(k)*=scaling;
    if(segment_list) Update_Segment_List();if(bounding_box) Update_Bounding_Box();
}
//#####################################################################
// Function Average_Edge_Length
//#####################################################################
template<class TV> typename TV::SCALAR SEGMENTED_CURVE<TV>::
Average_Edge_Length() const
{
    T average_edge_length=0;
    for(int s=0;s<mesh.elements.m;s++){
        int i,j;mesh.elements(s).Get(i,j);
        average_edge_length+=(particles.X(i)-particles.X(j)).Magnitude();}
    if(mesh.elements.m) average_edge_length/=mesh.elements.m;
    return average_edge_length;
}
//#####################################################################
// Function Total_Length
//#####################################################################
template<class TV> typename TV::SCALAR SEGMENTED_CURVE<TV>::
Total_Length() const
{
    T length=0;
    for(int t=0;t<mesh.elements.m;t++){
        int node1=mesh.elements(t)(0),node2=mesh.elements(t)(1);
        length+=(particles.X(node1)-particles.X(node2)).Magnitude();}
    return length;
}
template<class T> void Initialize_Boundary_Object_Helper(SEGMENTED_CURVE<VECTOR<T,1> >& sc)
{
    delete sc.point_simplices_1d;
    if(!sc.mesh.boundary_mesh) sc.mesh.Initialize_Boundary_Mesh();
    sc.point_simplices_1d=new POINT_SIMPLICES_1D<T>(*sc.mesh.boundary_mesh,sc.particles);
}
template<class T> void Initialize_Boundary_Object_Helper(SEGMENTED_CURVE<VECTOR<T,2> >& sc){PHYSBAM_NOT_IMPLEMENTED();}
template<class T> void Initialize_Boundary_Object_Helper(SEGMENTED_CURVE<VECTOR<T,3> >& sc){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Get_Boundary_Object
//#####################################################################
template<class TV> POINT_SIMPLICES_1D<typename TV::SCALAR>& SEGMENTED_CURVE<TV>::
Get_Boundary_Object()
{
    Initialize_Boundary_Object_Helper(*this);
    return *point_simplices_1d;
}
//#####################################################################
// Function Inside_Any_Simplex
//#####################################################################
template<class TV> bool SEGMENTED_CURVE<TV>::
Inside_Any_Simplex(const TV& location,int& segment_id,const T thickness_over_two) const
{
    assert(hierarchy);assert(segment_list);
    ARRAY<int> nearby_segments;hierarchy->Intersection_List(location,nearby_segments,thickness_over_two);
    for(int k=0;k<nearby_segments.m;k++){
        T_SEGMENT& segment=(*segment_list)(nearby_segments(k));
        if(segment.Inside(location,thickness_over_two)){segment_id=nearby_segments(k);return true;}}
    return false;
}
//#####################################################################
template class SEGMENTED_CURVE<VECTOR<float,1> >;
template class SEGMENTED_CURVE<VECTOR<float,2> >;
template class SEGMENTED_CURVE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SEGMENTED_CURVE<VECTOR<double,1> >;
template class SEGMENTED_CURVE<VECTOR<double,2> >;
template class SEGMENTED_CURVE<VECTOR<double,3> >;
#endif
}
