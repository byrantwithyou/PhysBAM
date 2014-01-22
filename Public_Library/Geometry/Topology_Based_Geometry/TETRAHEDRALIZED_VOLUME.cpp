//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ron Fedkiw, Eilene Hao, Geoffrey Irving, Neil Molino, Avi Robinson-Mosher, Mike Rodgers, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRALIZED_VOLUME
//#####################################################################
#include <Tools/Arrays/CONSTANT_ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/FRAME.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> TETRAHEDRALIZED_VOLUME<T>::
TETRAHEDRALIZED_VOLUME()
    :MESH_OBJECT<TV,TETRAHEDRON_MESH>(*new TETRAHEDRON_MESH,*new GEOMETRY_PARTICLES<TV>),tetrahedron_list(0),
    triangulated_surface(0),hierarchy(0),tetrahedron_volumes(0),nodal_volumes(0)
{
    this->need_destroy_mesh=this->need_destroy_particles=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> TETRAHEDRALIZED_VOLUME<T>::
TETRAHEDRALIZED_VOLUME(TETRAHEDRON_MESH& tetrahedron_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,TETRAHEDRON_MESH>(tetrahedron_mesh_input,particles_input),tetrahedron_list(0),
    triangulated_surface(0),hierarchy(0),tetrahedron_volumes(0),nodal_volumes(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TETRAHEDRALIZED_VOLUME<T>::
~TETRAHEDRALIZED_VOLUME()
{
    delete tetrahedron_list;
    delete triangulated_surface;
    delete hierarchy;
    delete tetrahedron_volumes;
    delete nodal_volumes;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Clean_Memory()
{
    MESH_OBJECT<TV,TETRAHEDRON_MESH>::Clean_Memory();
    delete tetrahedron_list;
    tetrahedron_list=0;
    delete triangulated_surface;
    triangulated_surface=0;
    delete hierarchy;
    hierarchy=0;
    delete tetrahedron_volumes;
    tetrahedron_volumes=0;
    delete nodal_volumes;
    nodal_volumes=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures_Helper
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Refresh_Auxiliary_Structures_Helper()
{
    if(tetrahedron_list) Update_Tetrahedron_List();
    if(hierarchy) Initialize_Hierarchy();
    if(triangulated_surface) Initialize_Triangulated_Surface();
    if(tetrahedron_volumes) Compute_Tetrahedron_Volumes();
    if(nodal_volumes) Compute_Nodal_Volumes();
}
//#####################################################################
// Function Update_Tetrahedron_List
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Update_Tetrahedron_List()
{
    if(!tetrahedron_list) tetrahedron_list=new ARRAY<TETRAHEDRON<T> >;
    tetrahedron_list->Resize(mesh.elements.m,false,false);
    for(int t=0;t<mesh.elements.m;t++){
        (*tetrahedron_list)(t).X=particles.X.Subset(mesh.elements(t));
        (*tetrahedron_list)(t).Create_Triangles();}
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Hierarchy(const bool update_boxes) // creates and updates the boxes as well
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    delete hierarchy;
    if(tetrahedron_list) hierarchy=new TETRAHEDRON_HIERARCHY<T>(mesh,particles,*tetrahedron_list,update_boxes);
    else hierarchy=new TETRAHEDRON_HIERARCHY<T>(mesh,particles,update_boxes);
}
//#####################################################################
// Function Initialize_Octahedron_Mesh_And_Particles
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Octahedron_Mesh_And_Particles(const GRID<TV>& grid)
{
    int i,j,k,m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z,particle=0;
    particles.Delete_All_Elements();
    mesh.Initialize_Octahedron_Mesh(m,n,mn);
    particles.Add_Elements(m*n*mn+(m+1)*(n+1)*(mn+1));
    for(k=0;k<mn;k++) for(j=0;j<n;j++) for(i=0;i<m;i++) particles.X(particle++)=grid.X(TV_INT(i,j,k));
    for(k=-1;k<mn;k++) for(j=-1;j<n;j++) for(i=-1;i<m;i++) particles.X(particle++)=grid.X(TV_INT(i,j,k))+(T).5*grid.dX;
}
//#####################################################################
// Function Initialize_Cube_Mesh_And_Particles
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Cube_Mesh_And_Particles(const GRID<TV>& grid)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z,particle=0;
    particles.Delete_All_Elements();
    mesh.Initialize_Cube_Mesh(m,n,mn);
    particles.Add_Elements(m*n*mn);
    for(int k=0;k<mn;k++) for(int j=0;j<n;j++) for(int i=0;i<m;i++) particles.X(particle++)=grid.X(TV_INT(i,j,k));
}
//#####################################################################
// Function Initialize_Prismatic_Cube_Mesh_And_Particles
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Prismatic_Cube_Mesh_And_Particles(const GRID<TV>& grid)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z,particle=0;
    particles.Delete_All_Elements();
    mesh.Initialize_Prismatic_Cube_Mesh(m,n,mn);
    particles.Add_Elements(m*n*mn);
    for(int k=0;k<mn;k++) for(int j=0;j<n;j++) for(int i=0;i<m;i++) particles.X(particle++)=grid.X(TV_INT(i,j,k));
}
//#####################################################################
// Function Initialize_Swept_Mesh_And_Particles
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Swept_Mesh_And_Particles(const TRIANGULATED_AREA<T>& ta,int layers,const FRAME<TV>& start_frame,const TV& sweep_offset)
{
    particles.Delete_All_Elements();
    particles.Add_Elements(ta.particles.X.m*(layers+1));
    TV dx=sweep_offset/layers;
    int m=ta.particles.X.m,n=layers*m;
    for(int i=0;i<m;i++) particles.X(i)=start_frame*TV(ta.particles.X(i));
    for(int i=0;i<n;i++) particles.X(i+m)=particles.X(i)+dx;
    mesh.Initialize_Swept_Mesh(ta.mesh,layers);
}
//#####################################################################
// Function Initialize_Cylinder_Mesh_And_Particles
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Cylinder_Mesh_And_Particles(const CYLINDER<T>& cyl,int num_elements_height,int num_elements_radius)
{
    TRIANGULATED_AREA<T>* ta=TESSELLATION::Generate_Triangles(SPHERE<VECTOR<T,2> >(VECTOR<T,2>(),cyl.radius),num_elements_radius);
    TV axis=cyl.plane2.x0-cyl.plane1.x0;
    Initialize_Swept_Mesh_And_Particles(*ta,num_elements_height,FRAME<TV>(cyl.plane1.x0,ROTATION<TV>::From_Rotated_Vector(TV(0,0,1),axis)),axis);
    delete ta;
}
//#####################################################################
// Function Check_Signed_Volumes_And_Make_Consistent
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Check_Signed_Volumes_And_Make_Consistent(bool verbose)
{
    for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        TV x0(particles.X(i)),x1(particles.X(j)),x2(particles.X(k)),x3(particles.X(l));
        T sign_of_volume=TV::Dot_Product(TV::Cross_Product(x1-x0,x2-x0),x3-x0); // left out division by 6
        if(sign_of_volume < 0){
            if(verbose) LOG::cout<<"tetrahedron number "<<t<<" is oriented improperly."<<std::endl;
            exchange(mesh.elements(t)(2),mesh.elements(t)(3));}}
    if(tetrahedron_list) Update_Tetrahedron_List();
}
//#####################################################################
// Function Initialize_Triangulated_Surface
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Initialize_Triangulated_Surface()
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    delete triangulated_surface;
    if(!mesh.boundary_mesh) mesh.Initialize_Boundary_Mesh();
    triangulated_surface=new TRIANGULATED_SURFACE<T>(*mesh.boundary_mesh,particles);
}
//#####################################################################
// Funcion Minimum_Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Volume(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Volume();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>::Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Signed_Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Signed_Volume(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Signed_Volume();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>::Signed_Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Total_Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Total_Volume() const
{
    T volume=0;
    for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        volume+=TETRAHEDRON<T>::Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));}
    return volume;
}
//#####################################################################
// Funcion Minimum_Angle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Angle(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Minimum_Angle();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Minimum_Angle();if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Maximum_Angle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Angle(int* index) const
{
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Maximum_Angle();if(temp > maximum){maximum=temp;t_save=t;}}
    else for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Maximum_Angle();if(temp > maximum){maximum=temp;t_save=t;}}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Minimum_Altitude
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Altitude(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Minimum_Altitude();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Minimum_Altitude();if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Minimum_Edge_Length
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Edge_Length(int* index) const
{
    int t_save=0;T minimum_squared=FLT_MAX;
    for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>::Minimum_Edge_Length_Squared(particles.X(i),particles.X(j),particles.X(k),particles.X(l));if(temp < minimum_squared){minimum_squared=temp;t_save=t;}}
    if(index) *index=t_save;
    return sqrt(minimum_squared);
}
//#####################################################################
// Funcion Maximum_Edge_Length
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Edge_Length(int* index) const
{
    int t_save=0;T maximum_squared=0;
    for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>::Maximum_Edge_Length_Squared(particles.X(i),particles.X(j),particles.X(k),particles.X(l));if(temp > maximum_squared){maximum_squared=temp;t_save=t;}}
    if(index) *index=t_save;
    return sqrt(maximum_squared);
}
//#####################################################################
// Funcion Maximum_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Aspect_Ratio(int* index) const
{
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}
    else for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Maximum_Interior_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Interior_Aspect_Ratio(int* index)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){if((*mesh.adjacent_elements)(t).m == 4){
        T temp=(*tetrahedron_list)(t).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}}
    else for(int t=0;t<mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m == 4){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Maximum_Boundary_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Boundary_Aspect_Ratio(int* index)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){if((*mesh.adjacent_elements)(t).m != 4){
        T temp=(*tetrahedron_list)(t).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}}
    else for(int t=0;t<mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m != 4){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();if(temp > maximum){maximum=temp;t_save=t;}}
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Average_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Average_Aspect_Ratio()
{
    int total=0;T sum=0;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){total++;sum+=(*tetrahedron_list)(t).Aspect_Ratio();}
    else for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        total++;sum+=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();}
    if(total != 0) sum/=total;
    return sum;
}
//#####################################################################
// Funcion Average_Interior_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Average_Interior_Aspect_Ratio()
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
    int total=0;T sum=0;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){if((*mesh.adjacent_elements)(t).m == 4){total++;sum+=(*tetrahedron_list)(t).Aspect_Ratio();}}
    else for(int t=0;t<mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m == 4){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        total++;sum+=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();}
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
    if(total != 0) sum/=total;
    return sum;
}
//#####################################################################
// Funcion Average_Boundary_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Average_Boundary_Aspect_Ratio()
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
    int total=0;T sum=0;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){if((*mesh.adjacent_elements)(t).m != 4){total++;sum+=(*tetrahedron_list)(t).Aspect_Ratio();}}
    else for(int t=0;t<mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m != 4){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        total++;sum+=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Aspect_Ratio();}
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
    if(total != 0) sum/=total;
    return sum;
}
//#####################################################################
// Funcion Minimum_Dihedral_Angle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Dihedral_Angle(int* index) const
{
    int t_save=0;T minimum=FLT_MAX;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Minimum_Dihedral_Angle();if(temp < minimum){minimum=temp;t_save=t;}}
    else for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Minimum_Dihedral_Angle();if(temp < minimum){minimum=temp;t_save=t;}}
    if(index) *index=t_save;
    return minimum;
}
//#####################################################################
// Funcion Maximum_Dihedral_Angle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Dihedral_Angle(int* index) const
{
    int t_save=0;T maximum=0;
    if(tetrahedron_list) for(int t=0;t<mesh.elements.m;t++){T temp=(*tetrahedron_list)(t).Maximum_Dihedral_Angle();if(temp > maximum){maximum=temp;t_save=t;}}
    else for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        T temp=TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Maximum_Dihedral_Angle();if(temp > maximum){maximum=temp;t_save=t;}}
    if(index) *index=t_save;
    return maximum;
}
//#####################################################################
// Funcion Maximum_Edge_Length
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Edge_Length(int* index)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool segment_mesh_defined=mesh.segment_mesh!=0;if(!segment_mesh_defined) mesh.Initialize_Segment_Mesh();
    int s_save=0;T maximum=0;
    for(int s=0;s<mesh.segment_mesh->elements.m;s++){int i,j;mesh.segment_mesh->elements(s).Get(i,j);
        T temp=(particles.X(i)-particles.X(j)).Magnitude();if(temp > maximum){maximum=temp;s_save=s;}}
    if(index) *index=s_save;
    if(!segment_mesh_defined){delete mesh.segment_mesh;mesh.segment_mesh=0;}
    return maximum;
}
//#####################################################################
// Funcion Minimum_Edge_Length
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Minimum_Edge_Length(int* index)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool segment_mesh_defined=mesh.segment_mesh!=0;if(!segment_mesh_defined) mesh.Initialize_Segment_Mesh();
    int s_save=0;T minimum=FLT_MAX;
    for(int s=0;s<mesh.segment_mesh->elements.m;s++){int i,j;mesh.segment_mesh->elements(s).Get(i,j);
        T temp=(particles.X(i)-particles.X(j)).Magnitude();if(temp < minimum){minimum=temp;s_save=s;}}
    if(index) *index=s_save;
    if(!segment_mesh_defined){delete mesh.segment_mesh;mesh.segment_mesh=0;}
    return minimum;
}
//#####################################################################
// Function Interior_Laplacian_Smoothing
//#####################################################################
// one step of mesh adjustment using Laplacian smoothing
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Advance_Interior_Laplacian_Smoothing()
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool neighbor_nodes_defined=mesh.neighbor_nodes!=0;if(!neighbor_nodes_defined) mesh.Initialize_Neighbor_Nodes();
    bool node_on_boundary_defined=mesh.node_on_boundary!=0;if(!node_on_boundary_defined) mesh.Initialize_Node_On_Boundary();
    // compute the centroid of the neighbors - if on the boundary, just use boundary neighbors
    for(int i=0;i<mesh.neighbor_nodes->m;i++){
        int number=0;TV target;
        for(int j=0;j<(*mesh.neighbor_nodes)(i).m;j++){
            int node=(*mesh.neighbor_nodes)(i)(j);
            if(!(*mesh.node_on_boundary)(i) || (*mesh.node_on_boundary)(node)){number++;target+=particles.X(node);}}
        if(number != 0){target/=(T)number;particles.X(i)=target;}}
    if(!neighbor_nodes_defined){delete mesh.neighbor_nodes;mesh.neighbor_nodes=0;}
    if(!node_on_boundary_defined){delete mesh.node_on_boundary;mesh.node_on_boundary=0;}
}
//#####################################################################
// Function Centroid
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRALIZED_VOLUME<T>::
Centroid(const int tetrahedron) const
{
    int i,j,k,l;mesh.elements(tetrahedron).Get(i,j,k,l);
    return (T).25*(particles.X(i)+particles.X(j)+particles.X(k)+particles.X(l));
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Rescale(const T scaling_factor)
{
    Rescale(scaling_factor,scaling_factor,scaling_factor);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Rescale(const T scaling_x,const T scaling_y,const T scaling_z)
{
    for(int k=0;k<particles.Size();k++) particles.X(k)*=TV(scaling_x,scaling_y,scaling_z);
}
//#####################################################################
// Function Signed_Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Signed_Volume(const int tetrahedron) const
{
    int i,j,k,l;mesh.elements(tetrahedron).Get(i,j,k,l);
    return TETRAHEDRON<T>::Signed_Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));
}
//#####################################################################
// Function Volume
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Volume(const int tetrahedron) const
{
    int i,j,k,l;mesh.elements(tetrahedron).Get(i,j,k,l);
    return TETRAHEDRON<T>::Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));
}
//#####################################################################
// Function Centroid_Of_Neighbors
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRALIZED_VOLUME<T>::
Centroid_Of_Neighbors(const int node) const
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool neighbor_nodes_defined=mesh.neighbor_nodes!=0;if(!neighbor_nodes_defined) mesh.Initialize_Neighbor_Nodes();
    TV target;
    int number_of_neighbors=(*mesh.neighbor_nodes)(node).m;
    for(int j=0;j<number_of_neighbors;j++) target+=particles.X((*mesh.neighbor_nodes)(node)(j));
    if(number_of_neighbors != 0) target/=(T)number_of_neighbors;
    else target=particles.X(node); // if no neighbors, return the current node location
    if(!neighbor_nodes_defined){delete mesh.neighbor_nodes;mesh.neighbor_nodes=0;}
    return target;
}
//#####################################################################
// Function Completely_Inside_Box
//#####################################################################
template<class T> bool TETRAHEDRALIZED_VOLUME<T>::
Completely_Inside_Box(const int tetrahedron,const RANGE<TV>& box) const
{
    int i,j,k,l;mesh.elements(tetrahedron).Get(i,j,k,l);
    return (box.Lazy_Inside(particles.X(i)) && box.Lazy_Inside(particles.X(j)) && box.Lazy_Inside(particles.X(k)) && box.Lazy_Inside(particles.X(l)));
}
//#####################################################################
// Function Discard_Spikes_From_Adjacent_Elements
//#####################################################################
// throws out all tetrahedrons with only one neighbor (i.e. a spike on the boundary)
// returns index of first discarded
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Spikes_From_Adjacent_Elements(ARRAY<int>* deletion_list)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
    if(deletion_list){
        deletion_list->Resize(0);
        for(int t=0;t<mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m == 1) deletion_list->Append(t);
        mesh.Delete_Sorted_Elements(*deletion_list);}
    else{
        ARRAY<int> list;
        for(int t=0;t<mesh.elements.m;t++) if((*mesh.adjacent_elements)(t).m == 1) list.Append(t);
        mesh.Delete_Sorted_Elements(list);}
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
}
//#####################################################################
// Function Interior_Edges_With_Boundary_Nodes
//#####################################################################
// throws out all tetrahedrons with all four nodes on boundary
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Interior_Edges_With_Boundary_Nodes(ARRAY<VECTOR<int,2> >* deletion_list)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool neighbor_nodes_defined=mesh.neighbor_nodes!=0;if(!neighbor_nodes_defined) mesh.Initialize_Neighbor_Nodes();
    bool node_on_boundary_defined=mesh.node_on_boundary!=0;if(!node_on_boundary_defined) mesh.Initialize_Node_On_Boundary();
    bool boundary_nodes_defined=mesh.boundary_nodes!=0;if(!boundary_nodes_defined) mesh.Initialize_Boundary_Nodes();
    bool boundary_mesh_defined=mesh.boundary_mesh!=0;if(!boundary_mesh_defined) mesh.Initialize_Boundary_Mesh();
    bool boundary_mesh_segment_mesh_defined=mesh.boundary_mesh->segment_mesh!=0;
    if(!boundary_mesh_segment_mesh_defined) mesh.boundary_mesh->Initialize_Segment_Mesh();
    bool segment_mesh_defined=mesh.segment_mesh!=0;if(!segment_mesh_defined) mesh.Initialize_Segment_Mesh();
    bool boundary_mesh_segment_mesh_incident_elements_defined=mesh.boundary_mesh->segment_mesh->incident_elements!=0;
    if(!boundary_mesh_segment_mesh_incident_elements_defined) mesh.boundary_mesh->segment_mesh->Initialize_Incident_Elements();

    deletion_list->Resize(0);
    for(int t=0;t<mesh.boundary_nodes->m;t++){
        int node1=(*mesh.boundary_nodes)(t);
        for(int i=0;i<(*mesh.neighbor_nodes)(node1).m;i++){
            int node2=(*mesh.neighbor_nodes)(node1)(i);
            if(node1 <node2 && (*mesh.node_on_boundary)(node2) && !mesh.boundary_mesh->segment_mesh->Segment(node1,node2))
                deletion_list->Append(VECTOR<int,2>(node1,node2));}}

    // delete node_on_boundary if defined in this function
    if(!boundary_mesh_segment_mesh_incident_elements_defined){
        delete mesh.boundary_mesh->segment_mesh->incident_elements;mesh.boundary_mesh->segment_mesh->incident_elements=0;}
    if(!segment_mesh_defined){delete mesh.segment_mesh;mesh.segment_mesh=0;}
    if(!boundary_mesh_segment_mesh_defined){delete mesh.boundary_mesh->segment_mesh;mesh.boundary_mesh->segment_mesh=0;}
    if(!boundary_mesh_defined){delete mesh.boundary_mesh;mesh.boundary_mesh=0;}
    if(!boundary_nodes_defined){delete mesh.boundary_nodes;mesh.boundary_nodes=0;}
    if(!node_on_boundary_defined){delete mesh.node_on_boundary;mesh.node_on_boundary=0;}
    if(!neighbor_nodes_defined){delete mesh.neighbor_nodes;mesh.neighbor_nodes=0;}
}
//#####################################################################
// Function Discard_Spikes
//#####################################################################
// throws out all tetrahedrons with all four nodes on boundary
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Spikes(ARRAY<int>* deletion_list)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool node_on_boundary_defined=mesh.node_on_boundary!=0;if(!node_on_boundary_defined) mesh.Initialize_Node_On_Boundary();
    if(deletion_list){
        deletion_list->Resize(0);
        for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
            if((*mesh.node_on_boundary)(i) && (*mesh.node_on_boundary)(j) && (*mesh.node_on_boundary)(k) && (*mesh.node_on_boundary)(l))
                deletion_list->Append(t);}
        mesh.Delete_Sorted_Elements(*deletion_list);}
    else{
        ARRAY<int> list;
        for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
            if((*mesh.node_on_boundary)(i) && (*mesh.node_on_boundary)(j) && (*mesh.node_on_boundary)(k) && (*mesh.node_on_boundary)(l))
                list.Append(t);}
        mesh.Delete_Sorted_Elements(list);}
    if(!node_on_boundary_defined){delete mesh.node_on_boundary;mesh.node_on_boundary=0;}
}
//#####################################################################
// Function Inverted_Tetrahedrons
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Inverted_Tetrahedrons(ARRAY<int>& inverted_tetrahedrons) const
{
    inverted_tetrahedrons.Resize(0);
    for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        TV x0(particles.X(i)),x1(particles.X(j)),x2(particles.X(k)),x3(particles.X(l));
        T sign_of_volume=TV::Dot_Product(TV::Cross_Product(x1-x0,x2-x0),x3-x0); // left out division by 6
        if(sign_of_volume < 0) inverted_tetrahedrons.Append(t);}
}
//#####################################################################
// Function Inside
//#####################################################################
// note: return the first tetrahedron that it is inside of (including boundary), otherwise returns 0
template<class T> int TETRAHEDRALIZED_VOLUME<T>::
Inside(const TV& location,const T thickness_over_two) const
{
    if(!hierarchy) PHYSBAM_FATAL_ERROR();
    if(hierarchy->box_hierarchy(hierarchy->root).Outside(location,thickness_over_two)) return -1;
    ARRAY<int> tetrahedrons_to_check;hierarchy->Intersection_List(location,tetrahedrons_to_check,thickness_over_two);
    if(tetrahedron_list) for(int p=0;p<tetrahedrons_to_check.m;p++){
        int t=tetrahedrons_to_check(p);if(!(*tetrahedron_list)(t).Outside(location,thickness_over_two)) return t;}
    else for(int p=0;p<tetrahedrons_to_check.m;p++){
        int t=tetrahedrons_to_check(p);int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        TETRAHEDRON<T> tetrahedron_to_check(particles.X(i),particles.X(j),particles.X(k),particles.X(l));
        if(!tetrahedron_to_check.Outside(location,thickness_over_two)) return t;}
    return -1;
}
//#####################################################################
// Function Inside_Any_Simplex
//#####################################################################
template<class T> bool TETRAHEDRALIZED_VOLUME<T>::
Inside_Any_Simplex(const TV& location,int& tetrahedron_id,const T thickness_over_two) const
{
    assert(hierarchy);assert(tetrahedron_list);
    ARRAY<int> nearby_tetrahedrons;
    hierarchy->Intersection_List(location,nearby_tetrahedrons,thickness_over_two);
    for(int k=0;k<nearby_tetrahedrons.m;k++){
        TETRAHEDRON<T>& tetrahedron=(*tetrahedron_list)(nearby_tetrahedrons(k));
        if(tetrahedron.Inside(location,thickness_over_two)){tetrahedron_id=nearby_tetrahedrons(k);return true;}}
    return false;
}
//#####################################################################
// Function Inside
//#####################################################################
// note: return the first tetrahedron that it is inside of.  If none, then the first including boundary.  Otherwise returns 0
template<class T> int TETRAHEDRALIZED_VOLUME<T>::
Find(const TV& location,const T thickness_over_two) const
{
    ARRAY<int> tetrahedrons_to_check;
    return Find(location,thickness_over_two,tetrahedrons_to_check);
}
//#####################################################################
// Function Inside
//#####################################################################
// note: return the first tetrahedron that it is inside of.  If none, then the first including boundary.  Otherwise returns 0
template<class T> int TETRAHEDRALIZED_VOLUME<T>::
Find(const TV& location,const T thickness_over_two,ARRAY<int>& scratch) const
{
    if(!tetrahedron_list || !hierarchy) PHYSBAM_FATAL_ERROR();
    scratch.Remove_All();
    hierarchy->Intersection_List(location,scratch,thickness_over_two);
    for(int p=0;p<scratch.m;p++){int t=scratch(p);if(!(*tetrahedron_list)(t).Outside(location,0)) return t;}
    for(int p=0;p<scratch.m;p++){int t=scratch(p);if(!(*tetrahedron_list)(t).Outside(location,thickness_over_two)) return t;}
    return -1;
}
//#####################################################################
// Function Discard_Tetrahedrons_Outside_Implicit_Surface
//#####################################################################
// uses Whitney-like criterion to discard only those tets that are for sure outside levelset (assuming accurate signed distance)
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Tetrahedrons_Outside_Implicit_Surface(IMPLICIT_OBJECT<TV>& implicit_surface)
{
    for(int t=mesh.elements.m-1;t>=0;t--){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        TV xi=particles.X(i),xj=particles.X(j),xk=particles.X(k),xl=particles.X(l);
        T max_length=TETRAHEDRON<T>::Maximum_Edge_Length(xi,xj,xk,xl);
        T min_phi=min(implicit_surface(xi),implicit_surface(xj),implicit_surface(xk),implicit_surface(xl));
        if(min_phi>max_length) mesh.elements.Remove_Index_Lazy(t);}
}
//#####################################################################
// Function Discard_Tetrahedrons_Outside_Implicit_Surface_Agressive
//#####################################################################
// uses Whitney-like criterion to discard only those tets that are not for sure inside levelset (assuming accurate signed distance)
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(IMPLICIT_OBJECT<TV>& implicit_surface)
{
    for(int t=mesh.elements.m-1;t>=0;t--){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        TV xi=particles.X(i),xj=particles.X(j),xk=particles.X(k),xl=particles.X(l);
        T max_length=TETRAHEDRON<T>::Maximum_Edge_Length(xi,xj,xk,xl);
        T max_phi=max(implicit_surface(xi),implicit_surface(xj),implicit_surface(xk),implicit_surface(xl));
        if(max_phi+max_length>0) mesh.elements.Remove_Index_Lazy(t);}
}
//#####################################################################
// Function Discard_Tetrahedrons_Outside_Implicit_Surface_Agressive
//#####################################################################
// same as above Discard_Tetrahedrons_Outside_Implicit_Surface_Agressive but restricted to tets fully contained inside one of the supplied boxes
// TODO: Need to decide what constitutes good criteria for tet to be contained in a box: for now demanding full containment, i.e. all four particles
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(IMPLICIT_OBJECT<TV>& implicit_surface,const ARRAY<RANGE<TV> >& bounding_boxes)
{
    for(int b=0;b<bounding_boxes.Size();b++){const RANGE<TV>& box=bounding_boxes(b);
        for(int t=mesh.elements.m-1;t>=0;t--){
            int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
            TV xi=particles.X(i),xj=particles.X(j),xk=particles.X(k),xl=particles.X(l);
            if(box.Lazy_Outside(xi)||box.Lazy_Outside(xj)||box.Lazy_Outside(xk)||box.Lazy_Outside(xl)) continue;
            T max_length=TETRAHEDRON<T>::Maximum_Edge_Length(xi,xj,xk,xl);
            T max_phi=max(implicit_surface(xi),implicit_surface(xj),implicit_surface(xk),implicit_surface(xl));
            if(max_phi+max_length>0) mesh.elements.Remove_Index_Lazy(t);}}
}
//#####################################################################
// Function Maximum_Magnitude_Phi_On_Boundary
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Maximum_Magnitude_Phi_On_Boundary(const IMPLICIT_OBJECT<TV>& implicit_surface,int* index)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool node_on_boundary_defined=mesh.node_on_boundary!=0;if(!node_on_boundary_defined) mesh.Initialize_Node_On_Boundary();
    T phi=0,max_phi=0;int p_save=0;
    for(int p=0;p<particles.Size();p++) if((*mesh.node_on_boundary)(p)){phi=abs(implicit_surface(particles.X(p)));if(phi > max_phi){max_phi=phi;p_save=p;}}
    if(index) *index=p_save;
    if(!node_on_boundary_defined){delete mesh.node_on_boundary;mesh.node_on_boundary=0;}
    return max_phi;
}
//#####################################################################
// Function Volume_Incident_On_A_Particle
//#####################################################################
template<class T> T TETRAHEDRALIZED_VOLUME<T>::
Volume_Incident_On_A_Particle(const int particle_index)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool incident_elements_defined=mesh.incident_elements!=0;if(!incident_elements_defined) mesh.Initialize_Incident_Elements();
    T total_incident_volume=0;
    for(int t=0;t<(*mesh.incident_elements)(particle_index).m;t++){int i,j,k,l;mesh.elements((*mesh.incident_elements)(particle_index)(t)).Get(i,j,k,l);
        total_incident_volume+=TETRAHEDRON<T>::Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));}
    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
    return total_incident_volume;
}
//#####################################################################
// Function Split_Along_Fracture_Plane
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Split_Along_Fracture_Plane(const PLANE<T>& plane,ARRAY<int>& particle_replicated)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    if(particle_replicated.m != particles.Size()) particle_replicated.Resize(particles.Size());
    bool incident_elements_defined=mesh.incident_elements!=0;if(!incident_elements_defined) mesh.Initialize_Incident_Elements();
    ARRAY<bool> positive_side(mesh.elements.m);int number_on_positive_side=0;
    int t;for(t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        TV x0=particles.X(i),x1=particles.X(j),x2=particles.X(k),x3=particles.X(l),centroid=(T).25*(x0+x1+x2+x3);
        if(plane.Signed_Distance(centroid) >= 0){positive_side(t)=true;number_on_positive_side++;}}
    int p;for(p=0;p<particles.Size();p++){
        bool seen_positive_side=false,seen_negative_side=false;
        for(t=0;t<(*mesh.incident_elements)(p).m;t++){
            if(positive_side((*mesh.incident_elements)(p)(t))) seen_positive_side=true;else seen_negative_side=true;}
        if(seen_positive_side && seen_negative_side) particle_replicated(p)=1;}
    int number_of_new_particles=0;
    for(p=0;p<particles.Size();p++) if(particle_replicated(p)){ // assumes we're storing mass (this is not set here!), position, & velocity
        int new_index=particles.Add_Element();particle_replicated(p)=new_index;number_of_new_particles++;
        particles.X(new_index)=particles.X(p);particles.V(new_index)=particles.V(p);}
    // loop through tets and change indices in negative_side to new_indices (from particle_replicated)
    for(t=0;t<mesh.elements.m;t++) if(!positive_side(t)){ int i,j,k,l;mesh.elements(t).Get(i,j,k,l); // replace indices with replicated_indices
        if(particle_replicated(i)) i=particle_replicated(i);if(particle_replicated(j)) j=particle_replicated(j);
        if(particle_replicated(k)) k=particle_replicated(k);if(particle_replicated(l)) l=particle_replicated(l);
        mesh.elements(t).Set(i,j,k,l);}
    mesh.number_nodes=particles.Size();
    if(incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
}
//#####################################################################
// Function Split_Node
//#####################################################################
// warning: will corrupt any aux structures aside from incident_elements
template<class T> int TETRAHEDRALIZED_VOLUME<T>::
Split_Node(const int particle_index,const TV& normal)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    bool incident_elements_defined=mesh.incident_elements!=0;if(!incident_elements_defined) mesh.Initialize_Incident_Elements();
    PLANE<T> plane(normal,particles.X(particle_index));ARRAY<int> tets_incident_on_old_particle,tets_incident_on_new_particle;
    int t;for(t=0;t<(*mesh.incident_elements)(particle_index).m;t++){
        int this_incident_tet=(*mesh.incident_elements)(particle_index)(t);int i,j,k,l;mesh.elements(this_incident_tet).Get(i,j,k,l);
        TV x0=particles.X(i),x1=particles.X(j),x2=particles.X(k),x3=particles.X(l),centroid=(T).25*(x0+x1+x2+x3);
        if(plane.Signed_Distance(centroid) < 0) tets_incident_on_new_particle.Append(this_incident_tet);
        else tets_incident_on_old_particle.Append(this_incident_tet);}
    int new_particle=0;
    if(tets_incident_on_old_particle.m != 0 && tets_incident_on_new_particle.m != 0){
        // new particle - assumes we're storing position, and velocity - user must fix mass outside this function call
        new_particle=particles.Add_Element();mesh.number_nodes=particles.Size();
        particles.X(new_particle)=particles.X(particle_index);particles.V(new_particle)=particles.V(particle_index);
        for(t=0;t<(*mesh.incident_elements)(particle_index).m;t++){
            int this_incident_tet=(*mesh.incident_elements)(particle_index)(t);int i,j,k,l;mesh.elements(this_incident_tet).Get(i,j,k,l);
            TV x0=particles.X(i),x1=particles.X(j),x2=particles.X(k),x3=particles.X(l),centroid=(T).25*(x0+x1+x2+x3);
            if(plane.Signed_Distance(centroid) < 0){ // relabel with duplicate node
                if(i == particle_index) i=new_particle;if(j == particle_index) j=new_particle;if(k == particle_index) k=new_particle;if(l == particle_index) l=new_particle;
                mesh.elements(this_incident_tet).Set(i,j,k,l);}}
        if(incident_elements_defined){ //repair incident tetrahedrons if necessary
            (*mesh.incident_elements)(particle_index).Clean_Memory();
            (*mesh.incident_elements)(particle_index).Append_Elements(tets_incident_on_old_particle);
            (*mesh.incident_elements).Append(tets_incident_on_new_particle);}}
    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
    return new_particle;
}
//#####################################################################
// Function Compute_Tetrahedron_Volumes
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Compute_Tetrahedron_Volumes()
{
    if(!tetrahedron_volumes) tetrahedron_volumes=new ARRAY<T>;
    tetrahedron_volumes->Resize(mesh.elements.m,false,false);
    for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        (*tetrahedron_volumes)(t)=TETRAHEDRON<T>::Signed_Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));}
}
//#####################################################################
// Function Compute_Nodal_Volumes
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Compute_Nodal_Volumes(bool save_tetrahedron_volumes)
{
    if(!nodal_volumes) nodal_volumes=new ARRAY<T>;
    *nodal_volumes=CONSTANT_ARRAY<T>(particles.Size(),0);
    if(save_tetrahedron_volumes){
        if(!tetrahedron_volumes) tetrahedron_volumes=new ARRAY<T>;
        tetrahedron_volumes->Resize(mesh.elements.m);
        for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
            (*tetrahedron_volumes)(t)=TETRAHEDRON<T>::Signed_Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));
            T volume=(T).25*(*tetrahedron_volumes)(t);
            (*nodal_volumes)(i)+=volume;(*nodal_volumes)(j)+=volume;(*nodal_volumes)(k)+=volume;(*nodal_volumes)(l)+=volume;}}
    else
        for(int t=0;t<mesh.elements.m;t++){int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
            T volume=(T).25*TETRAHEDRON<T>::Signed_Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l));
            (*nodal_volumes)(i)+=volume;(*nodal_volumes)(j)+=volume;(*nodal_volumes)(k)+=volume;(*nodal_volumes)(l)+=volume;}
}
//#####################################################################
// Function Print_Statistics
//#####################################################################
template<class T> void TETRAHEDRALIZED_VOLUME<T>::
Print_Statistics(std::ostream& output)
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR();
    Update_Bounding_Box();
    if(!mesh.segment_mesh) mesh.Initialize_Segment_Mesh();
    if(!mesh.incident_elements) mesh.Initialize_Incident_Elements();
    if(!mesh.adjacent_elements) mesh.Initialize_Adjacent_Elements();
    if(!mesh.boundary_mesh) mesh.Initialize_Boundary_Mesh();

    output<<"tetrahedrons = "<<mesh.elements.m<<std::endl;
    output<<"particles = "<<particles.Size()<<std::endl;
    {int particles_touched=0;for(int p=0;p<particles.Size();p++) if((*mesh.incident_elements)(p).m) particles_touched++;
    output<<"particles touched = "<<particles_touched<<std::endl;}
    output<<"bounding box = "<<*bounding_box<<std::endl;
    if(particles.store_velocity){
        int index=particles.V.Arg_Maximum_Magnitude();
        output<<"max_speed = "<<particles.V(index).Magnitude()<<" ("<<index<<")"<<std::endl;}
    int index;
    output<<"total volume = "<<Total_Volume()<<std::endl;
    output<<"max_aspect_ratio = "<<Maximum_Aspect_Ratio(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"max_boundary_aspect_ratio = "<<Maximum_Boundary_Aspect_Ratio(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"max_interior_aspect_ratio = "<<Maximum_Interior_Aspect_Ratio(&index);output<<" ("<<index<<")"<<std::endl;
    output<<"avg_boundary_aspect_ratio = "<<Average_Boundary_Aspect_Ratio()<<std::endl;
    output<<"avg_interior_aspect_ratio = "<<Average_Interior_Aspect_Ratio()<<std::endl;
    output<<"min_volume = "<<Minimum_Volume(&index)<<std::endl;
    output<<"min_angle = "<<180/pi*Minimum_Angle()<<std::endl;
    output<<"max_angle = "<<180/pi*Maximum_Angle()<<std::endl;
    output<<"min_dihedral_angle = "<<180/pi*Minimum_Dihedral_Angle()<<std::endl;
    output<<"max_dihedral_angle = "<<180/pi*Maximum_Dihedral_Angle()<<std::endl;
    output<<"min_edge_length = "<<Minimum_Edge_Length()<<std::endl;
    output<<"max_edge_length = "<<Maximum_Edge_Length()<<std::endl;
    output<<"min_altitude = "<<Minimum_Altitude()<<std::endl;

    ARRAY<int> nonmanifold_nodes;mesh.boundary_mesh->Non_Manifold_Nodes(nonmanifold_nodes);
    output<<nonmanifold_nodes.m<<" nonmanifold nodes = "<<nonmanifold_nodes;
}
//#####################################################################
namespace PhysBAM{
template class TETRAHEDRALIZED_VOLUME<float>;
template class TETRAHEDRALIZED_VOLUME<double>;
}
