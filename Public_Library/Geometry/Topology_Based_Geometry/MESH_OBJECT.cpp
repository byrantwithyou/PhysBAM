//#####################################################################
// Copyright 2004-2009, Geoffrey Irving, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/UNION_FIND.h>
#include <Core/Math_Tools/FACTORIAL.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Utilities/DEBUG_CAST.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_PARTITION.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T_MESH> MESH_OBJECT<TV,T_MESH>::
MESH_OBJECT(T_MESH& mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :mesh(mesh_input),particles(particles_input),bounding_box(0),particle_partition(0),need_destroy_mesh(false),need_destroy_particles(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T_MESH> MESH_OBJECT<TV,T_MESH>::
~MESH_OBJECT()
{
    if(need_destroy_mesh) delete &mesh;
    if(need_destroy_particles) delete &particles;
    delete bounding_box;
    delete particle_partition;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Clean_Memory()
{
    delete bounding_box;bounding_box=0;
    delete particle_partition;particle_partition=0;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,class T_MESH> typename MESH_TO_OBJECT<TV,T_MESH>::TYPE* MESH_OBJECT<TV,T_MESH>::
Create()
{
    T_DERIVED_OBJECT* object=new T_DERIVED_OBJECT(*(new T_MESH),*(new GEOMETRY_PARTICLES<TV>));
    object->need_destroy_mesh=object->need_destroy_particles=true;return object;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,class T_MESH> typename MESH_TO_OBJECT<TV,T_MESH>::TYPE* MESH_OBJECT<TV,T_MESH>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    T_DERIVED_OBJECT* object=new T_DERIVED_OBJECT(*(new T_MESH),particles);
    object->need_destroy_mesh=true;return object;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Refresh_Auxiliary_Structures()
{
    mesh.Refresh_Auxiliary_Structures();
    Refresh_Auxiliary_Structures_Helper();
    if(bounding_box) Update_Bounding_Box(); // this can depend on the hierarchy which by derived helper
    if(particle_partition) Initialize_Particle_Partition(particle_partition->grid.counts);
}
//#####################################################################
// Function Update_Number_Nodes
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Update_Number_Nodes()
{
    mesh.Set_Number_Nodes(particles.Size());
}
//#####################################################################
// Function Update_Bounding_Box
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Update_Bounding_Box()
{
    if(!bounding_box) bounding_box=new RANGE<TV>();
    const BOX_HIERARCHY<TV>* hierarchy=&*debug_cast<const typename MESH_TO_OBJECT<TV,T_MESH>::TYPE&>(*this).hierarchy;
    if(!mesh.elements.m) *bounding_box=RANGE<TV>::Bounding_Box(particles.X);
    else if(hierarchy) *bounding_box=hierarchy->box_hierarchy(hierarchy->root);
    else *bounding_box=RANGE<TV>::Bounding_Box(particles.X.Subset(mesh.elements.Flattened()));
}
//#####################################################################
// Function Initialize_Particle_Partition
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Initialize_Particle_Partition(const VECTOR<int,TV::m>& counts)
{
    PHYSBAM_ASSERT(bounding_box);VECTOR<int,TV::m> counts_new;
    for(int i=0;i<counts.m;i++) counts_new[i]=desired_particle_partition_counts[i]?desired_particle_partition_counts[i]:counts[i];
    PHYSBAM_ASSERT(counts_new.All_Greater(VECTOR<int,TV::m>()));
    delete particle_partition;particle_partition=new PARTICLE_PARTITION<TV>(*bounding_box,counts_new,particles);
}
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV,class T_MESH> STRUCTURE<TV>* MESH_OBJECT<TV,T_MESH>::
Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const // number_nodes must be set elsewhere
{
    typename MESH_TO_OBJECT<TV,T_MESH>::TYPE* object=Create(new_particles);
    int offset=new_particles.Size();
    new_particles.Append(particles);
    if(particle_indices) for(int p=0;p<particles.Size();p++) particle_indices->Append(p+offset);
    object->mesh.Initialize_Mesh_With_Particle_Offset(mesh,offset);
    return object;
}
//#####################################################################
// Function Discard_Valence_Zero_Particles_And_Renumber
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Discard_Valence_Zero_Particles_And_Renumber(ARRAY<int>& condensation_mapping)
{
    // mark which nodes are used
    ARRAY<bool> node_is_used(mesh.number_nodes);
    node_is_used.Subset(mesh.elements.Flattened()).Fill(true);
    
    // make condensation mapping
    condensation_mapping.Resize(mesh.number_nodes,init_all,-1);
    int counter=0;
    for(int t=0;t<mesh.number_nodes;t++) if(node_is_used(t)) condensation_mapping(t)=counter++;
    
    // make new triangle mesh
    mesh.number_nodes=counter;
    for(int t=0;t<mesh.elements.m;t++)
        mesh.elements(t)=condensation_mapping.Subset(mesh.elements(t));
    
    // do particles same way
    for(int p=0;p<condensation_mapping.m;p++) if(condensation_mapping(p)<0) particles.Add_To_Deletion_List(p);
    for(int p=condensation_mapping.m;p<particles.Size();p++) particles.Add_To_Deletion_List(p);
    particles.Delete_Elements_On_Deletion_List(true);
    particles.Compact();

    Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Discard_Valence_Zero_Particles_And_Renumber
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Merge_Overlapping_Particles_And_Renumber(ARRAY<int>& condensation_mapping,T tol)
{
    UNION_FIND<> uf(particles.X.m);
    PARTICLE_HIERARCHY<TV> ph(particles.X,true,0);
    ARRAY<int> intersection_list;
    for(int i=0;i<particles.X.m;i++){
        intersection_list.Remove_All();
        ph.Intersection_List(particles.X(i),intersection_list,tol);
        for(int j=0;j<intersection_list.m;j++)
            if((particles.X(i)-particles.X(intersection_list(j))).Magnitude_Squared()>=tol*tol)
                uf.Union(i,intersection_list(j));}

    int next=0;
    condensation_mapping.Resize(particles.X.m,init_all,-1);
    for(int i=0;i<particles.X.m;i++)
        if(condensation_mapping(i)==-1){
            int p=uf.Find(i);
            if(condensation_mapping(p)==-1){
                particles.X(next)=particles.X(i);
                condensation_mapping(p)=next++;}
            condensation_mapping(i)=condensation_mapping(p);}

    for(int p=condensation_mapping.m;p<particles.Size();p++)
        particles.Add_To_Deletion_List(p);
    particles.Delete_Elements_On_Deletion_List(true);

    // make new triangle mesh
    mesh.number_nodes=next;
    mesh.elements.Flattened()=condensation_mapping.Subset(mesh.elements.Flattened());
    Update_Number_Nodes();

    Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Union_Mesh_Objects_Relatively
//#####################################################################
template<class TV,class T_MESH> typename MESH_OBJECT<TV,T_MESH>::T_DERIVED_OBJECT* MESH_OBJECT<TV,T_MESH>::
Union_Mesh_Objects_Relatively(const ARRAY<T_DERIVED_OBJECT*>& object_list,const ARRAY<FRAME<TV> >& relative_frames)
{
    T_DERIVED_OBJECT* object=Create();
    Union_Mesh_Objects_Relatively(object,object_list,relative_frames);
    return object;
}
//#####################################################################
// Function Union_Mesh_Objects_Relatively
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Union_Mesh_Objects_Relatively(T_DERIVED_OBJECT *object,const ARRAY<T_DERIVED_OBJECT*>& object_list,const ARRAY<FRAME<TV> >& relative_frames)
{
    GEOMETRY_PARTICLES<TV>& particles=object->particles;
    particles.Clean_Memory();
    object->mesh.elements.Remove_All();
    // resize
    {int total_particles=0,total_elements=0;
    for(int i=0;i<object_list.m;i++){
        total_particles+=object_list(i)->particles.Size();
        total_elements+=object_list(i)->mesh.elements.m;}
    particles.Preallocate(total_particles);object->mesh.elements.Preallocate(total_elements);}
    // copy
    for(int i=0;i<object_list.m;i++){
        int particle_offset=particles.Size();
        particles.Add_Arrays(object_list(i)->particles);
        particles.Append(object_list(i)->particles);
        for(int p=0;p<object_list(i)->particles.Size();p++){int p2=p+particle_offset;
            particles.X(p2)=relative_frames(i)*particles.X(p2);
            if(particles.store_velocity) particles.V(p2)=relative_frames(i).r.Rotate(particles.V(p2));}
        object->mesh.elements.Append_Elements(object_list(i)->mesh.elements+particle_offset);}
    object->Update_Number_Nodes();
}
//#####################################################################
// Function Mark_Nodes_Referenced
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const
{
    mesh.Mark_Nodes_Referenced(marks,mark);
}
//#####################################################################
// Function Volume
//#####################################################################
struct UNUSABLE{};
template<class T_SIMPLICIAL_OBJECT> typename T_SIMPLICIAL_OBJECT::SCALAR
Filled_Volume_Helper(const T_SIMPLICIAL_OBJECT& object, enable_if_t<!(T_SIMPLICIAL_OBJECT::MESH::dimension==(T_SIMPLICIAL_OBJECT::VECTOR_T::m-1)),UNUSABLE> unused=UNUSABLE())
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class T_SIMPLICIAL_OBJECT> typename T_SIMPLICIAL_OBJECT::SCALAR
Filled_Volume_Helper(const T_SIMPLICIAL_OBJECT& object, enable_if_t<(T_SIMPLICIAL_OBJECT::MESH::dimension==(T_SIMPLICIAL_OBJECT::VECTOR_T::m-1)),UNUSABLE> unused=UNUSABLE())
{
    typedef typename T_SIMPLICIAL_OBJECT::SCALAR T;typedef typename T_SIMPLICIAL_OBJECT::VECTOR_T TV;//enum WORKAROUND{d=T_SIMPLICIAL_OBJECT::MESH::dimension};
    static const int d=T_SIMPLICIAL_OBJECT::MESH::dimension;
    if(d!=TV::m-1) PHYSBAM_FATAL_ERROR("only codimension 1 objects can be filled");
    const TV base=object.particles.X(object.mesh.elements(0)[0]);
    T scaled_volume=0; // (d+1)!*volume
    for(int t=0;t<object.mesh.elements.m;t++){const VECTOR<int,d+1>& nodes=object.mesh.elements(t);
        MATRIX<T,TV::m,d+1> DX;for(int i=0;i<nodes.m;i++) DX.Set_Column(i,object.particles.X(nodes[i])-base);
        scaled_volume+=DX.Parallelepiped_Measure();}
    return (T)1/factorial(d+1)*scaled_volume;
}

template<class T> T
Filled_Volume_Helper(const MESH_OBJECT<VECTOR<T,1>, POINT_SIMPLEX_MESH>& object)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class TV,class T_MESH>
typename TV::SCALAR MESH_OBJECT<TV,T_MESH>::Volumetric_Volume()
{
    if(!mesh.elements.m) PHYSBAM_FATAL_ERROR("mesh has no elements");
    return Filled_Volume_Helper(*this);
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Read(TYPED_ISTREAM input)
{
    int size;
    Read_Binary(input,mesh,size);
    particles.Clean_Memory(); // strip everything away except for position
    particles.Resize(size);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,particles.X.Get_Array_Pointer(),size);
    else Read_Binary_Array<float>(input.stream,particles.X.Get_Array_Pointer(),size);
    if(mesh.elements.m){
        int min_index=mesh.elements.Flattened().Min(),max_index=mesh.elements.Flattened().Max();
        if(min_index<0) throw READ_ERROR(LOG::sprintf("Invalid vertex index %d",min_index));
        if(max_index>=particles.Size()) throw READ_ERROR(LOG::sprintf("Read invalid vertex index %d (particles.Size() = %d)",max_index,particles.Size()));
        Update_Number_Nodes();}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Write(TYPED_OSTREAM output) const
{
    if(mesh.number_nodes!=particles.Size()) PHYSBAM_FATAL_ERROR("number_nodes mismatch");
    Write_Binary(output,mesh,particles.X);
}
//#####################################################################
// Function Read_Obj
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Read_Obj(const std::string& filename)
{
    Clean_Memory();
    std::ifstream fin(filename);
    char buffer[2048];
    double d;
    int n;
    while(fin){
        fin.getline(buffer,2048);
        if(buffer[0]==0 || buffer[0]=='#' || buffer[1]!=' ') continue;
        if(buffer[0]=='v'){
            TV X;
            char* s=buffer+2;
            for(int i=0;i<TV::m;i++){
                sscanf(s,"%lf%n",&d,&n);
                X(i)=d;
                s+=n;}
            int p=particles.Add_Element();
            particles.X(p)=X;}
        else if(buffer[0]=='f'){
            typename T_MESH::ELEMENT_TYPE E;
            char* s=buffer+2;
            for(int i=0;i<E.m-1;i++){
                sscanf(s,"%d%n",&E(i),&n);
                s+=n;
                if(*s=='/') s+=strcspn(s," \t\r\n");}
            while(sscanf(s,"%d%n",&E(E.m-1),&n)>0){
                mesh.elements.Append(E-1);
                s+=n;
                if(*s=='/') s+=strcspn(s," \t\r\n");
                E(E.m-2)=E(E.m-1);}}}
    Update_Number_Nodes();
}
//#####################################################################
// Function Write_Obj
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Write_Obj(const std::string& filename) const
{
    const char* coord_str_comma[4]={"","x","x,y","x,y,z"};
    const char* coord_str_space[4]={"","x","x y","x y z"};
    const char* elem_str_comma[9]={"","a","a,b","a,b,c","a,b,c,d",0,0,0,"a,b,c,d,e,f,g,h"};
    const char* elem_str_space[9]={"","a","a b","a b c","a b c d",0,0,0,"a b c d e f g h"};
    std::ofstream fout(filename);
    fout<<"# simple obj file format:\n";
    fout<<"#  vertex at coordinates ("<<coord_str_comma[TV::m]<<")\n";
    fout<<"#   v "<<coord_str_space[TV::m]<<"\n";
    fout<<"#  element with vertices "<<elem_str_comma[T_MESH::ELEMENT_TYPE::m]<<"\n";
    fout<<"#   f "<<elem_str_space[T_MESH::ELEMENT_TYPE::m]<<"\n";
    fout<<"# vertices are indexed starting from 1.\n";

    for(int i=0;i<particles.X.m;i++){
        fout<<"\nv ";
        particles.X(i).Write_Raw(fout);}

    for(int i=0;i<mesh.elements.m;i++){
        fout<<"\nf ";
        (mesh.elements(i)+1).Write_Raw(fout);}

    fout<<"\n";
}
//#####################################################################
// Function Compact_Copy
//#####################################################################
template<class TV,class T_MESH> void MESH_OBJECT<TV,T_MESH>::
Compact_Copy(const MESH_OBJECT& obj,ARRAY<int>* mapping)
{
    ARRAY<int> particle_map(obj.particles.X.m);
    particle_map.Subset(obj.mesh.elements.Flattened()).Fill(1);
    int next_p=0;
    for(int i=0;i<particle_map.m;i++)
    {
        if(particle_map(i)) particle_map(i)=next_p++;
        else particle_map(i)=-1;
    }
    mesh.elements.Resize(obj.mesh.elements.m);
    mesh.elements.Flattened()=particle_map.Subset(obj.mesh.elements.Flattened());

    particles.Add_Elements(next_p);
    for(int i=0;i<particle_map.m;i++)
        if(particle_map(i)>=0)
            particles.X(particle_map(i))=obj.particles.X(i);
    Update_Number_Nodes();
    if(mapping) mapping->Exchange(particle_map);
}
//#####################################################################
namespace PhysBAM{
template class MESH_OBJECT<VECTOR<float,1>,POINT_SIMPLEX_MESH>;
template class MESH_OBJECT<VECTOR<float,1>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<float,2>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<float,3>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<float,2>,TRIANGLE_MESH>;
template class MESH_OBJECT<VECTOR<float,3>,TRIANGLE_MESH>;
template class MESH_OBJECT<VECTOR<float,3>,TETRAHEDRON_MESH>;
template class MESH_OBJECT<VECTOR<float,3>,HEXAHEDRON_MESH>;
template class MESH_OBJECT<VECTOR<double,1>,POINT_SIMPLEX_MESH>;
template class MESH_OBJECT<VECTOR<double,1>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<double,2>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<double,3>,SEGMENT_MESH>;
template class MESH_OBJECT<VECTOR<double,2>,TRIANGLE_MESH>;
template class MESH_OBJECT<VECTOR<double,3>,TRIANGLE_MESH>;
template class MESH_OBJECT<VECTOR<double,3>,TETRAHEDRON_MESH>;
template class MESH_OBJECT<VECTOR<double,3>,HEXAHEDRON_MESH>;
}
