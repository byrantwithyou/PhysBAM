//#####################################################################
// Copyright 2005, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
using namespace PhysBAM;
//#################################################################
// Function Compute_Material_Tets_In_Embedding_Tets
//#################################################################
template <class T> void
Compute_Material_Tets_In_Embedding_Tets(const TETRAHEDRALIZED_VOLUME<T>& material_tetrahedralized_volume,const TETRAHEDRALIZED_VOLUME<T>& embedding_tetrahedralized_volume,
                                        ARRAY<ARRAY<int> >& material_tets_in_embedding_tets,const T thickness_over_two=0,const bool verbose=true)
{
    ARRAY<TETRAHEDRON<T> >& material_tetrahedron_list=*material_tetrahedralized_volume.tetrahedron_list;
    TETRAHEDRON_HIERARCHY<T>& material_tetrahedron_hierarchy=*material_tetrahedralized_volume.tetrahedron_hierarchy;
    ARRAY<TETRAHEDRON<T> >& embedding_tetrahedron_list=*embedding_tetrahedralized_volume.tetrahedron_list;

    material_tets_in_embedding_tets.Exact_Resize(embedding_tetrahedron_list.m);
    for(int et=1;et<=embedding_tetrahedron_list.m;et++){
        TETRAHEDRON<T> embedding_tet=embedding_tetrahedron_list(et).Thickened(thickness_over_two);
        ARRAY<int> intersection_list;
        material_tetrahedron_hierarchy.Intersection_List(embedding_tet.Bounding_Box(),intersection_list);
        for(int i=1;i<=intersection_list.m;i++){
            TETRAHEDRON<T>& material_tet=material_tetrahedron_list(intersection_list(i));
            if(SIMPLEX_INTERACTIONS<T>::Tetrahedron_Tetrahedron_Intersection(embedding_tet.x1,embedding_tet.x2,embedding_tet.x3,embedding_tet.x4,material_tet.x1,material_tet.x2,material_tet.x3,material_tet.x4))
                material_tets_in_embedding_tets(et).Append(intersection_list(i));}
        if(verbose&&et%10000==0)std::cout<<"Processed "<<et<<" tets"<<std::endl;}
}
//#################################################################
// Function Extend_Material_Tets_In_Embedding_Tets
//#################################################################
template <class T> void
Extend_Material_Tets_In_Embedding_Tets(const TETRAHEDRALIZED_VOLUME<T>& material_tetrahedralized_volume,ARRAY<ARRAY<int> >& material_tets_in_embedding_tets,const bool verbose=true)
{
    TETRAHEDRON_MESH& tetrahedron_mesh=material_tetrahedralized_volume.tetrahedron_mesh;
    ARRAY<ARRAY<int> >& neighbor_tetrahedrons=*tetrahedron_mesh.neighbor_tetrahedrons;
    OPERATION_HASH tet_present;tet_present.Initialize(tetrahedron_mesh.tetrahedrons.m);

    for(int et=1;et<=material_tets_in_embedding_tets.m;et++){
        int initial_material_tet_intersections=material_tets_in_embedding_tets(et).m;
        for(int i=1;i<=initial_material_tet_intersections;i++)tet_present.Mark(material_tets_in_embedding_tets(et)(i));
        for(int i=1;i<=initial_material_tet_intersections;i++){
            int mt=material_tets_in_embedding_tets(et)(i);
            for(int j=1;j<=neighbor_tetrahedrons(mt).m;j++){
                int nmt=neighbor_tetrahedrons(mt)(j);
                if(!tet_present.Is_Marked_Current(nmt)){material_tets_in_embedding_tets(et).Append(nmt);tet_present.Mark(nmt);}}}
        tet_present.Next_Operation();if(verbose&&et%10000==0)std::cout<<"Processed "<<et<<" tets"<<std::endl;}
}
//#################################################################
// Function Propagate_Material_To_Empty_Embedding_Tets
//#################################################################
template <class T> void
Propagate_Material_To_Empty_Embedding_Tets(const TETRAHEDRALIZED_VOLUME<T>& embedding_tetrahedralized_volume,ARRAY<ARRAY<int> >& material_tets_in_embedding_tets,const bool verbose=true)
{
    TETRAHEDRON_MESH& tetrahedron_mesh=embedding_tetrahedralized_volume.tetrahedron_mesh;
    ARRAY<ARRAY<int> >& adjacent_tetrahedrons=*embedding_tetrahedralized_volume.tetrahedron_mesh.adjacent_tetrahedrons;
    OPERATION_HASH tet_occupied;tet_occupied.Initialize(tetrahedron_mesh.tetrahedrons.m);

    for(int sweep=1;;sweep++){
        ARRAY<int> empty_tets;
        for(int et=1;et<=material_tets_in_embedding_tets.m;et++)if(material_tets_in_embedding_tets(et).m)tet_occupied.Mark(et);else empty_tets.Append(et);
        if(verbose)std::cout<<"Sweep "<<sweep<<" : "<<empty_tets.m<<" empty tets found"<<std::endl;
        if(!empty_tets.m)break;
        for(int i=1;i<=empty_tets.m;i++){
            int et=empty_tets(i);
            for(int j=1;j<=adjacent_tetrahedrons(et).m;j++){
                int aet=adjacent_tetrahedrons(et)(j);
                if(tet_occupied.Is_Marked_Current(aet))material_tets_in_embedding_tets(et).Append_Unique_Elements(material_tets_in_embedding_tets(aet));}}}
}
//#####################################################################
// Function Compute_Connected_Components_And_Split
//#####################################################################
template <class T>void 
Compute_Connected_Components_And_Split(const TETRAHEDRALIZED_VOLUME<T>& material_tetrahedralized_volume,ARRAY<ARRAY<int> >& material_tets_in_embedding_tets,
                                       ARRAY<ARRAY<int> >& material_tets_in_split_embedding_tets,ARRAY<ARRAY<int> >& split_embedding_tets_in_parent_embedding_tets,const bool verbose=true)
{
    TETRAHEDRON_MESH& tetrahedron_mesh=material_tetrahedralized_volume.tetrahedron_mesh;
    ARRAY<ARRAY<int> >& neighbor_tetrahedrons=*tetrahedron_mesh.neighbor_tetrahedrons;
    OPERATION_HASH tet_present;tet_present.Initialize(tetrahedron_mesh.tetrahedrons.m);
    int original_embedding_tets=material_tets_in_embedding_tets.m,split_embedding_tets=0;

    split_embedding_tets_in_parent_embedding_tets.Resize(original_embedding_tets);material_tets_in_split_embedding_tets.Resize(0);
    for(int et=1;et<=original_embedding_tets;et++){
        ARRAY<int> unassigned_material_tets=material_tets_in_embedding_tets(et);
        for(int i=1;i<=unassigned_material_tets.m;i++)tet_present.Mark(unassigned_material_tets(i));
        while(unassigned_material_tets.m){
            // Create new connected component
            split_embedding_tets_in_parent_embedding_tets(et).Append(++split_embedding_tets);
            material_tets_in_split_embedding_tets.Resize(split_embedding_tets);
            ARRAY<int>& current_component=material_tets_in_split_embedding_tets(split_embedding_tets);
            // Initialize with any unassigned tet
            current_component.Append(unassigned_material_tets(1));
            tet_present.Unmark(unassigned_material_tets(1));
            unassigned_material_tets.Remove_Index_Lazy(1);
            // Sweep component, adding neighbors as necessary
            for(int i=1;i<=current_component.m;i++){
                int mt=current_component(i);
                for(int j=1;j<=neighbor_tetrahedrons(mt).m;j++){
                    int nmt=neighbor_tetrahedrons(mt)(j);
                    if(tet_present.Is_Marked_Current(nmt)){
                        current_component.Append(nmt);
                        tet_present.Unmark(nmt);
                        int k=unassigned_material_tets.Find(nmt);if(!k){std::cerr<<"Mesh inconsistency"<<std::endl;exit(1);}
                        unassigned_material_tets.Remove_Index_Lazy(k);}}}}
        if(verbose&&et%10000==0)std::cout<<"Processed "<<et<<" tets"<<std::endl;}
}
//#####################################################################
// Function Merge_Contiguous_Split_Embedding_Tetrahedra
//#####################################################################
template <class T> void
Merge_Contiguous_Split_Embedding_Tetrahedra(const TETRAHEDRALIZED_VOLUME<T>& material_tetrahedralized_volume,
                                            TETRAHEDRALIZED_VOLUME<T>& embedding_tetrahedralized_volume,
                                            ARRAY<ARRAY<int> >& material_tets_in_split_embedding_tets,
                                            ARRAY<ARRAY<int> >& split_embedding_tets_in_parent_embedding_tets,
                                            const bool verbose=true)
{
    TETRAHEDRON_MESH &embedding_tetrahedron_mesh=embedding_tetrahedralized_volume.tetrahedron_mesh;
    TETRAHEDRON_MESH &material_tetrahedron_mesh=material_tetrahedralized_volume.tetrahedron_mesh;
    SOLIDS_PARTICLES<T,VECTOR_3D<T> > &embedding_particles=embedding_tetrahedralized_volume.particles;
    ARRAYS<int> &embedding_triangle_tetrahedrons=*embedding_tetrahedron_mesh.triangle_tetrahedrons;
    ARRAY<bool> tet_present(material_tetrahedron_mesh.tetrahedrons.m);
    ARRAY<int> new_particle_indices(4*material_tets_in_split_embedding_tets.m),parent_tets(material_tets_in_split_embedding_tets.m);
    ARRAY<ARRAY<int> > geometrically_distinct_nodes_to_duplicates(embedding_particles.number);
    SOLIDS_PARTICLES<T,VECTOR_3D<T> > new_particles;
    UNION_FIND vnodes;
    
    // Initialize lookup structures
    vnodes.Initialize(4*material_tets_in_split_embedding_tets.m);
    ARRAY<int>::copy(0,vnodes.parents);
    ARRAY<bool>::copy(false,tet_present);
    // Check for faces that should be identified
    for (int tri=1;tri<=embedding_triangle_tetrahedrons.m;tri++){
        int tet1=embedding_triangle_tetrahedrons(1,tri),tet2=embedding_triangle_tetrahedrons(2,tri);
        if (tet1==0||tet2==0) continue;
        // Test all geometrically adjacent tets
        for (int i=1;i<=split_embedding_tets_in_parent_embedding_tets(tet1).m;i++) for (int j=1;j<=split_embedding_tets_in_parent_embedding_tets(tet2).m;j++){
            int dup1=split_embedding_tets_in_parent_embedding_tets(tet1)(i),dup2=split_embedding_tets_in_parent_embedding_tets(tet2)(j);
            bool tets_contiguous=false;
            // Weaker version of material continuity used for robustness (for analytical test the common tet should intersect the common face of the embedding tets)
            for (int k=1;k<=material_tets_in_split_embedding_tets(dup1).m;k++) tet_present(material_tets_in_split_embedding_tets(dup1)(k))=true;
            for (int k=1;k<=material_tets_in_split_embedding_tets(dup2).m;k++) if (tet_present(material_tets_in_split_embedding_tets(dup2)(k))) tets_contiguous=true;
            for (int k=1;k<=material_tets_in_split_embedding_tets(dup1).m;k++) tet_present(material_tets_in_split_embedding_tets(dup1)(k))=false;
            if (tets_contiguous) for (int k=1;k<=4;k++) for (int l=1;l<=4;l++) 
                if (embedding_tetrahedron_mesh.tetrahedrons(k,tet1)==embedding_tetrahedron_mesh.tetrahedrons(l,tet2)) vnodes.Union(4*dup1+k-4,4*dup2+l-4);}}
    // Collapse equivalent particles
    int new_particle_count=0; 
    for (int i=1;i<=vnodes.parents.m;i++) if (vnodes.Is_Root(i)) new_particle_indices(i)=++new_particle_count;
    new_particles.Add_Particles(new_particle_count);
    for (int i=1;i<=split_embedding_tets_in_parent_embedding_tets.m;i++) for (int j=1;j<=split_embedding_tets_in_parent_embedding_tets(i).m;j++) for (int k=1;k<=4;k++) {
        int node_index=4*split_embedding_tets_in_parent_embedding_tets(i)(j)+k-4;
        if (!vnodes.Is_Root(node_index)) continue;
        new_particles.X(new_particle_indices(node_index))=embedding_particles.X(embedding_tetrahedron_mesh.tetrahedrons(k,i));
        geometrically_distinct_nodes_to_duplicates(embedding_tetrahedron_mesh.tetrahedrons(k,i)).Append_Unique(new_particle_indices(node_index));}
    // (For debugging purposes) display node multiplicity stats
    int max_node_multiplicity=0;
    for (int i=1;i<=geometrically_distinct_nodes_to_duplicates.m;i++) 
        if (max_node_multiplicity<geometrically_distinct_nodes_to_duplicates(i).m) max_node_multiplicity=geometrically_distinct_nodes_to_duplicates(i).m;
    for (int i=0;i<=max_node_multiplicity;i++){
        int count=0;for (int j=1;j<=geometrically_distinct_nodes_to_duplicates.m;j++) if (geometrically_distinct_nodes_to_duplicates(j).m==i) count++;
        std::cout<<"Nodes with multiplicity "<<i<<" : "<<count<<std::endl;}
    // Delete acceleration structures
    embedding_tetrahedron_mesh.Delete_Auxiliary_Structures();
    // Replace old tetrahedron mesh
    embedding_tetrahedron_mesh.tetrahedrons.Resize(4,material_tets_in_split_embedding_tets.m);
    for (int i=1;i<=material_tets_in_split_embedding_tets.m;i++) for (int j=1;j<=4;j++) embedding_tetrahedron_mesh.tetrahedrons(j,i)=new_particle_indices(vnodes.Find(4*i+j-4));
    // Copy particles array
    embedding_particles.Initialize_Particles(new_particles);embedding_tetrahedron_mesh.number_nodes=embedding_particles.number;
    // Remove duplicate tets (update split_embedding_tets_in_parent_embedding_tets in the process)
    for (int i=1;i<=split_embedding_tets_in_parent_embedding_tets.m;i++) for (int j=1;j<=split_embedding_tets_in_parent_embedding_tets(i).m;j++)
        parent_tets(split_embedding_tets_in_parent_embedding_tets(i)(j))=i;
    for (int i=1;i<=parent_tets.m;i++) if (!parent_tets(i)) {std::cerr<<"Dangling duplicate"<<std::endl;exit(1);}
    for (int i=1;i<=split_embedding_tets_in_parent_embedding_tets.m;i++) 
        for (int j=1;j<split_embedding_tets_in_parent_embedding_tets(i).m;j++) for (int k=j+1;k<=split_embedding_tets_in_parent_embedding_tets(i).m;k++)
            if (embedding_tetrahedron_mesh.tetrahedrons(1,j)==embedding_tetrahedron_mesh.tetrahedrons(1,k)&&
                embedding_tetrahedron_mesh.tetrahedrons(2,j)==embedding_tetrahedron_mesh.tetrahedrons(2,k)&&
                embedding_tetrahedron_mesh.tetrahedrons(3,j)==embedding_tetrahedron_mesh.tetrahedrons(3,k)&&
                embedding_tetrahedron_mesh.tetrahedrons(4,j)==embedding_tetrahedron_mesh.tetrahedrons(4,k)){
                parent_tets(split_embedding_tets_in_parent_embedding_tets(i)(k))=0;
                material_tets_in_split_embedding_tets(j).Append_Elements(material_tets_in_split_embedding_tets(k));
                std::cout<<"Identical tets found. Discarding duplicate ..."<<std::endl;}
    for (int i=1;i<=parent_tets.m;i++) if (!parent_tets(i)) {parent_tets.Remove_Index(i);embedding_tetrahedron_mesh.tetrahedrons.Remove_Index(i);}
    for (int i=1;i<=split_embedding_tets_in_parent_embedding_tets.m;i++) split_embedding_tets_in_parent_embedding_tets(i).Clean_Memory();
    for (int i=1;i<=parent_tets.m;i++) split_embedding_tets_in_parent_embedding_tets(parent_tets(i)).Append(i);
}
//#################################################################
// Function Compute_Nonmanifold_Embedding
//#################################################################
void Print_Histogram(const ARRAY<ARRAY<int> >& data)
{
    int max_size=0;
    for(int i=1;i<=data.m;i++)if(data(i).m>max_size)max_size=data(i).m;
    ARRAYS<VECTOR<int,1> > histogram(0,max_size);
    for(int i=0;i<=max_size;i++)histogram(i)=0;
    for(int i=1;i<=data.m;i++)histogram(data(i).m)++;
    for(int i=0;i<=max_size;i++)std::cout<<i<<" : "<<histogram(i)<<std::endl;
}

template <class T> void
Compute_Nonmanifold_Embedding(TETRAHEDRALIZED_VOLUME<T>& material_tetrahedralized_volume,TETRAHEDRALIZED_VOLUME<T>& embedding_tetrahedralized_volume)
{
    ARRAY<ARRAY<int> > material_tets_in_embedding_tets;
    ARRAY<ARRAY<int> > material_tets_in_split_embedding_tets;
    ARRAY<ARRAY<int> > split_embedding_tets_in_parent_embedding_tets;

    const T tetrahedron_intersection_thickness=0;

    // Initialize acceleration structures
    std::cout<<"Initializing acceleration structures"<<std::endl;
    if(!material_tetrahedralized_volume.tetrahedron_list)material_tetrahedralized_volume.Update_Tetrahedron_List();
    if(!material_tetrahedralized_volume.tetrahedron_hierarchy)material_tetrahedralized_volume.Initialize_Tetrahedron_Hierarchy();
    if(!material_tetrahedralized_volume.tetrahedron_mesh.incident_tetrahedrons)material_tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();
    if(!material_tetrahedralized_volume.tetrahedron_mesh.neighbor_tetrahedrons)material_tetrahedralized_volume.tetrahedron_mesh.Initialize_Neighbor_Tetrahedrons();
    if(!embedding_tetrahedralized_volume.tetrahedron_list)embedding_tetrahedralized_volume.Update_Tetrahedron_List();
    if(!embedding_tetrahedralized_volume.tetrahedron_mesh.triangle_mesh)embedding_tetrahedralized_volume.tetrahedron_mesh.Initialize_Triangle_Mesh();
    if(!embedding_tetrahedralized_volume.tetrahedron_mesh.incident_tetrahedrons)embedding_tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();
    if(!embedding_tetrahedralized_volume.tetrahedron_mesh.adjacent_tetrahedrons)embedding_tetrahedralized_volume.tetrahedron_mesh.Initialize_Adjacent_Tetrahedrons();
    if(!embedding_tetrahedralized_volume.tetrahedron_mesh.triangle_tetrahedrons)embedding_tetrahedralized_volume.tetrahedron_mesh.Initialize_Triangle_Tetrahedrons();
    
    std::cout<<"Computing material tets intersecting embedding tets"<<std::endl;
    Compute_Material_Tets_In_Embedding_Tets(material_tetrahedralized_volume,embedding_tetrahedralized_volume,material_tets_in_embedding_tets,tetrahedron_intersection_thickness);
    std::cout<<"Extending material tet intersection lists to their 1-ring"<<std::endl;
    Extend_Material_Tets_In_Embedding_Tets(material_tetrahedralized_volume,material_tets_in_embedding_tets);
    std::cout<<"Propagating material content to empty embedding tets"<<std::endl;
    Propagate_Material_To_Empty_Embedding_Tets(embedding_tetrahedralized_volume,material_tets_in_embedding_tets);

    FILE_UTILITIES::Write_To_File<T>("material_tets_in_embedding_tets.dat",material_tets_in_embedding_tets);
    //FILE_UTILITIES::Read_From_File<T>("material_tets_in_embedding_tets.dat",material_tets_in_embedding_tets);
    //std::cout<<"material_tets_in_embedding_tets"<<std::endl;Print_Histogram(material_tets_in_embedding_tets);

    std::cout<<"Computing connected components and splitting disconnected elements"<<std::endl;
    Compute_Connected_Components_And_Split(material_tetrahedralized_volume,material_tets_in_embedding_tets,material_tets_in_split_embedding_tets,split_embedding_tets_in_parent_embedding_tets);
    FILE_UTILITIES::Write_To_File<T>("material_tets_in_split_embedding_tets.dat",material_tets_in_split_embedding_tets);
    FILE_UTILITIES::Write_To_File<T>("split_embedding_tets_in_parent_embedding_tets.dat",split_embedding_tets_in_parent_embedding_tets);
    //FILE_UTILITIES::Read_From_File<T>("material_tets_in_split_embedding_tets.dat",material_tets_in_split_embedding_tets);
    //FILE_UTILITIES::Read_From_File<T>("split_embedding_tets_in_parent_embedding_tets.dat",split_embedding_tets_in_parent_embedding_tets);
    //std::cout<<"material_tets_in_split_embedding_tets"<<std::endl;Print_Histogram(material_tets_in_split_embedding_tets);
    std::cout<<"split_embedding_tets_in_parent_embedding_tets"<<std::endl;Print_Histogram(split_embedding_tets_in_parent_embedding_tets);

    std::cout<<"Merging contiguous tetrahedra"<<std::endl;
    Merge_Contiguous_Split_Embedding_Tetrahedra(material_tetrahedralized_volume,embedding_tetrahedralized_volume,material_tets_in_split_embedding_tets,split_embedding_tets_in_parent_embedding_tets);
}
//#################################################################
// Function main
//#################################################################
template <class T,class RW>
int main_templatized(int argc,char *argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-m","","material geometry file","material geometry file");
    parse_args.Add_String_Argument("-e","","embedding geometry file","embedding geometry file");
    parse_args.Add_String_Argument("-o","","output mesh file","output mesh file");
    parse_args.Parse(argc,argv);

    std::string material_filename,embedding_filename,output_filename;
    if(parse_args.Is_Value_Set("-m")) material_filename=parse_args.Get_String_Value("-m");else parse_args.Print_Usage(true);
    if(parse_args.Is_Value_Set("-e")) embedding_filename=parse_args.Get_String_Value("-e");else parse_args.Print_Usage(true);
    if(parse_args.Is_Value_Set("-o")) output_filename=parse_args.Get_String_Value("-o");else parse_args.Print_Usage(true);

    TETRAHEDRALIZED_VOLUME<T>* material_tetrahedralized_volume;FILE_UTILITIES::Create_From_File<RW>(material_filename,material_tetrahedralized_volume);
    TETRAHEDRALIZED_VOLUME<T>* embedding_tetrahedralized_volume;FILE_UTILITIES::Create_From_File<RW>(embedding_filename,embedding_tetrahedralized_volume);
    Compute_Nonmanifold_Embedding(*material_tetrahedralized_volume,*embedding_tetrahedralized_volume);
    FILE_UTILITIES::Write_To_File<RW>(output_filename,*embedding_tetrahedralized_volume);

    return 0;
}
int main(int argc,char *argv[]){
    return main_templatized<float,float>(argc,argv);
}
//#################################################################
