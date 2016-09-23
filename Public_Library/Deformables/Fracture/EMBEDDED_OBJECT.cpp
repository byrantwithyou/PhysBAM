//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_OBJECT
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/permutation.h>
#include <Core/Math_Tools/RANGE.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Fracture/EMBEDDED_OBJECT.h>
#include <Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Fracture/EMBEDDED_TRIANGULATED_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> EMBEDDED_OBJECT<TV,d>::
EMBEDDED_OBJECT(T_SIMPLICIAL_OBJECT& simplicial_object_input)
    :particles(simplicial_object_input.particles),embedded_particles(particles),interpolation_fraction_threshold((T).1),
    embedded_children_index(0),embedded_children(0),parents_to_embedded_particles_hash_table(0),hashtable_multiplier(12),
    embedded_subelements_in_parent_element_index(0),number_of_embedded_subelements_in_parent_element(0),
    embedded_subelements_in_parent_element(0),average_interpolation_fractions(false),
    simplicial_object(simplicial_object_input),embedded_object(embedded_mesh,particles),bounding_box(0),need_destroy_simplicial_object(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> EMBEDDED_OBJECT<TV,d>::
~EMBEDDED_OBJECT()
{
    if(need_destroy_simplicial_object) delete &simplicial_object;
    delete embedded_children_index;delete embedded_children;delete parents_to_embedded_particles_hash_table;delete embedded_subelements_in_parent_element_index;
    delete number_of_embedded_subelements_in_parent_element;delete embedded_subelements_in_parent_element;delete bounding_box;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Clean_Memory()
{
    embedded_particles.active_indices.Clean_Memory();parent_particles.Clean_Memory();interpolation_fraction.Clean_Memory();
    embedded_mesh.Clean_Memory();node_in_simplex_is_material.Clean_Memory();Delete_Auxiliary_Structures();
}
//#####################################################################
// Function Delete_Auxiliary_Structures
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Delete_Auxiliary_Structures()
{
    delete embedded_children_index;embedded_children_index=0;delete embedded_children;embedded_children=0;
    delete parents_to_embedded_particles_hash_table;parents_to_embedded_particles_hash_table=0;
    delete embedded_subelements_in_parent_element_index;embedded_subelements_in_parent_element_index=0;
    delete number_of_embedded_subelements_in_parent_element;number_of_embedded_subelements_in_parent_element=0;
    delete embedded_subelements_in_parent_element;embedded_subelements_in_parent_element=0;
    embedded_object.Clean_Memory();delete bounding_box;bounding_box=0;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT* EMBEDDED_OBJECT<TV,d>::
Create()
{
    T_EMBEDDED_OBJECT* object=new T_EMBEDDED_OBJECT(*T_SIMPLICIAL_OBJECT::Create());
    object->need_destroy_simplicial_object=true;return object;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT* EMBEDDED_OBJECT<TV,d>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    T_EMBEDDED_OBJECT* object=new T_EMBEDDED_OBJECT(*T_SIMPLICIAL_OBJECT::Create(particles));
    object->Initialize_Parents_To_Embedded_Particles_Hash_Table(15); // TODO: reconsider
    object->Set_Interpolation_Fraction_Threshold((T)1e-4); // TODO: reconsider
    object->need_destroy_simplicial_object=true;return object;
}
//#####################################################################
// Function Reset_Parents_To_Embedded_Particles_Hash_Table
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Reset_Parents_To_Embedded_Particles_Hash_Table()
{
    if(parents_to_embedded_particles_hash_table) parents_to_embedded_particles_hash_table->Remove_All();
}
//#####################################################################
// Function Copy_Then_Reset_Parents_To_Embedded_Particles_Hash_Table
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Copy_Then_Reset_Parents_To_Embedded_Particles_Hash_Table(HASHTABLE<VECTOR<int,2>,int>*& hash_table_copy)
{
    hash_table_copy=parents_to_embedded_particles_hash_table;parents_to_embedded_particles_hash_table=new HASHTABLE<VECTOR<int,2>,int>(hashtable_multiplier*particles.Size()); // TODO: Rethink the size
}
//#####################################################################
// Function Initialize_Parents_To_Embedded_Particles_Hash_Table
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Initialize_Parents_To_Embedded_Particles_Hash_Table(const int new_hashtable_multiplier)
{
    if(new_hashtable_multiplier) hashtable_multiplier=new_hashtable_multiplier;
    if(parents_to_embedded_particles_hash_table && parents_to_embedded_particles_hash_table->Max_Size() >= hashtable_multiplier*particles.Size()) return;
    delete parents_to_embedded_particles_hash_table;parents_to_embedded_particles_hash_table=new HASHTABLE<VECTOR<int,2>,int>(hashtable_multiplier*particles.Size());
    for(int p=0;p<embedded_particles.active_indices.m;p++)
        parents_to_embedded_particles_hash_table->Insert(parent_particles(p).Sorted(),p);
}
//#####################################################################
// Function Update_Embedded_Particle_Positions
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Update_Embedded_Particle_Positions()
{
    for(int p=0;p<embedded_particles.active_indices.m;p++){
        int i,j;parent_particles(p).Get(i,j);
        embedded_particles.X(p)=LINEAR_INTERPOLATION<typename TV::SCALAR,TV>::Linear(particles.X(i),particles.X(j),interpolation_fraction(p));}
}
//#####################################################################
// Function Initialize_Embedded_Children
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Initialize_Embedded_Children()
{
    delete embedded_children_index;delete embedded_children;
    embedded_children_index=new ARRAY<int>(particles.Size());embedded_children=new ARRAY<ARRAY<int> >();
    for(int p=0;p<embedded_particles.active_indices.m;p++) Add_Embedded_Particle_To_Embedded_Children(p);
}
//#####################################################################
// Function Add_Embedded_Particle_To_Embedded_Children
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Particle_To_Embedded_Children(const int embedded_particle)
{
    int parent1,parent2;parent_particles(embedded_particle).Get(parent1,parent2);
    int index1=(*embedded_children_index)(parent1),index2=(*embedded_children_index)(parent2);
    if(index1) (*embedded_children)(index1).Append(embedded_particle);
    else{
        embedded_children->Append(ARRAY<int>());(*embedded_children_index)(parent1)=embedded_children->m;
        (*embedded_children)(embedded_children->m).Append(embedded_particle);}
    if(index2) (*embedded_children)(index2).Append(embedded_particle);
    else{
        embedded_children->Append(ARRAY<int>());(*embedded_children_index)(parent2)=embedded_children->m;
        (*embedded_children)(embedded_children->m).Append(embedded_particle);}
}
//#####################################################################
// Function Add_Emdedded_Particle_If_Not_Already_There
//#####################################################################
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Particle_If_Not_Already_There(const VECTOR<int,2>& nodes,const T interpolation_fraction_input)
{
    int current_embedded_particle=Embedded_Particle_On_Segment(nodes);
    if(current_embedded_particle>=0) return current_embedded_particle;
    return Add_Embedded_Particle(nodes,interpolation_fraction_input);
}
//#####################################################################
// Function Add_Emdedded_Particle
//#####################################################################
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Particle(const VECTOR<int,2>& nodes,const T interpolation_fraction_input,const bool reinitialize_hash_table)
{
    assert(simplicial_object.mesh.number_nodes==particles.Size() && embedded_mesh.number_nodes==particles.Size() && Embedded_Particle_On_Segment(nodes)<0);
    int new_embedded_particle=embedded_particles.Add_Element();
    simplicial_object.mesh.Add_Nodes(1);embedded_mesh.Add_Nodes(1); // Need to update meshes and acceleration structures
    interpolation_fraction.Append(Clamp_Interpolation_Fraction(interpolation_fraction_input));
    parent_particles.Append(nodes);
    parents_to_embedded_particles_hash_table->Insert(nodes.Sorted(),new_embedded_particle);
    if(reinitialize_hash_table) Initialize_Parents_To_Embedded_Particles_Hash_Table();
    if(embedded_children_index) Add_Embedded_Particle_To_Embedded_Children(new_embedded_particle);
    T lambda=interpolation_fraction(new_embedded_particle);
    embedded_particles.X(new_embedded_particle)=LINEAR_INTERPOLATION<T,TV>::Linear(particles.X(nodes[0]),particles.X(nodes[1]),lambda);
    return new_embedded_particle;
}
//#####################################################################
// Function Element_Containing_Subelement
//#####################################################################
template<int d> static inline VECTOR<int,d+1> Compact_Parents(const VECTOR<VECTOR<int,2>,d>& parents)
{
    VECTOR<int,d+1> all_parents(parents[0]);int m=1;
    for(int i=1;i<d;i++)for(int j=0;j<2;j++)
        if(!all_parents.Contains(parents[i][j])) all_parents[++m]=parents[i][j];
    assert(m==d+1);
    return all_parents;
}
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Element_Containing_Subelement(const int embedded_subelement) const
{
    VECTOR<int,d> embedded_nodes(embedded_particles.subset_index_from_point_cloud_index.Subset(embedded_mesh.elements(embedded_subelement)));
    VECTOR<VECTOR<int,2>,d> parents(parent_particles.Subset(embedded_nodes));
    return simplicial_object.mesh.Simplex(Compact_Parents(parents));
}
//#####################################################################
// Function Add_Embedded_Subelement_To_Embedded_Subelements_In_Element
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Subelement_To_Embedded_Subelements_In_Element(const int subelement)
{
    int element=Element_Containing_Subelement(subelement);
    int index=(*embedded_subelements_in_parent_element_index)(element);
    if(index>=0){
        int embedded_subelement_number=++(*number_of_embedded_subelements_in_parent_element)(index);
        (*embedded_subelements_in_parent_element)(index)(embedded_subelement_number)=subelement;}
    else{
        VECTOR<int,EMBEDDED_OBJECT<TV,d>::max_subelements_per_element> subelements;subelements[0]=subelement;
        embedded_subelements_in_parent_element->Append(subelements);
        number_of_embedded_subelements_in_parent_element->Append(1);
        (*embedded_subelements_in_parent_element_index)(element)=embedded_subelements_in_parent_element->m;}
}
//#####################################################################
// Function Add_Embedded_Subelement_If_Not_Already_There
//#####################################################################
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Subelement_If_Not_Already_There(const VECTOR<int,d>& embedded_nodes)
{
    VECTOR<int,d> global_particles(embedded_particles.active_indices.Subset(embedded_nodes));
    int t=embedded_mesh.Simplex(global_particles);
    if(t>=0) return t;
    return Add_Embedded_Subelement(embedded_nodes);
}
//#####################################################################
// Function Add_Embedded_Subelement
//#####################################################################
// only updates embedded_mesh.incident_elements and embedded_subelements_in_parent_element
template<class TV,int d> int EMBEDDED_OBJECT<TV,d>::
Add_Embedded_Subelement(const VECTOR<int,d>& embedded_nodes)
{
    VECTOR<int,d> global_particles(embedded_particles.active_indices.Subset(embedded_nodes));
    assert(embedded_mesh.Simplex(global_particles)<0);
    int new_subelement=embedded_mesh.elements.Append(global_particles);
    if(embedded_mesh.incident_elements) // needs to be updated
        for(int i=0;i<global_particles.m;i++) (*embedded_mesh.incident_elements)(global_particles[i]).Append(new_subelement);
    if(embedded_subelements_in_parent_element) Add_Embedded_Subelement_To_Embedded_Subelements_In_Element(new_subelement);
    return new_subelement;
}
//#####################################################################
// Function Initialize_Embedded_Subelements_In_Parent_Element
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Initialize_Embedded_Subelements_In_Parent_Element()
{
    delete embedded_subelements_in_parent_element_index;delete number_of_embedded_subelements_in_parent_element;delete embedded_subelements_in_parent_element;
    embedded_subelements_in_parent_element_index=new ARRAY<int>(simplicial_object.mesh.elements.m);
    number_of_embedded_subelements_in_parent_element=new ARRAY<int>(0);
    embedded_subelements_in_parent_element=new ARRAY<VECTOR<int,EMBEDDED_OBJECT<TV,d>::max_subelements_per_element> >();
    for(int t=0;t<embedded_mesh.elements.m;t++) Add_Embedded_Subelement_To_Embedded_Subelements_In_Element(t);
}
template<class TV> static void Add_Levelset_Cuts(EMBEDDED_TRIANGULATED_OBJECT<TV>& embedded_object,const int positive_count,const VECTOR<int,3>& nodes)
{
    int i,j,k;nodes.Get(i,j,k);
    if(positive_count==1) embedded_object.Add_Embedded_Segment(i,k,j,k); // one segment
    else if(positive_count==2) embedded_object.Add_Embedded_Segment(i,j,i,k);
}
template<class T> static void Add_Levelset_Cuts(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object,const int positive_count,const VECTOR<int,4>& nodes)
{
    int i,j,k,l;nodes.Get(i,j,k,l);
    if(positive_count==1) embedded_object.Add_Embedded_Triangle(i,l,j,l,k,l); // one triangle
    else if(positive_count==2){ // two triangles
        int ik=embedded_object.Embedded_Particle_On_Segment(i,k),il=embedded_object.Embedded_Particle_On_Segment(i,l),
            jk=embedded_object.Embedded_Particle_On_Segment(j,k),jl=embedded_object.Embedded_Particle_On_Segment(j,l);
        embedded_object.Add_Embedded_Triangle(il,jl,jk);embedded_object.Add_Embedded_Triangle(jk,ik,il);}
    else if(positive_count==3) embedded_object.Add_Embedded_Triangle(i,j,i,k,i,l); // one triangle
}
//#####################################################################
// Function Calculate_Boundary_From_Levelset_On_Nodes
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Calculate_Boundary_From_Levelset_On_Nodes(ARRAY<T>& phi,const bool discard_elements_outside_levelset,const bool verbose)
{
    if(discard_elements_outside_levelset){ // TODO: Consider doing this before initializing the embedded object
        assert(embedded_particles.active_indices.m==0); // Discarding should be done before adding any embedded particles
        for(int t=simplicial_object.mesh.elements.m-1;t>=0;t--){
            VECTOR<T,d+1> phi_t(phi.Subset(simplicial_object.mesh.elements(t)));
            if(phi_t.Min()>0) simplicial_object.mesh.elements.Remove_Index_Lazy(t);
        }
        ARRAY<int> condensation_mapping;simplicial_object.Discard_Valence_Zero_Particles_And_Renumber(condensation_mapping);
        ARRAY<T> phi_new(particles.Size());for(int k=0;k<phi.m;k++) if(condensation_mapping(k)>=0) phi_new(condensation_mapping(k))=phi(k);
        phi.Exchange(phi_new);
    }

    if(embedded_particles.active_indices.m){
        LOG::cout<<"Calculate_Boundary_From_Levelset_On_Nodes cannot be called with existing embedded particles"<<std::endl;
        PHYSBAM_FATAL_ERROR();}

    Initialize_Parents_To_Embedded_Particles_Hash_Table();
    embedded_particles.Update_Subset_Index_From_Element_Index(); // TODO: move this somewhere else
    embedded_mesh.number_nodes=particles.Size();

    bool segment_mesh_defined=simplicial_object.mesh.segment_mesh!=0;if(!segment_mesh_defined) simplicial_object.mesh.Initialize_Segment_Mesh();
    bool embedded_incident_elements_defined=embedded_mesh.incident_elements!=0;if(!embedded_incident_elements_defined) embedded_mesh.Initialize_Incident_Elements();

    // calculate embedded particles
    for(int s=0;s<simplicial_object.mesh.segment_mesh->elements.m;s++){
        int n1,n2;simplicial_object.mesh.segment_mesh->elements(s).Get(n1,n2);T phi1=phi(n1),phi2=phi(n2);
        if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
            int inside,outside;if(phi1 <= 0){inside=n1;outside=n2;}else{inside=n2;outside=n1;exchange(phi1,phi2);}
            T theta=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
            Add_Embedded_Particle(VECTOR<int,2>(inside,outside),theta);}}

    // calculate embedded simplices
    for(int t=0;t<simplicial_object.mesh.elements.m;t++){
        VECTOR<int,d+1> nodes=simplicial_object.mesh.elements(t);
        {int i=0,j=nodes.m-1;while(i<j){if(phi(nodes[i])>0) exchange(nodes[i],nodes[j--]);else i++;} // move inside nodes before outside nodes
        if((nodes.m-j)&1){if(phi(nodes[1])<=0) exchange(nodes[0],nodes[1]);else exchange(nodes[d-1],nodes[d]);}} // one final swap to ensure an even permutation
        for(int i=0;i<nodes.m-1;i++) assert((phi(nodes[i])>0) <= (phi(nodes[i+1])>0));
        int positive_count=0;for(int i=0;i<nodes.m;i++) if(phi(nodes[i])>0) positive_count++;
        Add_Levelset_Cuts(dynamic_cast<typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT&>(*this),positive_count,nodes);}

    node_in_simplex_is_material.Resize(simplicial_object.mesh.elements.m);
    for(int t=0;t<simplicial_object.mesh.elements.m;t++){VECTOR<int,d+1>& element=simplicial_object.mesh.elements(t);
        for(int i=0;i<element.m;i++) node_in_simplex_is_material(t)(i)=phi(element[i])<=0;}

    if(!segment_mesh_defined){delete simplicial_object.mesh.segment_mesh;simplicial_object.mesh.segment_mesh=0;}
    if(!embedded_incident_elements_defined){delete embedded_mesh.incident_elements;embedded_mesh.incident_elements=0;}

    if(verbose){
        LOG::cout << "total particles = " << particles.Size() << std::endl;
        LOG::cout << "total elements = " << simplicial_object.mesh.elements.m << std::endl;
        LOG::cout << "total embedded subelements = " << embedded_mesh.elements.m << std::endl;}
}
//#####################################################################
// Function Read_Helper
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Read(TYPED_ISTREAM& input)
{
    Clean_Memory(); // reads entire simplicial_object instead of just simplicial_mesh
    int backward_compatible;
    Read_Binary(input,embedded_particles,backward_compatible,parent_particles,interpolation_fraction,simplicial_object,embedded_mesh,backward_compatible);
    Read_Binary(input,node_in_simplex_is_material,interpolation_fraction_threshold,orientation_index);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV,int d> void EMBEDDED_OBJECT<TV,d>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,embedded_particles,2,parent_particles,interpolation_fraction,simplicial_object,embedded_mesh,d+1,node_in_simplex_is_material);
    Write_Binary(output,interpolation_fraction_threshold,orientation_index);
}
//#####################################################################
namespace PhysBAM{
template class EMBEDDED_OBJECT<VECTOR<float,2>,2>;
template class EMBEDDED_OBJECT<VECTOR<float,3>,2>;
template class EMBEDDED_OBJECT<VECTOR<float,3>,3>;
template class EMBEDDED_OBJECT<VECTOR<double,2>,2>;
template class EMBEDDED_OBJECT<VECTOR<double,3>,2>;
template class EMBEDDED_OBJECT<VECTOR<double,3>,3>;
}
