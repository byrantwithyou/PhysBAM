//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/SCOPE.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/LINEAR_BINDING.h>
#include <Deformables/Collisions_And_Interactions/INTERSECTING_PAIRS_VISITOR.h>
#include <Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection)
    :mpi_solids(0),deformable_body_collection(deformable_body_collection)
{
    // set parameters
    Allow_Intersections(false);Set_Allow_Intersections_Tolerance();
    // output
    Output_Number_Checked(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
~TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY()
{
    structure_geometries.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters)
{
    LOG::Stat("self collision structures",structures.m);
    Build_Collision_Geometry();
    Allow_Intersections(triangle_collision_parameters.allow_intersections);
    Set_Allow_Intersections_Tolerance(triangle_collision_parameters.allow_intersections_tolerance);
    Set_Small_Number(triangle_collision_parameters.collisions_small_number);
    Output_Number_Checked(triangle_collision_parameters.collisions_output_number_checked);
    Set_Gauss_Jacobi(triangle_collision_parameters.use_gauss_jacobi);

    if(triangle_collision_parameters.collisions_disable_repulsions_based_on_proximity_factor)
        Save_Self_Collision_Free_State();
}
//#####################################################################
// Function Build_Collision_Geometry
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Build_Collision_Geometry()
{
    structure_geometries.Delete_Pointers_And_Clean_Memory();
    structure_geometries.Resize(structures.m);
    interacting_structure_pairs.Remove_All();
    for(int k=0;k<structures.m;k++){
        structure_geometries(k)=new STRUCTURE_INTERACTION_GEOMETRY<TV>(deformable_body_collection.particles);
        if(!structure_geometries(k)->Build_Collision_Geometry(*structures(k))){
            delete structure_geometries(k);
            structure_geometries(k)=0;}}
    for(int i=0;i<structures.m;i++)
        if(structure_geometries(i))
            for(int j=i;j<structures.m;j++)
                if(structure_geometries(j))
                    interacting_structure_pairs.Append(VECTOR<int,2>(i,j));
    intersecting_point_face_pairs.Remove_All();
    intersecting_edge_edge_pairs.Remove_All();
}
//#####################################################################
// Function Build_Topological_Structure_Of_Hierarchies
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Build_Topological_Structure_Of_Hierarchies()
{
    for(int k=0;k<structure_geometries.m;k++)
        if(structure_geometries(k)){
            structure_geometries(k)->Build_Topological_Structure_Of_Hierarchies();
            if(mpi_solids) structure_geometries(k)->Update_Processor_Masks(mpi_solids->Partition(),
                mpi_solids->partition_id_from_particle_index);}
}
//#####################################################################
// Function Allow_Intersections
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Allow_Intersections(const bool allow_intersections_input)
{
    allow_intersections=allow_intersections_input;
    if(allow_intersections)
        for(int k=0;k<structure_geometries.m;k++)
            if(structure_geometries(k))
                if(!structure_geometries(k)->triangulated_surface->mesh.element_edges) structure_geometries(k)->triangulated_surface->mesh.Initialize_Element_Edges();
}
//#####################################################################
// Function Check_For_Intersection
//#####################################################################
template<class T,class TV> bool Check_For_Intersection_Helper(const TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geo,const ARRAY<STRUCTURE_INTERACTION_GEOMETRY<TV>*>& structure_geometries,const VECTOR<int,2>& pair,
    ARRAY_VIEW<const VECTOR<T,1> > X,const bool grow_thickness_to_find_first_self_intersection,const T threshold)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class T,class TV> bool Check_For_Intersection_Helper(const TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geo,const ARRAY<STRUCTURE_INTERACTION_GEOMETRY<TV>*>& structure_geometries,const VECTOR<int,2>& pair,
    ARRAY_VIEW<const VECTOR<T,2> > X,const bool grow_thickness_to_find_first_self_intersection,const T threshold)
{
    if(!structure_geometries(pair[0]) || !structure_geometries(pair[1])) return false;
    SEGMENTED_CURVE_2D<T>* segmented_curve1=structure_geometries(pair[0])->segmented_curve;
    SEGMENTED_CURVE_2D<T>* segmented_curve2=structure_geometries(pair[1])->segmented_curve;
    if(segmented_curve1 && segmented_curve2){
        if(segmented_curve1->Segment_Segment_Intersection(segmented_curve2->mesh,X,threshold)){
            LOG::cout<<"intersections found, pair = "<<pair<<", threshold = "<<threshold<<std::endl;return true;}
        else if(grow_thickness_to_find_first_self_intersection) segmented_curve1->Find_First_Segment_Segment_Intersection(segmented_curve2->mesh,X,threshold,10);}
    return false;
}
template<class T,class TV> bool Check_For_Intersection_Helper(const TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geo,const ARRAY<STRUCTURE_INTERACTION_GEOMETRY<TV>*>& structure_geometries,const VECTOR<int,2>& pair,
    ARRAY_VIEW<const VECTOR<T,3> > X,const bool grow_thickness_to_find_first_self_intersection,const T threshold)
{
    ARRAY<VECTOR<int,2> > intersecting_segment_triangle_pairs;
    if(!structure_geometries(pair[0]) || !structure_geometries(pair[1])) return false;
    for(int i=0;i<2;i++){if(i==1 && pair[0]==pair[1]) break;
        SEGMENTED_CURVE<TV>* segmented_curve=structure_geometries(pair[i])->segmented_curve;
        TRIANGULATED_SURFACE<T>* triangulated_surface=structure_geometries(pair[1-i])->triangulated_surface;
        if(segmented_curve && triangulated_surface){
            if(triangulated_surface->Segment_Triangle_Intersection(segmented_curve->mesh,X,threshold,true,&intersecting_segment_triangle_pairs)){
                LOG::cout<<"intersections found, pair = "<<pair<<", threshold = "<<threshold<<std::endl;
                for(int k=0;k<intersecting_segment_triangle_pairs.m;k++){
                    int s,t;intersecting_segment_triangle_pairs(k).Get(s,t);
                    ARRAY<int> tmp;
                    geo.deformable_body_collection.binding_list.Flatten_Indices(tmp,segmented_curve->mesh.elements(s));
                    geo.deformable_body_collection.binding_list.Flatten_Indices(tmp,triangulated_surface->mesh.elements(t));
                    LOG::cout<<"segment "<<s<<", triangle "<<t<<", segment nodes = "<<segmented_curve->mesh.elements(s)<<", triangle nodes = "<<triangulated_surface->mesh.elements(t)<<"    real = "<<tmp<<std::endl;}
                return true;}
            else if(grow_thickness_to_find_first_self_intersection) triangulated_surface->Find_First_Segment_Triangle_Intersection(segmented_curve->mesh,X,threshold,10);}}
    return false;
}
template<class TV> bool TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Check_For_Intersection(const bool grow_thickness_to_find_first_self_intersection,const T thickness,VECTOR<int,2>* interaction_pair) const
{
    if(allow_intersections) return false;
    LOG::SCOPE scope("checking for self intersection");
    ARRAY_VIEW<const TV> X(deformable_body_collection.particles.X);
    T threshold=small_number;if(thickness) threshold=thickness;
    for(int k=0;k<interacting_structure_pairs.m;k++){const VECTOR<int,2>& pair=interacting_structure_pairs(k);
        if(Check_For_Intersection_Helper(*this,structure_geometries,pair,X,grow_thickness_to_find_first_self_intersection,threshold)){
            if(interaction_pair) *interaction_pair=pair;
            return true;}}
    return false;
}
template <class TV> bool TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Check_For_Intersection(const bool grow_thickness_to_find_first_self_intersection,const T thickness,ARRAY<VECTOR<int,2> >& interaction_pairs) const
{
    if(allow_intersections) return false;
    LOG::SCOPE scope("checking for self intersection");
    ARRAY_VIEW<const TV> X(deformable_body_collection.particles.X);
    T threshold=small_number;if(thickness) threshold=thickness;
    for(int k=0;k<interacting_structure_pairs.m;k++){const VECTOR<int,2>& pair=interacting_structure_pairs(k);
        if(Check_For_Intersection_Helper(*this,structure_geometries,pair,X,grow_thickness_to_find_first_self_intersection,threshold)){
            interaction_pairs.Append(pair);}}
    return interaction_pairs.Size()>0;
}
//#####################################################################
// Function Save_Self_Collision_Free_State
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Save_Self_Collision_Free_State() // assumes mass does not change
{
    DEFORMABLE_PARTICLES<TV>& full_particles=deformable_body_collection.particles;
    X_self_collision_free=full_particles.X;
    V_self_collision_free=full_particles.V;
}
//#####################################################################
// Function Restore_Self_Collision_Free_State
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Restore_Self_Collision_Free_State()
{
    DEFORMABLE_PARTICLES<TV>& full_particles=deformable_body_collection.particles;
    full_particles.X=X_self_collision_free;
    full_particles.V=V_self_collision_free;
}
//#####################################################################
// Function Compute_Intersecting_Segment_Face_Pairs
//#####################################################################
template<class TV> template<class S>  void TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Compute_Intersecting_Pairs_Helper(COMPUTE_INTERSECTING_PAIRS_HELPER_INPUT_WHEN_D_EQ_1 input)
{
    PHYSBAM_ASSERT(d==1);
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class TV> template<class S> void TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Compute_Intersecting_Pairs_Helper(COMPUTE_INTERSECTING_PAIRS_HELPER_INPUT_WHEN_D_NE_1 input)
{
    LOG::cout<<"allowing intersections!!!"<<std::endl;
    intersecting_point_face_pairs.Remove_All();intersecting_edge_edge_pairs.Remove_All();

    // update hierarchies
    for(int k=0;k<structure_geometries.m;k++)
        if(structure_geometries(k)){
            STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*structure_geometries(k);
            if(structure.Face_Mesh_Object()) structure.Face_Hierarchy().Update_Boxes(X_self_collision_free);
            if(d==3 && structure.segmented_curve) structure.segmented_curve->hierarchy->Update_Boxes(X_self_collision_free);}

    // find intersecting pairs
    for(int pair_index=0;pair_index<interacting_structure_pairs.m;pair_index++){const VECTOR<int,2>& pair=interacting_structure_pairs(pair_index);
        if(!structure_geometries(pair[0]) || !structure_geometries(pair[1])) continue;
        for(int i=0;i<2;i++){if(i==1 && (d==2 || pair[0]==pair[1])) break;
            STRUCTURE_INTERACTION_GEOMETRY<TV>& segment_structure=*structure_geometries(pair[i]);
            STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure=*structure_geometries(pair[1-i]);
            if(!segment_structure.segmented_curve || !face_structure.Face_Mesh_Object()) continue;

            int count=0;
            INTERSECTING_PAIRS_VISITOR<TV> visitor(intersecting_point_face_pairs,intersecting_edge_edge_pairs,segment_structure,face_structure,X_self_collision_free,
                allow_intersections_tolerance,count);
            if(mpi_solids){
                BOX_VISITOR_MPI<INTERSECTING_PAIRS_VISITOR<TV> > mpi_visitor(visitor,segment_structure.segmented_curve_processor_masks,face_structure.Face_Processor_Masks());
                segment_structure.segmented_curve->hierarchy->Intersection_List(face_structure.Face_Hierarchy(),mpi_visitor,allow_intersections_tolerance);}
            else segment_structure.segmented_curve->hierarchy->Intersection_List(face_structure.Face_Hierarchy(),visitor,allow_intersections_tolerance);
            if(count) LOG::cout<<"pair "<<pair<<" has "<<count<<" intersecting segment triangle pairs"<<std::endl;}} // TODO: correct for mpi

    // synchronize if necessary
    if(mpi_solids) mpi_solids->All_Gather_Intersecting_Pairs(intersecting_point_face_pairs,intersecting_edge_edge_pairs);
}
template<class TV> void TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>::
Compute_Intersecting_Segment_Face_Pairs()
{
    Compute_Intersecting_Pairs_Helper<void>((TV*)NULL);
}
//####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(T,d) \
    template class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<VECTOR<T,d> >;
INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
}
