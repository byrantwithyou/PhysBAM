//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_GEOMETRY_SPATIAL_PARTITION
//##################################################################### 
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T_ARRAY,class ID> COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
COLLISION_GEOMETRY_SPATIAL_PARTITION(T_ARRAY& collision_bodies_input,T collision_body_thickness_input)
    :collision_bodies(collision_bodies_input),voxel_range(collision_bodies.Size()),hashtable(10*Value(collision_bodies.Size())),reinitialize_counter(-1),
    collision_body_thickness(collision_body_thickness_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T_ARRAY,class ID> COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
~COLLISION_GEOMETRY_SPATIAL_PARTITION()
{
    hashtable.Delete_Pointers_Stored_In_Table();
}
//#####################################################################
// Function Scene_Bounding_Box_Size
//#####################################################################
template<class TV,class T_ARRAY,class ID> typename TV::SCALAR COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Scene_Bounding_Box_Size()
{
    TV lengths=Scene_Bounding_Box().Edge_Lengths();
    return lengths.Average()+2*collision_body_thickness;
}
//#####################################################################
// Function Scene_Bounding_Box
//#####################################################################
template<class TV,class T_ARRAY,class ID> RANGE<TV> COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Scene_Bounding_Box()
{
    RANGE<TV> scene_bounding_box;
    for(ID i(0);i<collision_bodies.Size();i++){
        if(RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(collision_bodies(i)))
            if(!rigid_collision_geometry->rigid_body.rigid_body_collection.Is_Active(rigid_collision_geometry->rigid_body.particle_index)) continue;
        if(collision_bodies(i) && collision_bodies(i)->add_to_spatial_partition)
            scene_bounding_box=RANGE<TV>::Combine(scene_bounding_box,collision_bodies(i)->Axis_Aligned_Bounding_Box());}
    return scene_bounding_box;
}
//#####################################################################
// Function Average_Bounding_Box_Size
//#####################################################################
template<class TV,class T_ARRAY,class ID> typename TV::SCALAR COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Average_Bounding_Box_Size()
{
    T average_size=0;int count=0;
    for(ID i(0);i<collision_bodies.Size();i++) if(collision_bodies(i)){
        if(RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(collision_bodies(i)))
            if(!rigid_collision_geometry->rigid_body.rigid_body_collection.Is_Active(rigid_collision_geometry->rigid_body.particle_index)) continue;
        if(collision_bodies(i)->add_to_spatial_partition){count++;TV size=collision_bodies(i)->Axis_Aligned_Bounding_Box().Edge_Lengths();average_size+=size.Sum();}}
    if(!count) return 0;
    return average_size/(TV::dimension*count)+2*collision_body_thickness;
}
//#####################################################################
// Function Maximum_Bounding_Box_Size
//#####################################################################
template<class TV,class T_ARRAY,class ID> typename TV::SCALAR COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Maximum_Bounding_Box_Size()
{
    T max_size=0;
    for(ID i(0);i<collision_bodies.Size();i++) if(collision_bodies(i)){
        if(RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(collision_bodies(i)))
            if(!rigid_collision_geometry->rigid_body.rigid_body_collection.Is_Active(rigid_collision_geometry->rigid_body.particle_index)) continue;
        if(collision_bodies(i)->add_to_spatial_partition){TV size=collision_bodies(i)->Axis_Aligned_Bounding_Box().Edge_Lengths();max_size=max(max_size,size.Max());}}
    return max_size+2*collision_body_thickness;
}
//#####################################################################
// Function Compute_Voxel_Size
//#####################################################################
template<class TV,class T_ARRAY,class ID> void COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Compute_Voxel_Size(const SPATIAL_PARTITION_VOXEL_SIZE_HEURISTIC heuristic,const int number_of_boxes,const T voxel_size_scale_factor)
{
    reinitialize_counter=max(reinitialize_counter,0);
    for(ID i(0);i<collision_bodies.Size();i++) if(collision_bodies(i)){
        if(RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(collision_bodies(i)))
            if(!rigid_collision_geometry->rigid_body.rigid_body_collection.Is_Active(rigid_collision_geometry->rigid_body.particle_index)) continue;
        collision_bodies(i)->Update_Bounding_Box();}
    if(heuristic==SPATIAL_PARTITION_SCENE_SIZE) voxel_size=Scene_Bounding_Box_Size()/number_of_boxes;
    else if(heuristic==SPATIAL_PARTITION_MAX_BOX_SIZE) voxel_size=voxel_size_scale_factor*Maximum_Bounding_Box_Size();
    else if(heuristic==SPATIAL_PARTITION_AVERAGE_BOX_SIZE) voxel_size=voxel_size_scale_factor*Average_Bounding_Box_Size();
    else PHYSBAM_FATAL_ERROR();
    if(voxel_size==0) voxel_size=(T)1; // doesn't matter what value
    one_over_voxel_size=(T)1/voxel_size;
    LOG::cout<<"Collision Body Spatial Partition voxel size="<<voxel_size<<" heuristic="<<heuristic<<" # of boxes="<<number_of_boxes<<" voxel_size_scale_factor="<<voxel_size_scale_factor
        <<std::endl;
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class TV,class T_ARRAY,class ID> void COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Reinitialize()
{
    assert(reinitialize_counter>=0);
    if(reinitialize_counter%10 == 0){hashtable.Delete_Pointers_Stored_In_Table();hashtable.Remove_All();}
    else hashtable.Reset_List_Arrays_Stored_In_Table();
    reinitialize_counter++;
    voxel_range.Resize(collision_bodies.Size());
    bodies_not_in_partition.Remove_All();
    for(ID i(0);i<collision_bodies.Size();i++) if(collision_bodies(i)){
        if(RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(collision_bodies(i)))
            if(!rigid_collision_geometry->rigid_body.rigid_body_collection.Is_Active(rigid_collision_geometry->rigid_body.particle_index)) continue;
        if(collision_bodies(i)->add_to_spatial_partition){
            collision_bodies(i)->Update_Bounding_Box();
            voxel_range(i)=Voxel_Range(i);
            Insert_Into_Hashtable(i);}
        else bodies_not_in_partition.Append(i);}
}
//#####################################################################
// Function Print_Initial_Statistics
//#####################################################################
template<class TV,class T_ARRAY,class ID> void COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Print_Initial_Statistics() const
{
    int min_size=INT_MAX,max_size=0,number=0;T average_size=0;
    for(ID i(0);i<collision_bodies.Size();i++) if(collision_bodies(i) && collision_bodies(i)->add_to_spatial_partition){number++;
        int size=Number_Of_Voxels_Occupied(i);min_size=min(min_size,size);max_size=max(max_size,size);average_size+=size;}
    average_size/=number;
    PHYSBAM_DEBUG_PRINT("Spatial partition statistics",voxel_size,min_size,average_size,max_size);
    PHYSBAM_DEBUG_PRINT("Spatial partition bodies",bodies_not_in_partition.m,hashtable.Size());
}
//#####################################################################
// Function Update_Body
//#####################################################################
template<class TV,class T_ARRAY,class ID> void COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Update_Body(const ID index)
{
    if(collision_bodies(index)->add_to_spatial_partition){
        collision_bodies(index)->Update_Bounding_Box();
        RANGE<TV_INT> new_range(Voxel_Range(index));
        if(voxel_range(index)==new_range) return;
        Remove_If_Not_Still_Present(voxel_range(index),new_range,index);
        Add_If_Newly_Present(voxel_range(index),new_range,index);
        voxel_range(index)=new_range;}
}
//#####################################################################
// Function Insert_Into_Hashtable
//#####################################################################
template<class TV,class T_ARRAY,class ID> void COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Insert_Into_Hashtable(const ID index)
{
    assert(collision_bodies(index)->add_to_spatial_partition);
    GRID<TV> unused;
    for(CELL_ITERATOR<TV> iterator(unused,voxel_range(index));iterator.Valid();iterator.Next()) Add_To_Cell(iterator.Cell_Index(),index);
}
//#####################################################################
// Function Remove_From_Cell
//#####################################################################
template<class TV,class T_ARRAY,class ID> bool COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Remove_From_Cell(const TV_INT& voxel,const ID index)
{
    ARRAY<ID>* list=0;
    if(!hashtable.Get(voxel,list)) return false;
    int found=list->Find(index);
    if(found==-1) return false;
    list->Remove_Index_Lazy(found);
    return true;
}
//#####################################################################
// Function Add_To_Cell
//#####################################################################
template<class TV,class T_ARRAY,class ID> void COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Add_To_Cell(const TV_INT& voxel,const ID index)
{
    ARRAY<ID>* occupancy_list=0;
    if(!hashtable.Get(voxel,occupancy_list)){occupancy_list=new ARRAY<ID>();occupancy_list->Preallocate(5);hashtable.Insert(voxel,occupancy_list);}
    occupancy_list->Append(index);
}
//#####################################################################
// Function Remove_If_Not_Still_Present
//#####################################################################
template<class TV,class T_ARRAY,class ID> void COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Remove_If_Not_Still_Present(const RANGE<TV_INT>& old_locations,const RANGE<TV_INT>& new_locations,const ID index)
{
    if(new_locations.Contains(old_locations)) return;GRID<TV> unused;
    for(CELL_ITERATOR<TV> iterator(unused,old_locations);iterator.Valid();iterator.Next()){TV_INT voxel=iterator.Cell_Index();
        if(new_locations.Lazy_Outside(voxel)) Remove_From_Cell(voxel,index);}
}
//#####################################################################
// Function Add_If_Newly_Present
//#####################################################################
template<class TV,class T_ARRAY,class ID> void COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Add_If_Newly_Present(const RANGE<TV_INT>& old_locations,const RANGE<TV_INT>& new_locations,const ID index)
{
    if(old_locations.Contains(new_locations)) return;GRID<TV> unused;
    for(CELL_ITERATOR<TV> iterator(unused,new_locations);iterator.Valid();iterator.Next()){TV_INT voxel=iterator.Cell_Index();
        if(old_locations.Lazy_Outside(voxel)) Add_To_Cell(voxel,index);}
}
//#####################################################################
// Function Get_Potential_Collisions
//#####################################################################
template<class TV,class T_ARRAY,class ID> void COLLISION_GEOMETRY_SPATIAL_PARTITION<TV,T_ARRAY,ID>::
Get_Potential_Collisions(const ID index,const RANGE<TV_INT>& range,ARRAY<ID>& object_indices,bool only_higher_index) const
{
    object_indices.Remove_All();
    already_added.Initialize(collision_bodies.Size());
    GRID<TV> unused;
    for(CELL_ITERATOR<TV> iterator(unused,range);iterator.Valid();iterator.Next()){TV_INT voxel=iterator.Cell_Index();
        ARRAY<ID>* occupancy_list=0;
        if(hashtable.Get(voxel,occupancy_list)) for(int t=0;t<occupancy_list->m;t++){
            ID k=(*occupancy_list)(t);
            if(k>index || (!only_higher_index && k!=index))
                if(!already_added.Is_Marked_Current(k)){already_added.Mark(k);object_indices.Append(k);}}}
    for(int i=0;i<bodies_not_in_partition.m;i++){ // add bodies that aren't in spatial partition - don't need append unique
        ID k=bodies_not_in_partition(i);if(k>index || (!only_higher_index && k!=index)) object_indices.Append(k);}
}
//##################################################################### 
namespace PhysBAM{
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<float,1>,const ARRAY<COLLISION_GEOMETRY<VECTOR<float,1> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<float,2>,const ARRAY<COLLISION_GEOMETRY<VECTOR<float,2> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<float,3>,const ARRAY<COLLISION_GEOMETRY<VECTOR<float,3> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<float,1>,ARRAY<COLLISION_GEOMETRY<VECTOR<float,1> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<float,2>,ARRAY<COLLISION_GEOMETRY<VECTOR<float,2> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<float,3>,ARRAY<COLLISION_GEOMETRY<VECTOR<float,3> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<double,1>,const ARRAY<COLLISION_GEOMETRY<VECTOR<double,1> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<double,2>,const ARRAY<COLLISION_GEOMETRY<VECTOR<double,2> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<double,3>,const ARRAY<COLLISION_GEOMETRY<VECTOR<double,3> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<double,1>,ARRAY<COLLISION_GEOMETRY<VECTOR<double,1> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<double,2>,ARRAY<COLLISION_GEOMETRY<VECTOR<double,2> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
template class COLLISION_GEOMETRY_SPATIAL_PARTITION<VECTOR<double,3>,ARRAY<COLLISION_GEOMETRY<VECTOR<double,3> >*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>;
}
