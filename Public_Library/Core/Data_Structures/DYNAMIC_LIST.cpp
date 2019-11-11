//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYNAMIC_LIST
//#####################################################################
#include <Core/Data_Structures/DYNAMIC_LIST.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
DYNAMIC_LIST_CORE::
DYNAMIC_LIST_CORE(void (*deleter)(void*))
    :deleter(deleter),next_unique_id(0)
{}
//#####################################################################
// Destructor
//#####################################################################
DYNAMIC_LIST_CORE::
~DYNAMIC_LIST_CORE()
{
    Delete_All();
}
//#####################################################################
// Function Delete_And_Clear
//#####################################################################
void DYNAMIC_LIST_CORE::
Delete_And_Clear(void* pointer)
{
    deleter(pointer);pointer=0;
}
//#####################################################################
// Function Delete_All
//#####################################################################
void DYNAMIC_LIST_CORE::
Delete_All()
{
    for(int i=0;i<array.m;i++) deleter(array(i));
    array.Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
void DYNAMIC_LIST_CORE::
Clean_Memory()
{
    Delete_All();
    pointer_to_id_map.Clean_Memory();
    needs_write.Clean_Memory();index_to_id_map.Clean_Memory();id_to_index_map.Clean_Memory();next_unique_id=0;deletion_list.Clean_Memory();
}
//#####################################################################
// Function Remove_All
//#####################################################################
void DYNAMIC_LIST_CORE::
Remove_All()
{
    Delete_All();
    pointer_to_id_map.Clean_Memory();
    index_to_id_map.Remove_All();id_to_index_map.Remove_All();
    needs_write.Remove_All();deletion_list.Remove_All();next_unique_id=0;
}
//#####################################################################
// Function Add_Element
//#####################################################################
int DYNAMIC_LIST_CORE::
Add_Element(void* element)
{
    int id;
    if(pointer_to_id_map.Get(element,id)){assert(array(id_to_index_map(id))==element);return id;}
    int index=array.Append(element);
    if(deletion_list.m){id=deletion_list.Pop_Value();assert(id_to_index_map(id)<0);id_to_index_map(id)=index;}
    else{id=next_unique_id++;id_to_index_map.Append(index);assert(id_to_index_map.Size()==id+1);}
    index_to_id_map.Append(id);
    needs_write.Append(id);
    pointer_to_id_map.Set(element,id);
    return id;
}
//#####################################################################
// Function Reactivate_Element
//#####################################################################
int DYNAMIC_LIST_CORE::
Reactivate_Element(void* element,const int id_number)
{
    assert(id_to_index_map(id_number)<0);
    id_to_index_map(id_number)=array.Append(element);
    needs_write.Append(id_number);
    pointer_to_id_map.Set(element,id_number);
    index_to_id_map.Append(id_number);
    PHYSBAM_ASSERT(index_to_id_map.Size()==array.Size());
    return id_number;
}
//#####################################################################
// Function Deactivate_Element
//#####################################################################
void DYNAMIC_LIST_CORE::
Deactivate_Element(const int id,const bool delete_element)
{
    int index=id_to_index_map(id);assert(index>=0);
    pointer_to_id_map.Delete(array(index));
    if(delete_element) Delete_And_Clear(array(index));
    id_to_index_map(id)=-1;
    array.Remove_Index_Lazy(index);
    index_to_id_map.Remove_Index_Lazy(index);
    if(index<array.Size()) id_to_index_map(index_to_id_map(index))=index;
}
//#####################################################################
// Function Swap_Elements
//#####################################################################
int DYNAMIC_LIST_CORE::
Swap_Elements(void* element,const int id_number,const int id)
{
    int index=id_to_index_map(id);
    assert(index>=0);
    pointer_to_id_map.Delete(array(index));
    id_to_index_map(id)=-1;
    assert(id_to_index_map(id_number)==-1);
    array(index)=element;
    needs_write.Append(id_number);
    pointer_to_id_map.Set(element,id_number);
    index_to_id_map(index)=id_number;
    id_to_index_map(id_number)=index;
    assert(index_to_id_map.Size()==array.Size());
    return id_number;
}
//#####################################################################
// Function Remove_Elements
//#####################################################################
void DYNAMIC_LIST_CORE::
Remove_Element(const int id,const bool delete_element=true,const bool allow_id_reuse=true)
{
    Deactivate_Element(id,delete_element);if(allow_id_reuse)deletion_list.Append(id);
}
//#####################################################################
// Function Remove_Elements
//#####################################################################
void DYNAMIC_LIST_CORE::
Purge_Element(const int id)
{
    assert(id_to_index_map(id)<0);
    deletion_list.Append(id);
}
//#####################################################################
// Function Fill_Needs_Write
//#####################################################################
void DYNAMIC_LIST_CORE::
Fill_Needs_Write()
{
    needs_write=index_to_id_map;
}
//#####################################################################
// Constructor
//#####################################################################
void DYNAMIC_LIST_CORE::
Read(const VIEWER_DIR& viewer_dir,const std::string& list_name,ARRAY<int>& needs_init)
{
    pointer_to_id_map.Clean_Memory();
    needs_init.Remove_All();
    needs_write.Remove_All();
    ARRAY<int> active_ids;
    char version;
    Read_From_File(viewer_dir.current_directory+"/"+list_name+"_active_ids",version,next_unique_id,active_ids);
    PHYSBAM_ASSERT(version==1);
    ARRAY<void*> new_array;
    ARRAY<bool> element_copied(array.Size());
    id_to_index_map.Resize(next_unique_id,use_init,-1);
    for(int i=0;i<active_ids.Size();i++){
        int index=id_to_index_map(active_ids(i));
        if(index>=0){
            new_array.Append(array(index));
            element_copied(index)=true;}
        else{
            new_array.Append((void*)(0));
            needs_init.Append(active_ids(i));}}
    for(int i=0;i<array.Size();i++)
        if(!element_copied(i))
            Delete_And_Clear(array(i));
    index_to_id_map.Resize(active_ids.Size());
    id_to_index_map.Fill(-1);
    for(int i=0;i<new_array.Size();i++){
        pointer_to_id_map.Set(new_array(i),active_ids(i));
        index_to_id_map(i)=active_ids(i);
        id_to_index_map(active_ids(i))=i;}
    array.Exchange(new_array);
}
//#####################################################################
// Constructor
//#####################################################################
void DYNAMIC_LIST_CORE::
Write(STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir,const std::string& list_name) const
{
    const char version=1;
    Write_To_File(stream_type,viewer_dir.current_directory+"/"+list_name+"_active_ids",version,next_unique_id,index_to_id_map);
    for(int i=needs_write.m-1;i>=0;i--) if(id_to_index_map(needs_write(i))<0) needs_write.Remove_Index_Lazy(i);
    // handle case of new element which was removed without being written
}
