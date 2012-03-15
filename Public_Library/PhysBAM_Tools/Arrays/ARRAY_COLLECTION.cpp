//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <sstream>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
ARRAY_COLLECTION::
ARRAY_COLLECTION()
    :number(0),buffer_size(0),delete_data(true)
{}
//#####################################################################
// Constructor
//#####################################################################
ARRAY_COLLECTION::
~ARRAY_COLLECTION()
{
    if(delete_data) for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++) delete arrays(i);
}
//#####################################################################
// Function Initialize
//#####################################################################
void ARRAY_COLLECTION::
Initialize(const ARRAY_COLLECTION& elements)
{
    Clean_Memory();
    Add_Arrays(elements);
    Append(elements);
}
//#####################################################################
// Function Initialize
//#####################################################################
void ARRAY_COLLECTION::
Initialize(ARRAY_VIEW<const ARRAY_COLLECTION* const> elements_per_cell)
{
    Clean_Memory();
    int total_number=0;for(int c=0;c<elements_per_cell.Size();c++) if(elements_per_cell(c)){
        total_number+=elements_per_cell(c)->number;
        Add_Arrays(*elements_per_cell(c));} // include arrays that occur on any of the cell elements
    Preallocate(total_number);
    for(int c=0;c<elements_per_cell.Size();c++) if(elements_per_cell(c)) Append(*elements_per_cell(c));
}
//#####################################################################
// Function Add_Arrays
//#####################################################################
void ARRAY_COLLECTION::
Add_Arrays(const ARRAY_COLLECTION& collection)
{
    ATTRIBUTE_INDEX i(0),j(0);
    for(;i<arrays.m && j<collection.arrays.m;i++){
        if(arrays(i)->id<collection.arrays(j)->id) continue;
        if(arrays(i)->id>collection.arrays(j)->id) Add_Array(collection.arrays(j)->id,collection.arrays(j)->Clone_Default());
        j++;}
    for(;j<collection.arrays.m;j++)
        Add_Array(collection.arrays(j)->id,collection.arrays(j)->Clone_Default());
}
//#####################################################################
// Function Add_Elements_From_Deletion_List
//#####################################################################
void ARRAY_COLLECTION::
Add_Elements_From_Deletion_List(const int count,ARRAY<int>& added_indices)
{
    added_indices.Preallocate(added_indices.Size()+count);
    int added=min(deletion_list.m,count);
    added_indices.Append_Elements(deletion_list.Pop_Elements(added));
    added_indices.Append_Elements(Add_Elements(count-added));
}
//#####################################################################
// Function Delete_Elements_On_Deletion_List
//#####################################################################
void ARRAY_COLLECTION::
Delete_Elements_On_Deletion_List(const bool preserve_order)
{
    Sort(deletion_list);
    if(preserve_order){
        for(int k=0;k<deletion_list.m;k++){
            int next=k<deletion_list.m-1?deletion_list(k+1):number;
            for(int i=deletion_list(k)+1;i<next;i++) Copy_Element_Helper(i,i-k-1);}}
    else{
        int last=number;
        for(int k=deletion_list.m-1;k>=0;k--)
            Copy_Element_Helper(--last,deletion_list(k));}
    Resize(number-deletion_list.m);
    deletion_list.Remove_All();
}
//#####################################################################
// Function Copy_Element_Helper
//#####################################################################
void ARRAY_COLLECTION::
Copy_Element(const ARRAY_COLLECTION& from_collection,const int from,const int to)
{
    ATTRIBUTE_INDEX i(0),j(0);
    while(i<arrays.m && j<from_collection.arrays.m){
        if(arrays(i)->id<from_collection.arrays(j)->id) arrays(i++)->Clear(to);
        else if(arrays(i)->id>from_collection.arrays(j)->id) j++;
        else arrays(i++)->Copy_Element(*from_collection.arrays(j++),from,to);}
    for(;i<arrays.m;i++) arrays(i)->Clear(to);
}
//#####################################################################
// Function Copy_All_Elements_Helper
//#####################################################################
void ARRAY_COLLECTION::
Copy_All_Elements_Helper(const ARRAY_COLLECTION& from_collection,const int offset)
{
    PHYSBAM_ASSERT(this!=&from_collection);
    ATTRIBUTE_INDEX i(0),j(0);
    while(i<arrays.m && j<from_collection.arrays.m){
        if(arrays(i)->id<from_collection.arrays(j)->id) arrays(i++)->Clear_Range(offset+1,offset+from_collection.number);
        else if(arrays(i)->id>from_collection.arrays(j)->id) j++;
        else arrays(i++)->Copy_With_Offset(*from_collection.arrays(j++),offset);}
    for(;i<arrays.m;i++) arrays(i)->Clear_Range(offset+1,offset+from_collection.number);
}
//#####################################################################
// Function Get_Attribute_Index
//#####################################################################
ATTRIBUTE_INDEX ARRAY_COLLECTION::
Find_Attribute_Index(const ATTRIBUTE_ID attribute_id) const
{
    ATTRIBUTE_INDEX first(0),last(arrays.m);
    while(first<last){
        ATTRIBUTE_INDEX middle((Value(first)+Value(last))/2);
        if(arrays(middle)->id<attribute_id) first=middle+1;
        else last=middle;}
    return first;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
void ARRAY_COLLECTION::
Clean_Memory()
{
    number=buffer_size=0;
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Clean_Memory();
}
//#####################################################################
// operator==
//#####################################################################
bool ARRAY_COLLECTION::
operator==(const ARRAY_COLLECTION& collection) const
{
    if(this==&collection) return true;
    if(arrays.m!=collection.arrays.m) return false;
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++){
        if(arrays(i)->id!=collection.arrays(i)->id) return false;
        assert(typeid(arrays(i))==typeid(collection.arrays(i)));
        if(arrays(i)!=collection.arrays(i)) return false;}
    return true;
}
//#####################################################################
// Function Resize
//#####################################################################
void ARRAY_COLLECTION::
Resize(const int new_size)
{
    if(buffer_size<new_size) Reallocate_Buffer(max(4*number/3+2,new_size));
    number=new_size;
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++) arrays(i)->Set_Size(number);
}
//#####################################################################
// Function Add_Array
//#####################################################################
ATTRIBUTE_INDEX ARRAY_COLLECTION::
Add_Array(const ATTRIBUTE_ID attribute_id,ARRAY_COLLECTION_ELEMENT_BASE* array)
{
    ATTRIBUTE_INDEX index=Find_Attribute_Index(attribute_id);
    if(index<arrays.m && arrays(index)->id==attribute_id){
        PHYSBAM_ASSERT(array==arrays(index));
        return index;}
    array->id=attribute_id;
    array->Reallocate(buffer_size);
    array->Set_Size(number);
    arrays.Insert(array,index);
    return index;
}
//#####################################################################
// Function Pack_Size
//#####################################################################
int ARRAY_COLLECTION::
Pack_Size() const
{
    int pack_size=0;
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        pack_size+=arrays(i)->Pack_Size();
    return pack_size;
}
//#####################################################################
// Function Pack
//#####################################################################
void ARRAY_COLLECTION::
Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const
{
    assert((unsigned)p<(unsigned)number);
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Pack(buffer,position,p);
}
//#####################################################################
// Function Unpack
//#####################################################################
void ARRAY_COLLECTION::
Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p)
{
    assert((unsigned)p<(unsigned)number);
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Unpack(buffer,position,p);
}
//#####################################################################
// Function Copy_Element_Helper
//#####################################################################
void ARRAY_COLLECTION::
Copy_Element_Helper(const int from,const int to)
{
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Copy_Element(from,to);
}
//#####################################################################
// Function Reallocate_Buffer
//#####################################################################
void ARRAY_COLLECTION::
Reallocate_Buffer(int new_size)
{
    buffer_size=new_size;
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Reallocate(buffer_size);
}
//#####################################################################
// Function Read_Write_Array_Collection_Registry
//#####################################################################
HASHTABLE<ATTRIBUTE_ID,READ_WRITE_ARRAY_COLLECTION_FUNCTIONS>& ARRAY_COLLECTION::
Read_Write_Array_Collection_Registry()
{
    static HASHTABLE<ATTRIBUTE_ID,READ_WRITE_ARRAY_COLLECTION_FUNCTIONS> read_write_array_collection_registry;
    return read_write_array_collection_registry;
}
//#####################################################################
// Function Attribute_Names_Registry
//#####################################################################
static HASHTABLE<ATTRIBUTE_ID,const char*>& Attribute_Names_Registry()
{
    static HASHTABLE<ATTRIBUTE_ID,const char*> names_registry;
    return names_registry;
}
//#####################################################################
// Function Register_Attribute_Name
//#####################################################################
void PhysBAM::Register_Attribute_Name(const ATTRIBUTE_ID id,const char* name)
{
    PHYSBAM_ASSERT(Attribute_Names_Registry().Set(id,name));
}
//#####################################################################
// Function Register_Attribute_Name
//#####################################################################
const char* PhysBAM::Get_Attribute_Name(const ATTRIBUTE_ID id)
{
    if(const char** name=Attribute_Names_Registry().Get_Pointer(id)) return *name;
    return 0;
}
//#####################################################################
// Function Read_Arrays
//#####################################################################
void ARRAY_COLLECTION::
Read(TYPED_ISTREAM& input)
{
    int size;
    Read_Binary(input,size);
    if(size<0) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Invalid negative size %d",size));
    Clean_Memory();
    Resize(size);
    ATTRIBUTE_INDEX num_attributes;
    Read_Binary(input,size,num_attributes);
    Resize(size);

    for(ATTRIBUTE_INDEX i(0);i<num_attributes;i++){
        ATTRIBUTE_ID hashed_id;int read_size;
        Read_Binary(input,hashed_id,read_size);

        READ_WRITE_ARRAY_COLLECTION_FUNCTIONS* read_write_functions=Read_Write_Array_Collection_Registry().Get_Pointer(Type_Only(hashed_id));
        if(!read_write_functions){input.stream.ignore(read_size);continue;}

        ATTRIBUTE_INDEX index=Get_Attribute_Index(Id_Only(hashed_id));
        if(index<ATTRIBUTE_INDEX()) index=Add_Array(Id_Only(hashed_id),read_write_functions->sample_attribute->Clone_Default());
        // TODO: this really ought to know whether we're running in float or double
        arrays(index)->Read(input);
        arrays(index)->id=Id_Only(hashed_id);}
}
//#####################################################################
// Function Write_Arrays
//#####################################################################
void ARRAY_COLLECTION::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,Size());
    Write_Binary(output,number,arrays.m);
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++){
        const ARRAY_COLLECTION_ELEMENT_BASE* entry=arrays(i);
        int calculated_write_size=entry->Write_Size(output.type.use_doubles);
        Write_Binary(output,output.type.use_doubles?entry->Typed_Hashed_Id(0.):entry->Typed_Hashed_Id(0.f),calculated_write_size);
        if(calculated_write_size) entry->Write(output);}
}
//#####################################################################
// Function Print
//#####################################################################
void ARRAY_COLLECTION::
Print(std::ostream& output,const int p) const
{
    if(p<0 || p>=number) throw INDEX_ERROR("Index out of range");
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++) arrays(i)->Print(output,p);
}
//#####################################################################

