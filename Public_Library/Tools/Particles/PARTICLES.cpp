//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY_VIEW.h>
#include <Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <Tools/Log/LOG.h>
#include <Tools/Particles/PARTICLES.h>
#include <sstream>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLES<TV>::
PARTICLES()
    :number(0),buffer_size(0),delete_data(true)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLES<TV>::
~PARTICLES()
{
    if(delete_data) for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++) delete arrays(i);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void PARTICLES<TV>::
Initialize(const PARTICLES& elements)
{
    Clean_Memory();
    Add_Arrays(elements);
    Append(elements);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void PARTICLES<TV>::
Initialize(ARRAY_VIEW<const PARTICLES* const> elements_per_cell)
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
template<class TV> void PARTICLES<TV>::
Add_Arrays(const PARTICLES& particles)
{
    ATTRIBUTE_INDEX i(0),j(0);
    for(;i<arrays.m && j<arrays.m;i++){
        if(arrays(i)->id<arrays(j)->id) continue;
        if(arrays(i)->id>arrays(j)->id) Add_Array(arrays(j)->id,arrays(j)->Clone_Default());
        j++;}
    for(;j<arrays.m;j++)
        Add_Array(arrays(j)->id,arrays(j)->Clone_Default());
}
//#####################################################################
// Function Add_Elements_From_Deletion_List
//#####################################################################
template<class TV> void PARTICLES<TV>::
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
template<class TV> void PARTICLES<TV>::
Delete_Elements_On_Deletion_List(const bool preserve_order)
{
    deletion_list.Sort();
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
template<class TV> void PARTICLES<TV>::
Copy_Element(const PARTICLES& from_particles,const int from,const int to)
{
    ATTRIBUTE_INDEX i(0),j(0);
    while(i<arrays.m && j<from_particles.arrays.m){
        if(arrays(i)->id<from_particles.arrays(j)->id) arrays(i++)->Clear(to);
        else if(arrays(i)->id>from_particles.arrays(j)->id) j++;
        else arrays(i++)->Copy_Element(*from_particles.arrays(j++),from,to);}
    for(;i<arrays.m;i++) arrays(i)->Clear(to);
}
//#####################################################################
// Function Copy_All_Elements_Helper
//#####################################################################
template<class TV> void PARTICLES<TV>::
Copy_All_Elements_Helper(const PARTICLES& from_particles,const int offset)
{
    PHYSBAM_ASSERT(this!=&from_particles);
    ATTRIBUTE_INDEX i(0),j(0);
    while(i<arrays.m && j<from_particles.arrays.m){
        if(arrays(i)->id<from_particles.arrays(j)->id) arrays(i++)->Clear_Range(offset+1,offset+from_particles.number);
        else if(arrays(i)->id>from_particles.arrays(j)->id) j++;
        else arrays(i++)->Copy_With_Offset(*from_particles.arrays(j++),offset);}
    for(;i<arrays.m;i++) arrays(i)->Clear_Range(offset+1,offset+from_particles.number);
}
//#####################################################################
// Function Get_Attribute_Index
//#####################################################################
template<class TV> ATTRIBUTE_INDEX PARTICLES<TV>::
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
template<class TV> void PARTICLES<TV>::
Clean_Memory()
{
    number=buffer_size=0;
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Clean_Memory();
}
//#####################################################################
// operator==
//#####################################################################
template<class TV> bool PARTICLES<TV>::
operator==(const PARTICLES& particles) const
{
    if(this==&particles) return true;
    if(arrays.m!=arrays.m) return false;
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++){
        if(arrays(i)->id!=arrays(i)->id) return false;
        assert(typeid(arrays(i))==typeid(arrays(i)));
        if(arrays(i)!=arrays(i)) return false;}
    return true;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void PARTICLES<TV>::
Resize(const int new_size)
{
    if(buffer_size<new_size) Reallocate_Buffer(max(4*number/3+2,new_size));
    number=new_size;
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++) arrays(i)->Set_Size(number);
}
//#####################################################################
// Function Add_Array
//#####################################################################
template<class TV> ATTRIBUTE_INDEX PARTICLES<TV>::
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
template<class TV> int PARTICLES<TV>::
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
template<class TV> void PARTICLES<TV>::
Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const
{
    assert((unsigned)p<(unsigned)number);
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Pack(buffer,position,p);
}
//#####################################################################
// Function Unpack
//#####################################################################
template<class TV> void PARTICLES<TV>::
Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p)
{
    assert((unsigned)p<(unsigned)number);
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Unpack(buffer,position,p);
}
//#####################################################################
// Function Copy_Element_Helper
//#####################################################################
template<class TV> void PARTICLES<TV>::
Copy_Element_Helper(const int from,const int to)
{
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Copy_Element(from,to);
}
//#####################################################################
// Function Reallocate_Buffer
//#####################################################################
template<class TV> void PARTICLES<TV>::
Reallocate_Buffer(int new_size)
{
    buffer_size=new_size;
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++)
        arrays(i)->Reallocate(buffer_size);
}
namespace PhysBAM{
//#####################################################################
// Function Attribute_Sample_Registry
//#####################################################################
static HASHTABLE<ATTRIBUTE_ID,ARRAY_COLLECTION_ELEMENT_BASE*>&
Attribute_Sample_Registry(int type=0) // 0 = actual type, 1 = float version, 2 = double version
{
    static HASHTABLE<ATTRIBUTE_ID,ARRAY_COLLECTION_ELEMENT_BASE*> registry[3];
    return registry[type];
}
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
inline ATTRIBUTE_ID Type_Only(ATTRIBUTE_ID id)
{
    return ATTRIBUTE_ID(Value(id)&0xFFFF0000);
}
inline ATTRIBUTE_ID Id_Only(ATTRIBUTE_ID id)
{
    return ATTRIBUTE_ID(Value(id)&0x0000FFFF);
}
//#####################################################################
// Function Read_Arrays
//#####################################################################
template<class TV> void PARTICLES<TV>::
Read(TYPED_ISTREAM& input)
{
    int version;
    Read_Binary(input,version);
    if(version!=1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized particle version %d",(int)version));

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

        ARRAY_COLLECTION_ELEMENT_BASE* sample_attribute=0;
        if(!Attribute_Sample_Registry(sizeof(T)==sizeof(float)?1:2).Get(Type_Only(hashed_id),sample_attribute)){
            input.stream.ignore(read_size);continue;}

        ATTRIBUTE_INDEX index=Get_Attribute_Index(Id_Only(hashed_id));
        if(index<ATTRIBUTE_INDEX()) index=Add_Array(Id_Only(hashed_id),sample_attribute->Clone_Default());
        // TODO: this really ought to know whether we're running in float or double
        arrays(index)->Read(input);
        arrays(index)->id=Id_Only(hashed_id);}
}
//#####################################################################
// Function Write_Arrays
//#####################################################################
template<class TV> void PARTICLES<TV>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,1,Size(),number,arrays.m);
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++){
        const ARRAY_COLLECTION_ELEMENT_BASE* entry=arrays(i);
        int calculated_write_size=entry->Write_Size(output.type.use_doubles);
        Write_Binary(output,output.type.use_doubles?entry->Typed_Hashed_Id(0.):entry->Typed_Hashed_Id(0.f),calculated_write_size);
        if(calculated_write_size) entry->Write(output);}
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void PARTICLES<TV>::
Print(std::ostream& output,const int p) const
{
    if(p<0 || p>=number) throw INDEX_ERROR("Index out of range");
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++) arrays(i)->Print(output,p);
}
//#####################################################################
// Function Remove_Array_Using_Index
//#####################################################################
template<class TV> void PARTICLES<TV>::
Remove_Array_Using_Index(const ATTRIBUTE_INDEX attribute_index)
{
    delete arrays(attribute_index);
    arrays.Remove_Index(attribute_index);
}
struct ELEMENT_SAMPLES_HELPER
{
    ARRAY<ARRAY_COLLECTION_ELEMENT_BASE*> samples;
    ~ELEMENT_SAMPLES_HELPER()
    {samples.Delete_Pointers_And_Clean_Memory();}
};
void PhysBAM::Register_Attribute_Sample(ARRAY_COLLECTION_ELEMENT_BASE* element)
{
    static ELEMENT_SAMPLES_HELPER sample_helper;
    element->id=ATTRIBUTE_ID();
    PHYSBAM_ASSERT(Attribute_Sample_Registry(0).Set(Type_Only(element->Hashed_Id()),element));
    PHYSBAM_ASSERT(Attribute_Sample_Registry(1).Set(Type_Only(element->Hashed_Id()),element));
    PHYSBAM_ASSERT(Attribute_Sample_Registry(2).Set(Type_Only(element->Hashed_Id()),element));
    sample_helper.samples.Append(element);
}
void PhysBAM::Register_Attribute_Sample(ARRAY_COLLECTION_ELEMENT_BASE* element_float,ARRAY_COLLECTION_ELEMENT_BASE* element_double)
{
    static ELEMENT_SAMPLES_HELPER sample_helper;
    element_float->id=ATTRIBUTE_ID();
    element_double->id=ATTRIBUTE_ID();
    PHYSBAM_ASSERT(Attribute_Sample_Registry(0).Set(Type_Only(element_float->Hashed_Id()),element_float));
    PHYSBAM_ASSERT(Attribute_Sample_Registry(1).Set(Type_Only(element_float->Hashed_Id()),element_float));
    PHYSBAM_ASSERT(Attribute_Sample_Registry(2).Set(Type_Only(element_float->Hashed_Id()),element_double));
    PHYSBAM_ASSERT(Attribute_Sample_Registry(0).Set(Type_Only(element_double->Hashed_Id()),element_double));
    PHYSBAM_ASSERT(Attribute_Sample_Registry(1).Set(Type_Only(element_double->Hashed_Id()),element_float));
    PHYSBAM_ASSERT(Attribute_Sample_Registry(2).Set(Type_Only(element_double->Hashed_Id()),element_double));
    sample_helper.samples.Append(element_float);
    sample_helper.samples.Append(element_double);
}
//#####################################################################
namespace PhysBAM{
template class PARTICLES<VECTOR<float,1> >;
template class PARTICLES<VECTOR<float,2> >;
template class PARTICLES<VECTOR<float,3> >;
template class PARTICLES<VECTOR<double,1> >;
template class PARTICLES<VECTOR<double,2> >;
template class PARTICLES<VECTOR<double,3> >;
}
