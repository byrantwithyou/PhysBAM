//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/TWIST.h>
#include <Tools/Particles/PARTICLES.h>
#include <sstream>
using namespace PhysBAM;
namespace{
struct ATTRIBUTE_INFO
{
    int id;
    int scalar;
    int m,n;
    int size[2];
    int Encode(bool use_doubles) const
    {
        int sc=scalar<2?use_doubles:scalar;
        return (id<<20)|(sc<<16)|(m<<12)|(n<<8)|size[use_doubles];
    }
};
}
inline int Adjust_Scalar(int id,bool use_doubles)
{
    if((id&0xf)>=2) return id;
    return (id&0xfffffff0)|use_doubles;
}
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
        if(arrays(i)->name<arrays(j)->name) continue;
        if(arrays(i)->name>arrays(j)->name) Add_Array(arrays(j)->name,arrays(j)->Clone_Default());
        j++;}
    for(;j<arrays.m;j++)
        Add_Array(arrays(j)->name,arrays(j)->Clone_Default());
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
    int p=Add_Elements(count-added);
    for(int i=p;i<count-added+p;i++) added_indices.Append(i);
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
        if(arrays(i)->name<from_particles.arrays(j)->name) arrays(i++)->Clear(to);
        else if(arrays(i)->name>from_particles.arrays(j)->name) j++;
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
        if(arrays(i)->name<from_particles.arrays(j)->name) arrays(i++)->Clear_Range(offset+1,offset+from_particles.number);
        else if(arrays(i)->name>from_particles.arrays(j)->name) j++;
        else arrays(i++)->Copy_With_Offset(*from_particles.arrays(j++),offset);}
    for(;i<arrays.m;i++) arrays(i)->Clear_Range(offset+1,offset+from_particles.number);
}
//#####################################################################
// Function Get_Attribute_Index
//#####################################################################
template<class TV> ATTRIBUTE_INDEX PARTICLES<TV>::
Find_Attribute_Index(const std::string& name) const
{
    ATTRIBUTE_INDEX first(0),last(arrays.m);
    while(first<last){
        ATTRIBUTE_INDEX middle((Value(first)+Value(last))/2);
        if(arrays(middle)->name<name) first=middle+1;
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
        if(arrays(i)->name!=arrays(i)->name) return false;
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
Add_Array(const std::string& name,ARRAY_COLLECTION_ELEMENT_BASE* array)
{
    ATTRIBUTE_INDEX index=Find_Attribute_Index(name);
    if(index<arrays.m && arrays(index)->name==name){
        PHYSBAM_ASSERT(array==arrays(index));
        return index;}
    array->name=name;
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
static HASHTABLE<int,VECTOR<ARRAY_COLLECTION_ELEMENT_BASE*,2> >&
Attribute_Sample_Registry()
{
    static HASHTABLE<int,VECTOR<ARRAY_COLLECTION_ELEMENT_BASE*,2> > registry;
    return registry;
}
//#####################################################################
// Function Attribute_Sample_Registry
//#####################################################################
static HASHTABLE<std::string,ATTRIBUTE_INFO>&
Attribute_Id_Lookup() // 0 = actual type, 1 = float version, 2 = double version
{
    static HASHTABLE<std::string,ATTRIBUTE_INFO> hash;
    return hash;
}
}
//#####################################################################
// Function Legacy_Read
//#####################################################################
template<class TV> void PARTICLES<TV>::
Legacy_Read(TYPED_ISTREAM& input)
{
    static const char* attribute_table[]={
        "","X","V","frame","twist","rigid_body","structure_ids","mass",
        "rigid_mass","angular_momentum","E","rho","age","phi","grad_phi",
        "radius","vorticity","material_volume","quantized_collision_distance",
        "one_over_mass","id","effective_mass","one_over_effective_mass",
        "rigid_inertia_tensor","kinematic","","","","","","collidable","color",
        "display_size","","","","","","","","F","volume","B","valid","S","C",
        "","","","Fp","mu","lambda","mu0","lambda0","plastic_deformation",
        "dp_rho_f","dp_cohesion","viscosity","phase"
    };
    
    int size;
    ATTRIBUTE_INDEX num_attributes;
    Read_Binary(input,size,size,num_attributes);
    if(size<0) throw READ_ERROR(LOG::sprintf("Invalid negative size %d",size));
    Clean_Memory();
    Resize(size);

    bool use_doubles=sizeof(T)==sizeof(double);
    for(ATTRIBUTE_INDEX i(0);i<num_attributes;i++){
        int hashed_id,read_size;
        Read_Binary(input,hashed_id,read_size);
        int type=hashed_id&0xFFFF0000,id=hashed_id&0x0000FFFF;

        VECTOR<ARRAY_COLLECTION_ELEMENT_BASE*,2> sample_attribute;
        for(const auto& it:Attribute_Id_Lookup())
            if(Hash(it.key)*0x10000==type){
                int coded_id=it.data.Encode(input.type.use_doubles);
                Attribute_Sample_Registry().Get(coded_id,sample_attribute);
                break;}
        PHYSBAM_ASSERT(sample_attribute(use_doubles));
        std::string name=attribute_table[id];
        ATTRIBUTE_INDEX index=Get_Attribute_Index(name);
        if(index<ATTRIBUTE_INDEX())
            index=Add_Array(name,sample_attribute(use_doubles)->Clone_Default());
        arrays(index)->Read(input);
        arrays(index)->name=name;}
}
//#####################################################################
// Function Read_Arrays
//#####################################################################
template<class TV> void PARTICLES<TV>::
Read(TYPED_ISTREAM& input)
{
    // bit 0 => doubles
    int version;
    Read_Binary(input,version);
    if(version==1) return Legacy_Read(input);
    if(version!=2) throw READ_ERROR(LOG::sprintf("Unrecognized particle version %d",(int)version));

    int size;
    ATTRIBUTE_INDEX num_attributes;
    Read_Binary(input,size,num_attributes);

    if(size<0) throw READ_ERROR(LOG::sprintf("Invalid negative size %d",size));
    Clean_Memory();
    Resize(size);

    bool use_doubles=sizeof(T)==sizeof(double);
    for(ATTRIBUTE_INDEX i(0);i<num_attributes;i++){
        std::string name;
        int coded_id;
        Read_Binary(input,coded_id,name);
        VECTOR<ARRAY_COLLECTION_ELEMENT_BASE*,2> sample_attribute;
        if(!Attribute_Sample_Registry().Get(coded_id,sample_attribute)){
            int skip_size=(coded_id&0xff)*number+sizeof(int);
            input.stream.ignore(skip_size);
            continue;}

        ATTRIBUTE_INDEX index=Get_Attribute_Index(name);
        if(index<ATTRIBUTE_INDEX())
            index=Add_Array(name,sample_attribute(use_doubles)->Clone_Default());
        int scalar=(coded_id>>16)&0xf;
        bool read_doubles=scalar<2?scalar:input.type.use_doubles;
        TYPED_ISTREAM actual_input(input.stream,STREAM_TYPE(read_doubles));
        arrays(index)->Read(actual_input);
        arrays(index)->name=name;}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void PARTICLES<TV>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,2,number,arrays.m);
    for(ATTRIBUTE_INDEX i(0);i<arrays.m;i++){
        const ARRAY_COLLECTION_ELEMENT_BASE* entry=arrays(i);
        const ATTRIBUTE_INFO& ai=Attribute_Id_Lookup().Get(entry->Type_Name());
        int coded_id=ai.Encode(output.type.use_doubles);
        Write_Binary(output,coded_id,entry->name);
        arrays(i)->Write(output);}
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
void PhysBAM::Register_Attribute_Sample(int id,int scalar,int m,int n,
    int size_scalar,int size,ARRAY_COLLECTION_ELEMENT_BASE* float_element,
    ARRAY_COLLECTION_ELEMENT_BASE* double_element)
{
    static ELEMENT_SAMPLES_HELPER sample_helper;
    float_element->name="";
    double_element->name="";
    int sf=size_scalar*sizeof(float)+size;
    int sd=size_scalar*sizeof(double)+size;
    ATTRIBUTE_INFO ai={id,scalar,m,n,{sf,sd}};
    VECTOR<ARRAY_COLLECTION_ELEMENT_BASE*,2> e={float_element,double_element};
    PHYSBAM_ASSERT(Attribute_Id_Lookup().Set(float_element->Type_Name(),ai));
    PHYSBAM_ASSERT(Attribute_Sample_Registry().Set(ai.Encode(false),e));
    sample_helper.samples.Append(float_element);
    if(double_element!=float_element){
        PHYSBAM_ASSERT(Attribute_Id_Lookup().Set(double_element->Type_Name(),ai));
        PHYSBAM_ASSERT(Attribute_Sample_Registry().Set(ai.Encode(true),e));
        sample_helper.samples.Append(double_element);}
}
template<class T,int d> void Register_Attributes_ts(int s)
{
}
template<int d> void Register_Attributes_s()
{
    Register_Attribute_Sample<VECTOR<float,d>,VECTOR<double,d> >(2,0,d,d,0);
    Register_Attribute_Sample<MATRIX<float,d,d>,MATRIX<double,d,d> >(3,d,d,d*d,0);
    Register_Attribute_Sample<DIAGONAL_MATRIX<float,d>,DIAGONAL_MATRIX<double,d> >(4,d,d,d,0);
    Register_Attribute_Sample<VECTOR<int,d> >(2,2,0,d,d*sizeof(int));
}
template<class T,int d> void Register_Attributes_ts(int s,int fs,int ts)
{
}
template<int d> void Register_Attributes_s(int fs,int ts)
{
    Register_Attribute_Sample<FRAME<VECTOR<float,d> >,FRAME<VECTOR<double,d> > >(5,0,d,fs,0);
    Register_Attribute_Sample<TWIST<VECTOR<float,d> >,TWIST<VECTOR<double,d> > >(6,0,d,ts,0);
}
//#####################################################################
// Function Register_Attributes
//#####################################################################
static int Register_Attributes()
{
    Register_Attributes_s<0>();
    Register_Attributes_s<1>();
    Register_Attributes_s<2>();
    Register_Attributes_s<3>();
    Register_Attributes_s<1>(1,1);
    Register_Attributes_s<2>(4,3);
    Register_Attributes_s<3>(7,6);

    Register_Attribute_Sample<float,double>(1,0,0,1,0);
    Register_Attribute_Sample<int>(1,2,0,0,sizeof(int));
    Register_Attribute_Sample<bool>(1,3,0,0,sizeof(bool));
    Register_Attribute_Sample<unsigned short>(1,4,0,0,sizeof(unsigned short));

    return 1;
}
int register_attributes=Register_Attributes();
//#####################################################################
namespace PhysBAM{
template class PARTICLES<VECTOR<float,1> >;
template class PARTICLES<VECTOR<float,2> >;
template class PARTICLES<VECTOR<float,3> >;
template class PARTICLES<VECTOR<double,1> >;
template class PARTICLES<VECTOR<double,2> >;
template class PARTICLES<VECTOR<double,3> >;
}
