//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION
//#####################################################################
#ifndef __ARRAY_COLLECTION__
#define __ARRAY_COLLECTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION_ELEMENT.h>
#include <PhysBAM_Tools/Arrays/ATTRIBUTE_ID.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <memory>
namespace PhysBAM{

struct READ_WRITE_ARRAY_COLLECTION_FUNCTIONS_HELPER
{
    ARRAY<ARRAY_COLLECTION_ELEMENT_BASE*> samples;
    ~READ_WRITE_ARRAY_COLLECTION_FUNCTIONS_HELPER()
    {samples.Delete_Pointers_And_Clean_Memory();}
};

struct READ_WRITE_ARRAY_COLLECTION_FUNCTIONS
{
    std::string name;
    ARRAY_COLLECTION_ELEMENT_BASE* sample_attribute;
};

void Register_Attribute_Name(const ATTRIBUTE_ID id,const char* name);
const char* Get_Attribute_Name(const ATTRIBUTE_ID id);

class ARRAY_COLLECTION:public CLONEABLE<ARRAY_COLLECTION>,public NONCOPYABLE
{
public:
    typedef int HAS_TYPED_READ_WRITE;
    static HASHTABLE<ATTRIBUTE_ID,READ_WRITE_ARRAY_COLLECTION_FUNCTIONS>& Read_Write_Array_Collection_Registry();
    int number;
    int buffer_size;
    ARRAY<ARRAY_COLLECTION_ELEMENT_BASE*,ATTRIBUTE_INDEX> arrays;
    ARRAY<int> deletion_list;
    bool delete_data;

    ARRAY_COLLECTION();
    virtual ~ARRAY_COLLECTION();

    void Clone_Helper(const ARRAY_COLLECTION& collection)
    {Initialize(collection);}

    int Size() const
    {return number;}

    bool operator!=(const ARRAY_COLLECTION& collection) const
    {return !(*this==collection);}

    void Compact()
    {if(number!=buffer_size) Reallocate_Buffer(number);}

    void Preallocate(const int max_size)
    {if(buffer_size<max_size) Reallocate_Buffer(max_size);}

    int Add_Element()
    {Resize(number+1);return number-1;}

    ARRAY_PLUS_SCALAR<int,IDENTITY_ARRAY<> > Add_Elements(const int new_element)
    {int old_number=number;Resize(number+new_element);
    return IDENTITY_ARRAY<>(new_element)+old_number;}

    int Add_Element_From_Deletion_List()
    {return deletion_list.m?deletion_list.Pop():Add_Element();}

    void Delete_Element(const int p)
    {Copy_Element_Helper(number-1,p);Resize(number-1);}

    void Delete_All_Elements()
    {Resize(0);}

    void Add_To_Deletion_List(const int p)
    {assert((unsigned)p<(unsigned)number);deletion_list.Append(p);}

    ATTRIBUTE_INDEX Get_Attribute_Index(const ATTRIBUTE_ID attribute_id) const
    {ATTRIBUTE_INDEX index=Find_Attribute_Index(attribute_id);if(index<arrays.m && arrays(index)->id==attribute_id) return index;return ATTRIBUTE_INDEX(-1);}

    template<class T> ARRAY_VIEW<T>* Get_Array(const ATTRIBUTE_ID attribute_id)
    {ATTRIBUTE_INDEX index=Get_Attribute_Index(attribute_id); if(index!=ATTRIBUTE_INDEX(-1)) return Get_Array_From_Index<T>(index);return 0;}

    template<class T> const ARRAY_VIEW<T>* Get_Array(const ATTRIBUTE_ID attribute_id) const
    {ATTRIBUTE_INDEX index=Get_Attribute_Index(attribute_id); if(index!=ATTRIBUTE_INDEX(-1)) return Get_Array_From_Index<T>(index);return 0;}

    template<class T> ARRAY_VIEW<T>* Get_Array_From_Index(const ATTRIBUTE_INDEX attribute_index)
    {return dynamic_cast<ARRAY_COLLECTION_ELEMENT<T>*>(arrays(attribute_index))->array;}

    template<class T> const ARRAY_VIEW<T>* Get_Array_From_Index(const ATTRIBUTE_INDEX attribute_index) const
    {return dynamic_cast<ARRAY_COLLECTION_ELEMENT<T>*>(arrays(attribute_index))->array;}

    template<class T> ATTRIBUTE_INDEX Add_Array(const ATTRIBUTE_ID attribute_id,ARRAY_VIEW<T>* array)
    {if(ARRAY_VIEW<T>* existing=Get_Array<T>(attribute_id)){PHYSBAM_ASSERT(array==existing);return Get_Attribute_Index(attribute_id);}
    return Add_Array(attribute_id,new ARRAY_COLLECTION_ELEMENT<T>(array));}

    template<class T> ATTRIBUTE_INDEX Add_Array(const ATTRIBUTE_ID attribute_id)
    {ATTRIBUTE_INDEX index=Get_Attribute_Index(attribute_id);
    if(index!=ATTRIBUTE_INDEX(-1)) return index;
    return Add_Array(attribute_id,new ARRAY_COLLECTION_ELEMENT<T>);}

    virtual void Remove_Array_Using_Index(const ATTRIBUTE_INDEX attribute_index)
    {delete arrays(attribute_index);arrays.Remove_Index(attribute_index);}

    void Remove_Array(const ATTRIBUTE_ID attribute_id)
    {ATTRIBUTE_INDEX index=Get_Attribute_Index(attribute_id);if(index!=ATTRIBUTE_INDEX(-1)) Remove_Array_Using_Index(index);}

    ATTRIBUTE_INDEX Number_Of_Arrays() const
    {return arrays.m;}

    template<class T_ARRAY_COLLECTION> void
    Initialize(const ARRAY_VIEW<T_ARRAY_COLLECTION*>& elements_per_cell)
    {PHYSBAM_ASSERT(static_cast<void*>(static_cast<T_ARRAY_COLLECTION*>(0))==static_cast<ARRAY_COLLECTION*>(0)); // make sure the following cast is valid
    Initialize(ARRAY_VIEW<const ARRAY_COLLECTION* const>(elements_per_cell.Size(),reinterpret_cast<const ARRAY_COLLECTION* const*>(elements_per_cell.Get_Array_Pointer())));}

    template<class T_ARRAY_COLLECTION> void
    Initialize(const ARRAY<T_ARRAY_COLLECTION*>& elements_per_cell)
    {Initialize(static_cast<const ARRAY_VIEW<T_ARRAY_COLLECTION*>&>(elements_per_cell));}

    template<class T_ARRAY_COLLECTION,class T_ARRAYS> void
    Initialize(const ARRAY_BASE<T_ARRAY_COLLECTION*,T_ARRAYS,typename T_ARRAYS::INDEX>& elements_per_cell)
    {Initialize(elements_per_cell.array);}

    int Append(const ARRAY_COLLECTION& source,int from)
    {Copy_Element(source,from,Add_Element());return number-1;}

    void Append(const ARRAY_COLLECTION& from_elements)
    {int offset=number;Add_Elements(from_elements.number);Copy_All_Elements_Helper(from_elements,offset);}

    int Take(ARRAY_COLLECTION& source,int from)
    {Append(source,from);source.Delete_Element(from);return number;}

    void Take(ARRAY_COLLECTION& source)
    {Append(source);source.Delete_All_Elements();}

    static ATTRIBUTE_ID Type_Only(ATTRIBUTE_ID id)
    {return ATTRIBUTE_ID(Value(id)&0xFFFF0000);}

    static ATTRIBUTE_ID Id_Only(ATTRIBUTE_ID id)
    {return ATTRIBUTE_ID(Value(id)&0x0000FFFF);}

    template<class E> static void Register_Read_Write()
    {static READ_WRITE_ARRAY_COLLECTION_FUNCTIONS_HELPER sample_helper;
    READ_WRITE_ARRAY_COLLECTION_FUNCTIONS functions;
    functions.sample_attribute=new ARRAY_COLLECTION_ELEMENT<E>();
    printf("Register: %s\n", typeid(E).name());
    functions.sample_attribute->id=ATTRIBUTE_ID();
    PHYSBAM_ASSERT(Read_Write_Array_Collection_Registry().Set(Type_Only(functions.sample_attribute->Hashed_Id()),functions));
    sample_helper.samples.Append(functions.sample_attribute);}

//#####################################################################
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
    void Print(std::ostream& output,const int p) const;
    virtual void Clean_Memory();
    bool operator==(const ARRAY_COLLECTION& collection) const;
    virtual void Resize(const int new_size);
    ATTRIBUTE_INDEX Add_Array(const ATTRIBUTE_ID attribute_id,ARRAY_COLLECTION_ELEMENT_BASE* array);
    void Initialize(const ARRAY_COLLECTION& elements);
    void Initialize(ARRAY_VIEW<const ARRAY_COLLECTION* const> elements_per_cell);
    void Add_Arrays(const ARRAY_COLLECTION& elements);
    void Add_Elements_From_Deletion_List(const int count,ARRAY<int>& added_indices);
    void Delete_Elements_On_Deletion_List(const bool preserve_order=false);
    void Copy_Element(const ARRAY_COLLECTION& from_elements,const int from,const int to);
    void Copy_All_Elements_Helper(const ARRAY_COLLECTION& from_elements,const int offset);
    ATTRIBUTE_INDEX Find_Attribute_Index(const ATTRIBUTE_ID attribute_id) const;
    int Pack_Size() const;
    void Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const;
    void Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p);
protected:
    void Copy_Element_Helper(const int from,const int to);
    virtual void Reallocate_Buffer(int new_size);
//#####################################################################
};
}
#endif
