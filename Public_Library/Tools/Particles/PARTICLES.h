//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLES
//#####################################################################
#ifndef __PARTICLES__
#define __PARTICLES__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Arrays_Nd/ARRAYS_ND_BASE.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Tools/Clone/CLONEABLE.h>
#include <Tools/Particles/ARRAY_COLLECTION_ELEMENT.h>
#include <Tools/Particles/ARRAY_COLLECTION_ELEMENT_BASE.h>
#include <Tools/Particles/PARTICLES_FORWARD.h>
#include <Grid_Tools/Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <memory>
namespace PhysBAM{

void Register_Attribute_Name(const ATTRIBUTE_ID id,const char* name);
const char* Get_Attribute_Name(const ATTRIBUTE_ID id);

void Register_Attribute_Sample(ARRAY_COLLECTION_ELEMENT_BASE* element);

void Register_Attribute_Sample(ARRAY_COLLECTION_ELEMENT_BASE* element_float,ARRAY_COLLECTION_ELEMENT_BASE* element_double);

template<class E> void Register_Attribute_Sample() {Register_Attribute_Sample(new ARRAY_COLLECTION_ELEMENT<E>);}
template<class F,class D> void Register_Attribute_Sample() {Register_Attribute_Sample(new ARRAY_COLLECTION_ELEMENT<F>,new ARRAY_COLLECTION_ELEMENT<D>);}

template<class TV>
class PARTICLES:public CLONEABLE<PARTICLES<TV> >,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    typedef T SCALAR;
    typedef TV VECTOR_T;
    typedef int HAS_TYPED_READ_WRITE;
    int number;
    int buffer_size;
    ARRAY<ARRAY_COLLECTION_ELEMENT_BASE*,ATTRIBUTE_INDEX> arrays;
    ARRAY<int> deletion_list;
    bool delete_data;

    PARTICLES();
    virtual ~PARTICLES();

    void Clone_Helper(const PARTICLES& particles)
    {Initialize(particles);}

    bool operator!=(const PARTICLES<TV>& particles) const
    {return !(*this==particles);}

    int Size() const
    {return number;}

    void Compact()
    {if(number!=buffer_size) Reallocate_Buffer(number);}

    void Preallocate(const int max_size)
    {if(buffer_size<max_size) Reallocate_Buffer(max_size);}

    int Add_Element()
    {Resize(number+1);return number-1;}

    int Add_Elements(const int new_element)
    {int old_number=number;Resize(number+new_element);return old_number;}

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

    template<class E> ARRAY_VIEW<E>* Get_Array(const ATTRIBUTE_ID attribute_id)
    {ATTRIBUTE_INDEX index=Get_Attribute_Index(attribute_id); if(index!=ATTRIBUTE_INDEX(-1)) return Get_Array_From_Index<E>(index);return 0;}

    template<class E> const ARRAY_VIEW<E>* Get_Array(const ATTRIBUTE_ID attribute_id) const
    {ATTRIBUTE_INDEX index=Get_Attribute_Index(attribute_id); if(index!=ATTRIBUTE_INDEX(-1)) return Get_Array_From_Index<E>(index);return 0;}

    template<class E> ARRAY_VIEW<E>* Get_Array_From_Index(const ATTRIBUTE_INDEX attribute_index)
    {return dynamic_cast<ARRAY_COLLECTION_ELEMENT<E>*>(arrays(attribute_index))->array;}

    template<class E> const ARRAY_VIEW<E>* Get_Array_From_Index(const ATTRIBUTE_INDEX attribute_index) const
    {return dynamic_cast<ARRAY_COLLECTION_ELEMENT<E>*>(arrays(attribute_index))->array;}

    template<class E> ATTRIBUTE_INDEX Add_Array(const ATTRIBUTE_ID attribute_id,ARRAY_VIEW<E>* array)
    {if(ARRAY_VIEW<E>* existing=Get_Array<E>(attribute_id)){PHYSBAM_ASSERT(array==existing);return Get_Attribute_Index(attribute_id);}
    return Add_Array(attribute_id,new ARRAY_COLLECTION_ELEMENT<E>(array));}

    template<class E> ATTRIBUTE_INDEX Add_Array(const ATTRIBUTE_ID attribute_id)
    {ATTRIBUTE_INDEX index=Get_Attribute_Index(attribute_id);
    if(index!=ATTRIBUTE_INDEX(-1)) return index;
    return Add_Array(attribute_id,new ARRAY_COLLECTION_ELEMENT<E>);}

    void Remove_Array(const ATTRIBUTE_ID attribute_id)
    {ATTRIBUTE_INDEX index=Get_Attribute_Index(attribute_id);if(index!=ATTRIBUTE_INDEX(-1)) Remove_Array_Using_Index(index);}

    ATTRIBUTE_INDEX Number_Of_Arrays() const
    {return arrays.m;}

    template<class T_PARTICLES> void
    Initialize(const ARRAY_VIEW<T_PARTICLES*>& elements_per_cell)
    {PHYSBAM_ASSERT(static_cast<void*>(static_cast<T_PARTICLES*>(0))==static_cast<PARTICLES*>(0)); // make sure the following cast is valid
    Initialize(ARRAY_VIEW<const PARTICLES* const>(elements_per_cell.Size(),reinterpret_cast<const PARTICLES* const*>(elements_per_cell.Get_Array_Pointer())));}

    template<class T_PARTICLES> void
    Initialize(const ARRAY<T_PARTICLES*>& elements_per_cell)
    {Initialize(static_cast<const ARRAY_VIEW<T_PARTICLES*>&>(elements_per_cell));}

    template<class T_PARTICLES,class T_ARRAYS> void
    Initialize(const ARRAY_BASE<T_PARTICLES*,T_ARRAYS,typename T_ARRAYS::INDEX>& elements_per_cell)
    {Initialize(elements_per_cell.array);}

    int Append(const PARTICLES& source,int from)
    {Copy_Element(source,from,Add_Element());return number-1;}

    void Append(const PARTICLES& from_elements)
    {int offset=number;Add_Elements(from_elements.number);Copy_All_Elements_Helper(from_elements,offset);}

    int Take(PARTICLES& source,int from)
    {Append(source,from);source.Delete_Element(from);return number;}

    void Take(PARTICLES& source)
    {Append(source);source.Delete_All_Elements();}

//#####################################################################
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
    void Print(std::ostream& output,const int p) const;
    virtual void Clean_Memory();
    bool operator==(const PARTICLES& particles) const;
    virtual void Resize(const int new_size);
    ATTRIBUTE_INDEX Add_Array(const ATTRIBUTE_ID attribute_id,ARRAY_COLLECTION_ELEMENT_BASE* array);
    void Initialize(const PARTICLES& elements);
    void Initialize(ARRAY_VIEW<const PARTICLES* const> elements_per_cell);
    void Add_Arrays(const PARTICLES& elements);
    void Add_Elements_From_Deletion_List(const int count,ARRAY<int>& added_indices);
    void Delete_Elements_On_Deletion_List(const bool preserve_order=false);
    void Copy_Element(const PARTICLES& from_elements,const int from,const int to);
    void Copy_All_Elements_Helper(const PARTICLES& from_elements,const int offset);
    ATTRIBUTE_INDEX Find_Attribute_Index(const ATTRIBUTE_ID attribute_id) const;
    int Pack_Size() const;
    void Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const;
    void Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p);
    virtual void Remove_Array_Using_Index(const ATTRIBUTE_INDEX attribute_index);
protected:
    void Copy_Element_Helper(const int from,const int to);
    virtual void Reallocate_Buffer(int new_size);
};
}
#endif
