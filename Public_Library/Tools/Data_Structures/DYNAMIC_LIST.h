//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYNAMIC_LIST
//#####################################################################
// Elements in this list are identified by persistent id's (assigned sequentially as elements are added to the list). 
// At any given time, the set of active elements is a subset of all of the elements previously allocated. 
// One typically accesses elements by their unique id using Element(id), but it is also possible to acccess them using their index in the list of currently active elements using
// Active_Element(index).  The latter is less preferrable since this index is not persistent. 
// The state of the dynamic list can be read/written at any time using Read/Write. 
//#####################################################################
#ifndef __DYNAMIC_LIST__
#define __DYNAMIC_LIST__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/ELEMENT_ID.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

class DYNAMIC_LIST_CORE:public NONCOPYABLE
{
private:

    void (*deleter)(void*);
public:
    typedef int HAS_UNTYPED_READ_WRITE;

    ARRAY<void*> array; // real type: ID2 -> T*
    HASHTABLE<void*,int> pointer_to_id_map; // T* -> ID
    ARRAY<int> index_to_id_map; // ID2 -> ID
    ARRAY<int> id_to_index_map; // ID -> ID2
    ARRAY<int> deletion_list; // ID
    int next_unique_id; // ID
    mutable ARRAY<int> needs_write; // ID

    DYNAMIC_LIST_CORE(void (*deleter)(void*));
    ~DYNAMIC_LIST_CORE();

    bool Is_Active(const int id) const
    {return id_to_index_map.Valid_Index(id) && id_to_index_map(id)>=0;}

    // Reads the set of active id's and updates the index<->id maps.
    // needs_init is filled with id's of those elements which have become newly active and need to be initialized.
    template<class RW> void Read(const std::string& prefix,ARRAY<int>& needs_init);

    // Writes the set of active id's.
    // needs_write indicates which elements have not been written since their creation, so derived classes should write those out (and reset the needs_write list).
    template<class RW> void Write(const std::string& prefix) const;

    void Clean_Memory();
    void Remove_All();
    int Add_Element(void* element);
    int Reactivate_Element(void* element,const int id_number);
    void Deactivate_Element(const int id,const bool delete_element);
    int Swap_Elements(void* element,const int id_number,const int id);
    void Remove_Element(const int id,const bool delete_element,const bool allow_id_reuse);
    void Purge_Element(const int id);
    void Fill_Needs_Write();
    void Delete_And_Clear(void* pointer);
private:
    void Delete_All();
};

template<class T,class ID,class ID2=int>
class DYNAMIC_LIST:public NONCOPYABLE
{
public:
    DYNAMIC_LIST_CORE core;

    DYNAMIC_LIST()
        :core(&Delete)
    {}

    void Set_Element(const ID id,const T* element)
    {core.array(core.id_to_index_map(Value(id)))=(void*)(element);}

    T* Element(const ID id) const
    {return Cast(core.array(core.id_to_index_map(Value(id))));}

    ID2 Element_Index(const ID id) const
    {return ID2(core.id_to_index_map(Value(id)));}

    int Number_Of_Elements() const // including inactive elements
    {return core.next_unique_id-core.deletion_list.m;}

    ID Size() const
    {return ID(core.next_unique_id);}

    bool Is_Active(const ID id) const
    {return core.Is_Active(Value(id));}

    void Set_Active_Element(const ID2 index,const T* element)
    {core.array(Value(index))=(void*)(element);}

    T* Active_Element(const ID2 index) const
    {return Cast(core.array(Value(index)));}

    ID Active_Element_Id(const ID2 index) const
    {return ID(core.index_to_id_map(Value(index)));}

    ID2 Number_Of_Active_Elements() const
    {return ID2(core.array.Size());}

    void Clean_Memory()
    {core.Clean_Memory();}

    void Remove_All()
    {core.Remove_All();}

    ID Add_Element(T* element)
    {return ID(core.Add_Element(element));}

    ID Reactivate_Element(T* element,const ID id_number)
    {return ID(core.Reactivate_Element((void*)(element),Value(id_number)));}

    void Deactivate_Element(const ID id,const bool delete_element=true)
    {core.Deactivate_Element(Value(id),delete_element);}

    ID Swap_Elements(T* element,const ID id_number,const ID id)
    {return ID(core.Swap_Elements((void*)(element),Value(id_number),Value(id)));}

    void Remove_Element(const ID id,const bool delete_element=true,const bool allow_id_reuse=true)
    {core.Remove_Element(Value(id),delete_element,allow_id_reuse);}

    void Purge_Element(const ID id)
    {core.Purge_Element(Value(id));}

    ARRAY<ID>& Needs_Write() const
    {return reinterpret_cast<ARRAY<ID>&>(core.needs_write);}

    void Fill_Needs_Write() // all elements will be written during the next Write()
    {core.Fill_Needs_Write();}

    template<class RW> void Read(const std::string& prefix,ARRAY<ID>& needs_init)
    {core.template Read<RW>(prefix,reinterpret_cast<ARRAY<int>&>(needs_init));}

    template<class RW> void Write(const std::string& prefix) const
    {core.template Write<RW>(prefix);}

private:
    static T* Cast(void* pointer)
    {return (T*)pointer;}

    static void Delete(void* pointer)
    {delete (T*)pointer;}

//#####################################################################
};
}
#endif
