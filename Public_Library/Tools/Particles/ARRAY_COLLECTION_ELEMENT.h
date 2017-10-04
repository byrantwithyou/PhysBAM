//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION_ELEMENT
//#####################################################################
#ifndef __ARRAY_COLLECTION_ELEMENT__
#define __ARRAY_COLLECTION_ELEMENT__

#include <Core/Arrays/ARRAY_VIEW.h>
#include <Tools/Clone/CLONEABLE.h>
#include <Tools/Particles/ARRAY_COLLECTION_ELEMENT_BASE.h>
namespace PhysBAM{

template<class T>
class ARRAY_COLLECTION_ELEMENT:public CLONEABLE<ARRAY_COLLECTION_ELEMENT<T>,ARRAY_COLLECTION_ELEMENT_BASE>
{
    typedef CLONEABLE<ARRAY_COLLECTION_ELEMENT<T>,ARRAY_COLLECTION_ELEMENT_BASE> BASE;
public:
    typedef int HAS_TYPED_READ_WRITE;
    using BASE::owns_data;using BASE::name;
    ARRAY_VIEW<T>* array;
    int buffer_size;

    ARRAY_COLLECTION_ELEMENT(ARRAY_VIEW<T>* array_input)
        :array(array_input),buffer_size(0)
    {owns_data=false;}

    ARRAY_COLLECTION_ELEMENT()
        :array(new ARRAY_VIEW<T>(0,0))
    {}

    ~ARRAY_COLLECTION_ELEMENT()
    {ARRAY_VIEW<T> temp(0,0);temp.Exchange(*array);delete[] temp.Get_Array_Pointer();if(owns_data) delete array;}

    void Clone_Helper(const ARRAY_COLLECTION_ELEMENT& element)
    {array=element.array;owns_data=false;}

    void Clean_Memory() override
    {Reallocate(0);}

    void Clear(const int p) override
    {(*array)(p)=T();}

    void Clear_Range(const int start,const int end) override
    {for(int i=start;i<end;i++) (*array)(i)=T();}

    void Set_Size(const int new_size) override
    {assert(buffer_size>=new_size);int n=array->m;array->m=new_size;Clear_Range(n,array->m);}

    void Reallocate(const int new_size) override
    {if(new_size<array->m) array->m=new_size;
    buffer_size=new_size;
    assert(array->m<=buffer_size);
    ARRAY_VIEW<T> temp(array->m,new T[Value(new_size)]);
    temp=*array;
    temp.Exchange(*array);
    delete[] temp.Get_Array_Pointer();}

    void Copy_Element(const int from,const int to) override
    {(*array)(to)=(*array)(from);}

    void Copy_Element(const ARRAY_COLLECTION_ELEMENT_BASE& from_attribute,const int from,const int to) override
    {(*array)(to)=(*dynamic_cast<const ARRAY_COLLECTION_ELEMENT<T>&>(from_attribute).array)(from);}

    void Copy_With_Offset(const ARRAY_COLLECTION_ELEMENT_BASE& from_attribute,const int offset) override
    {array->Copy_With_Offset(*dynamic_cast<const ARRAY_COLLECTION_ELEMENT<T>&>(from_attribute).array,offset);}

    int Pack_Size() const override
    {return array->Pack_Size();}

    void Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const override
    {array->Pack(buffer,position,p);}

    void Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p) override
    {array->Unpack(buffer,position,p);}
 
    virtual const char* Type_Name() const override
    {return typeid(ARRAY_COLLECTION_ELEMENT<T>).name();}

    virtual void Read(TYPED_ISTREAM& input) override
    {Read_Binary(input,*array);}

    virtual void Write(TYPED_OSTREAM& output) const override
    {Write_Binary(output,*array);}

    virtual void Print(std::ostream& output,const int p) const override
    {output<<name<<" = "<<(*array)(p)<<std::endl;}
};
}
#endif
