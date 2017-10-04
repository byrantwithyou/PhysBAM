//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION_ELEMENT_BASE
//#####################################################################
#ifndef __ARRAY_COLLECTION_ELEMENT_BASE__
#define __ARRAY_COLLECTION_ELEMENT_BASE__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Arrays/ATTRIBUTE_ID.h>
#include <Tools/Clone/CLONEABLE.h>
#include <string>
namespace PhysBAM{

class ARRAY_COLLECTION_ELEMENT_BASE:public CLONEABLE_ABSTRACT<ARRAY_COLLECTION_ELEMENT_BASE>
{
public:
    typedef int HAS_TYPED_READ_WRITE;
    std::string name;
    bool owns_data;

    ARRAY_COLLECTION_ELEMENT_BASE();
    virtual ~ARRAY_COLLECTION_ELEMENT_BASE();

    virtual void Clean_Memory()=0;
    virtual void Set_Size(const int new_size)=0;
    virtual void Reallocate(const int new_size)=0;
    virtual void Clear(const int p)=0;
    virtual void Clear_Range(const int start,const int end)=0;
    virtual void Copy_Element(const int from,const int to)=0;
    virtual void Copy_Element(const ARRAY_COLLECTION_ELEMENT_BASE& from_attribute,const int from,const int to)=0;
    virtual void Copy_With_Offset(const ARRAY_COLLECTION_ELEMENT_BASE& from_attribute,const int offset)=0;
    virtual int Pack_Size() const=0;
    virtual void Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const=0;
    virtual void Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p)=0;
    virtual const char* Type_Name() const=0;
    virtual void Read(TYPED_ISTREAM& input)=0;
    virtual void Write(TYPED_OSTREAM& output) const=0;
    virtual void Print(std::ostream& output,const int p) const=0;
};
}
#endif
