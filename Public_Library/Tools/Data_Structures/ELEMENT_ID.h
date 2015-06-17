//#####################################################################
// Copyright 2008, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELEMENT_ID
//#####################################################################
#ifndef __ELEMENT_ID__
#define __ELEMENT_ID__

#include <Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Tools/Data_Structures/HASH_REDUCE.h>
#include <Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

namespace ELEMENT_ID_HELPER{
enum {none=0,equality=1,compare=2,increment=4,add_T=8,to_bool=16,negate=32,for_loop=compare|increment,logical=equality|to_bool};
}

template<class ID,class T,int flags>
class ELEMENT_ID
{
    T id_value;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T VALUE;
    enum WORKAROUND{capability_flags=flags};

    ELEMENT_ID()
        :id_value(0)
    {}

    explicit ELEMENT_ID(T n)
        :id_value(n)
    {}

    T Value() const
    {return id_value;}

    bool operator==(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::equality);return id_value==id.id_value;}

    bool operator!=(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::equality);return id_value!=id.id_value;}

    bool operator<(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::compare);return id_value<id.id_value;}

    bool operator>(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::compare);return id_value>id.id_value;}

    bool operator<=(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::compare);return id_value<=id.id_value;}

    bool operator>=(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::compare);return id_value>=id.id_value;}

    ID operator++()
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::increment);id_value++;return static_cast<ID&>(*this);}

    ID operator++(int)
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::increment);return ID(id_value++);}

    ID operator--()
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::increment);id_value--;return static_cast<ID&>(*this);}

    ID operator--(int)
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::increment);return ID(id_value--);}

    ID operator+(T i) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::add_T);return ID(id_value+i);}

    ID& operator+=(T i)
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::add_T);id_value+=i;return (ID&)*this;}

    ID operator-(T i) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::add_T);return ID(id_value-i);}

    ID& operator-=(T i)
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::add_T);id_value-=i;return (ID&)*this;}

private:
    struct UNUSABLE{void F(){}};
    typedef void (UNUSABLE::*SAFE_BOOL)();
public:

    operator SAFE_BOOL() const // allow conversion to bool without allowing conversion to T
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::to_bool);return id_value?&UNUSABLE::F:0;}

    ID operator-() const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::negate);return ID(-id_value);}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,id_value);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,id_value);}
//#####################################################################
};

template<class T> inline typename enable_if<is_fundamental<T>::value,T>::type
Value(T i)
{return i;}

template<class ID> inline typename ID::VALUE
Value(ID i)
{return i.Value();}

//#####################################################################

#define PHYSBAM_DECLARE_ELEMENT_ID(ID,T,flags)       \
    struct ID:public ELEMENT_ID<ID,T,flags>          \
    {                                                \
        ID(){}                                       \
        explicit ID(T n):ELEMENT_ID<ID,T,flags>(n){} \
        typedef int ELEMENT_ID_TAG;                  \
    };

PHYSBAM_DECLARE_ELEMENT_ID(INITIAL_SIZE,int,ELEMENT_ID_HELPER::equality);

//#####################################################################
template<class ID,class T,int flags> inline std::ostream& operator<<(std::ostream& output,const ELEMENT_ID<ID,T,flags> id)
{return output<<id.Value();}

template<class ID,class T,int flags> inline std::istream& operator>>(std::istream& input,ELEMENT_ID<ID,T,flags>& id)
{T i;input>>i;id=ID(i);return input;}
template<class ID_TYPE> struct HASH_REDUCE<ID_TYPE,typename FIRST<void,typename ID_TYPE::ELEMENT_ID_TAG>::TYPE>
{static int H(ID_TYPE id){return id.Value();}};
}
#endif
