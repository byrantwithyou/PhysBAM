//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTED_ARRAY
//#####################################################################
#ifndef __PROJECTED_ARRAY__
#define __PROJECTED_ARRAY__

#include <Core/Arrays/ARRAY_BASE.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T_ARRAY,class T_PROJECTOR> struct IS_ARRAY<PROJECTED_ARRAY<T_ARRAY,T_PROJECTOR> > {static const bool value=true;};
template<class T_ARRAY,class T_PROJECTOR> struct IS_ARRAY_VIEW<PROJECTED_ARRAY<T_ARRAY,T_PROJECTOR> > {static const bool value=true;};

template<class T_ARRAY> struct PROJECTED_ARRAY_BASE{};
template<int d> struct PROJECTED_ARRAY_BASE<VECTOR<int,d>&>{enum {m=VECTOR<int,d>::m};};
template<class T_ARRAY,class T_PROJECTOR> struct PROJECTED_ARRAY_BASE<PROJECTED_ARRAY<T_ARRAY,T_PROJECTOR>&>:public PROJECTED_ARRAY_BASE<T_ARRAY>{};

template<class T_ARRAY,class T_PROJECTOR> struct PROJECTED_ARRAY_ELEMENT
{
    typedef typename remove_const<typename remove_reference<typename T_PROJECTOR::template RESULT<typename T_ARRAY::ELEMENT>::TYPE>::type>::type TYPE;
};

//#####################################################################
// Class PROJECTED_ARRAY
//#####################################################################
template<class T_ARRAY,class T_PROJECTOR>
class PROJECTED_ARRAY:public PROJECTED_ARRAY_BASE<T_ARRAY>,public ARRAY_BASE<typename PROJECTED_ARRAY_ELEMENT<T_ARRAY,T_PROJECTOR>::TYPE,PROJECTED_ARRAY<T_ARRAY,T_PROJECTOR> >,
    private T_PROJECTOR
{
    typedef typename PROJECTED_ARRAY_ELEMENT<T_ARRAY,T_PROJECTOR>::TYPE T;
    typedef typename T_PROJECTOR::template RESULT<typename T_ARRAY::ELEMENT>::TYPE RESULT_NONCONST;
    typedef typename T_PROJECTOR::template RESULT<const typename T_ARRAY::ELEMENT>::TYPE RESULT_CONST;
public:
    typedef SIMPLE_ITERATOR<PROJECTED_ARRAY> iterator;
    typedef SIMPLE_ITERATOR<const PROJECTED_ARRAY> const_iterator;
    typedef T ELEMENT;typedef typename T_ARRAY::INDEX INDEX;

    T_ARRAY& array;

    PROJECTED_ARRAY(T_ARRAY& array)
        :array(array)
    {
        STATIC_ASSERT((is_empty<T_PROJECTOR>::value));
    }

    PROJECTED_ARRAY(T_ARRAY& array,const T_PROJECTOR& projector)
        :T_PROJECTOR(projector),array(array)
    {}

    PROJECTED_ARRAY(const PROJECTED_ARRAY<typename remove_const<T_ARRAY>::type,T_PROJECTOR>& projected_array)
        :T_PROJECTOR(projected_array.Projector()),array(projected_array.array)
    {}

    const T_PROJECTOR& Projector() const
    {return *this;}

    INDEX Size() const
    {return array.Size();}

    typename conditional<is_const<T_ARRAY>::value,RESULT_CONST,RESULT_NONCONST>::type operator()(const INDEX i)
    {return T_PROJECTOR::operator()(array(i));}

    RESULT_CONST operator()(const INDEX i) const
    {return T_PROJECTOR::operator()(array(i));}

    PROJECTED_ARRAY operator=(const PROJECTED_ARRAY& source)
    {return ARRAY_BASE<T,PROJECTED_ARRAY>::operator=(source);}

    template<class T_ARRAY1>
    PROJECTED_ARRAY operator=(const T_ARRAY1& source)
    {return ARRAY_BASE<T,PROJECTED_ARRAY>::operator=(source);}

    SIMPLE_ITERATOR<PROJECTED_ARRAY> begin()
    {return SIMPLE_ITERATOR<PROJECTED_ARRAY>(*this,0);}

    SIMPLE_ITERATOR<const PROJECTED_ARRAY> begin() const
    {return SIMPLE_ITERATOR<const PROJECTED_ARRAY>(*this,0);}

    SIMPLE_ITERATOR<PROJECTED_ARRAY> end()
    {return SIMPLE_ITERATOR<PROJECTED_ARRAY>(*this,Size());}

    SIMPLE_ITERATOR<const PROJECTED_ARRAY> end() const
    {return SIMPLE_ITERATOR<const PROJECTED_ARRAY>(*this,Size());}
};
//#####################################################################
// Class FIELD_PROJECTOR
//#####################################################################
template<class T_STRUCT,class T_FIELD,T_FIELD T_STRUCT::* field>
struct FIELD_PROJECTOR
{
    template<class U> struct RESULT{typedef typename conditional<is_const<U>::value,const T_FIELD&,T_FIELD&>::type TYPE;};

    const T_FIELD& operator()(const T_STRUCT& element) const
    {return element.*field;}

    T_FIELD& operator()(T_STRUCT& element) const
    {return element.*field;}
};
//#####################################################################
// Class INDEX_PROJECTOR
//#####################################################################
struct INDEX_PROJECTOR
{
    template<class T_ARRAY> struct RESULT{typedef typename conditional<is_const<T_ARRAY>::value,const typename T_ARRAY::ELEMENT&,typename T_ARRAY::ELEMENT&>::type TYPE;};
    int index;

    INDEX_PROJECTOR(const int index)
        :index(index)
    {}

    template<class T_ARRAY> const typename T_ARRAY::ELEMENT& operator()(const T_ARRAY& array) const
    {return array(index);}

    template<class T_ARRAY> typename T_ARRAY::ELEMENT& operator()(T_ARRAY& array) const
    {return array(index);}
};
//#####################################################################
}
#endif
