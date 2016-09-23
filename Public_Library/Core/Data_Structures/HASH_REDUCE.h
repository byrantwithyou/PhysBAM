//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Geoffrey Irving, Mike Rodgers, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function HASH_REDUCE
//#####################################################################
#ifndef __HASH_REDUCE__
#define __HASH_REDUCE__

#include <Core/Utilities/TYPE_UTILITIES.h>
#include <cstring>
#include <string>
namespace PhysBAM{

static const int missing_element_hash=32138912;

inline unsigned int int_hash(unsigned int key) 
{
    STATIC_ASSERT(sizeof(int)==4);
    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;
}

inline unsigned int int_hash(unsigned int a,unsigned int b,unsigned int c)
{
    STATIC_ASSERT(sizeof(int)==4);
    a-=b;a-=c;a^=(c>>13);
    b-=c;b-=a;b^=(a<<8);
    c-=a;c-=b;c^=(b>>13);
    a-=b;a-=c;a^=(c>>12);
    b-=c;b-=a;b^=(a<<16);
    c-=a;c-=b;c^=(b>>5);
    a-=b;a-=c;a^=(c>>3);
    b-=c;b-=a;b^=(a<<10);
    c-=a;c-=b;c^=(b>>15);
    return c;
}

inline unsigned int int_hash(unsigned int a,unsigned int b)
{return int_hash(missing_element_hash,a,b);}

template<class T,class ENABLE=void> struct HASH_REDUCE;

template<> struct HASH_REDUCE<bool>{static int H(bool key){return key;}};
template<> struct HASH_REDUCE<char>{static int H(char key){return key;}};
template<> struct HASH_REDUCE<unsigned char>{static int H(unsigned char key){return key;}};
template<> struct HASH_REDUCE<short>{static int H(short key){return key;}};
template<> struct HASH_REDUCE<unsigned short>{static int H(unsigned short key){return key;}};
template<> struct HASH_REDUCE<int>{static int H(int key){return key;}};
template<> struct HASH_REDUCE<unsigned int>{static int H(unsigned int key){return key;}};

template<> struct HASH_REDUCE<float>
{
    static int H(float key)
    {
        STATIC_ASSERT(sizeof(float)==sizeof(int));
        union {float f;int i;} raw;
        raw.f=key;
        return raw.i;
    }
};

template<> struct HASH_REDUCE<double>
{
    static int H(double key)
    {
        STATIC_ASSERT(sizeof(double)==2*sizeof(int));
        union {double d;int i[2];} raw;
        raw.d=key;
        return int_hash(1278312,raw.i[0],raw.i[1]);
    }
};

template<int s> inline int Hash_Reduce_Helper(const void* key);
template<> inline int Hash_Reduce_Helper<1>(const void* key){union {const void* p;int i;} raw;raw.p=key;return raw.i;}
template<> inline int Hash_Reduce_Helper<2>(const void* key){union {const void* p;int i[2];} raw;raw.p=key;return int_hash(raw.i[0],raw.i[1]);}
template<class T> struct HASH_REDUCE<T*>{static int H(const T* key){STATIC_ASSERT((!is_same<T,char>::value && !is_same<T,const char>::value));return Hash_Reduce_Helper<sizeof(T*)/sizeof(int)>(key);}};

//#####################################################################
// Function Hash
//#####################################################################
template<class T> inline int Hash(const T& key){return int_hash(HASH_REDUCE<T>::H(key));}
//#####################################################################
}
#endif
