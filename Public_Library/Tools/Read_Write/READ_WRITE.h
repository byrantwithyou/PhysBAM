//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class Read_Write
//#####################################################################
// Functions for reading and writing which do the correct thing for objects, pointers, primitive types, etc. In general, use Read/Write_Binary (and Read/Write_Binary_Array) using T for the type
// of the object you're reading/writing and RW the underlying floating point scalar type (float/double).
//#####################################################################
#ifndef __Read_Write__
#define __Read_Write__

#include <Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Tools/Utilities/STATIC_ASSERT.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/SCALAR_POLICY.h>
#include <cassert>
#include <iostream>
namespace PhysBAM{

#ifdef PHYSBAM_BIG_ENDIAN
static const bool big_endian=true;
#else
static const bool big_endian=false;
#endif

template<> struct PLATFORM_INDEPENDENT_SIZE<bool>{static const int value=1;};
template<> struct PLATFORM_INDEPENDENT_SIZE<char>{static const int value=1;};
template<> struct PLATFORM_INDEPENDENT_SIZE<unsigned char>{static const int value=1;};
template<> struct PLATFORM_INDEPENDENT_SIZE<short>{static const int value=2;};
template<> struct PLATFORM_INDEPENDENT_SIZE<unsigned short>{static const int value=2;};
template<> struct PLATFORM_INDEPENDENT_SIZE<int>{static const int value=4;};
template<> struct PLATFORM_INDEPENDENT_SIZE<unsigned int>{static const int value=4;};
template<> struct PLATFORM_INDEPENDENT_SIZE<float>{static const int value=4;};
template<> struct PLATFORM_INDEPENDENT_SIZE<double>{static const int value=8;};
template<class T> struct PLATFORM_INDEPENDENT_SIZE<T,typename enable_if<is_enum<T>::value>::type>{static const int value=4;};

template<class T> struct IS_PRIMITIVE_BINARY_IO_SAFE {static const bool value=!big_endian && (sizeof(T)==PLATFORM_INDEPENDENT_SIZE<T>::value);};

// classify types which can be written directly to files without conversion
template<class T,class RW,class ENABLER> struct IS_BINARY_IO_SAFE {static const bool value=false;};
template<class T,class RW> struct IS_BINARY_IO_SAFE<T,RW,typename enable_if<is_integral<T>::value || is_enum<T>::value>::type>:public IS_PRIMITIVE_BINARY_IO_SAFE<T>{};
template<> struct IS_BINARY_IO_SAFE<float,float>:public IS_PRIMITIVE_BINARY_IO_SAFE<float>{};
template<> struct IS_BINARY_IO_SAFE<double,double>:public IS_PRIMITIVE_BINARY_IO_SAFE<double>{};

//#####################################################################
// Function Swap_Endianity
//#####################################################################
template<class T>
inline void Swap_Endianity(T& x)
{assert(sizeof(T)<=8);
if(sizeof(T)>1) {T old=x;for(unsigned int k=0;k<sizeof(T);k++) ((char*)&x)[k]=((char*)&old)[sizeof(T)-k-1];}}

//#####################################################################
// Read and Write for primitives
//#####################################################################
template<class T> inline void
Read_Primitive(std::istream& input,T& d)
{
    STATIC_ASSERT((sizeof(T)==PLATFORM_INDEPENDENT_SIZE<T>::value));
    input.read((char*)&d,sizeof(T));
    if(big_endian) Swap_Endianity(d); // convert to little endian if necessary
}

template<class T> inline void
Write_Primitive(std::ostream& output,const T& d)
{
    STATIC_ASSERT((sizeof(T)==PLATFORM_INDEPENDENT_SIZE<T>::value));
    if(big_endian){T d2=d;Swap_Endianity(d2);output.write((const char*)&d2,sizeof(T));} // convert to big endian if necessary
    else output.write((const char*)&d,sizeof(T));
}
}
#endif
