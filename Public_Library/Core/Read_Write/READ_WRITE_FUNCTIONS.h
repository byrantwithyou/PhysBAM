//#####################################################################
// Copyright 2004-2006, Eran Guendelman, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FUNCTIONS
//#####################################################################
// Functions for reading and writing which do the correct thing for objects, pointers, primitive types, etc. In general, use Read/Write_Binary (and Read/Write_Binary_Array) using T for the type
// of the object you're reading/writing and RW the underlying floating point scalar type (float/double).
//#####################################################################
#ifndef __READ_WRITE_FUNCTIONS__
#define __READ_WRITE_FUNCTIONS__

#include <Core/Read_Write/READ_WRITE_FORWARD.h>
#include <Core/Read_Write/TYPED_STREAM.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <complex>
namespace PhysBAM{

template<class R,class A>
inline std::ostream& operator<<(std::ostream& output_stream,R (*func)(A))
{output_stream<<(void*)func;return output_stream;}
template<class R,class A,class B>
inline std::ostream& operator<<(std::ostream& output_stream,R (*func)(A,B))
{output_stream<<(void*)func;return output_stream;}
template<class R,class A,class B,class C>
inline std::ostream& operator<<(std::ostream& output_stream,R (*func)(A,B,C))
{output_stream<<(void*)func;return output_stream;}
template<class R,class A,class B,class C,class D>
inline std::ostream& operator<<(std::ostream& output_stream,R (*func)(A,B,C,D))
{output_stream<<(void*)func;return output_stream;}
//#####################################################################
// Read_Binary
//#####################################################################
template<class RW,class T> inline typename enable_if<!HAS_TYPED_READ<T>::value && !IS_BINARY_IO_SAFE<T,RW>::value>::type
Read_Binary(std::istream& input,T& d)
{d.template Read<RW>(input);}

template<class RW,class T> inline typename enable_if<IS_BINARY_IO_SAFE<T,RW>::value>::type
Read_Binary(std::istream& input,T& d)
{input.read(reinterpret_cast<char*>(&d),sizeof(T));}

template<class RW,class T> inline typename enable_if<HAS_TYPED_READ<T>::value && !IS_BINARY_IO_SAFE<T,RW>::value>::type
Read_Binary(std::istream& input,T& d)
{TYPED_ISTREAM typed_input(input,STREAM_TYPE(RW()));d.Read(typed_input);}

template<class T> inline void Read_Primitive(std::istream& input,T& d);

template<class RW> inline void
Read_Binary(std::istream& input,float& d)
{RW tmp;Read_Primitive(input,tmp);d=(float)tmp;}

template<class RW> inline void
Read_Binary(std::istream& input,double& d)
{RW tmp;Read_Primitive(input,tmp);d=(double)tmp;}

template<class RW,class T> inline void
Read_Binary(std::istream& input,std::complex<T>& d)
{RW r,c;Read_Primitive(input,r);Read_Primitive(input,c);d=std::complex<T>(r,c);}

template<class RW,class T> inline void
Read_Binary(std::istream& input,T*& d)
{bool data_exists;Read_Primitive(input,data_exists);if(data_exists){d=new T();Read_Binary<RW>(input,*d);}else d=0;}

template<class RW> inline void
Read_Binary(std::istream& input,std::string& d)
{int n;Read_Primitive(input,n);char* buffer=new char[n];input.read(buffer,n);d.assign(buffer,buffer+n);delete[] buffer;}

template<class T> inline typename enable_if<!HAS_TYPED_READ<T>::value>::type
Read_Binary(TYPED_ISTREAM& input,T& d)
{
    if(input.type.use_doubles)
        Read_Binary<double>(input.stream,d);
    else
        Read_Binary<float>(input.stream,d);
}

template<class T> inline typename enable_if<HAS_TYPED_READ<T>::value>::type
Read_Binary(TYPED_ISTREAM& input,T& d)
{d.Read(input);}

//#####################################################################
// Write_Binary
//#####################################################################
template<class RW,class T> inline typename enable_if<!HAS_TYPED_WRITE<T>::value && !IS_BINARY_IO_SAFE<T,RW>::value && !is_pointer<T>::value>::type
Write_Binary(std::ostream& output,const T& d)
{d.template Write<RW>(output);}

template<class RW,class T> inline typename enable_if<IS_BINARY_IO_SAFE<T,RW>::value>::type
Write_Binary(std::ostream& output,const T& d)
{output.write(reinterpret_cast<const char*>(&d),sizeof(T));}

template<class RW,class T> inline typename enable_if<HAS_TYPED_WRITE<T>::value && !IS_BINARY_IO_SAFE<T,RW>::value>::type
Write_Binary(std::ostream& output,const T& d)
{TYPED_OSTREAM typed_output(output,STREAM_TYPE(RW()));d.Write(typed_output);}

template<class T> inline void Write_Primitive(std::ostream& output,const T& d);

template<class RW> inline void
Write_Binary(std::ostream& output,const float& d)
{Write_Primitive(output,(RW)d);}

template<class RW> inline void
Write_Binary(std::ostream& output,const double& d)
{Write_Primitive(output,(RW)d);}

template<class RW,class T> inline void
Write_Binary(std::ostream& output,const std::complex<T>& d)
{Write_Primitive(output,(RW)d.real());Write_Primitive(output,(RW)d.imag());}

template<class RW,class T> inline void
Write_Binary(std::ostream& output,const T* d)
{Write_Primitive(output,d!=0);if(d) Write_Binary<RW>(output,*d);} // write a bool tag indicating whether pointer's data follows

template<class RW> inline void
Write_Binary(std::ostream& output,const std::string& d)
{int n=int(d.size());Write_Primitive(output,n);const char* s=d.c_str();output.write(s,n);}

template<class T> inline typename enable_if<!HAS_TYPED_WRITE<T>::value>::type
Write_Binary(TYPED_OSTREAM& output,const T& d)
{
    if(output.type.use_doubles)
        Write_Binary<double>(output.stream,d);
    else
        Write_Binary<float>(output.stream,d);
}

template<class T> inline typename enable_if<HAS_TYPED_WRITE<T>::value>::type
Write_Binary(TYPED_OSTREAM& output,const T& d)
{
    d.Write(output);
}

//#####################################################################
// Multiple Argument Read_Binary
//#####################################################################
template<class RW,class T1,class T2,class ...Args>
inline void Read_Binary(std::istream& input,T1& d1,T2& d2,Args&& ...args)
{Read_Binary<RW>(input,d1);Read_Binary<RW>(input,d2,args...);}

template<class T1,class T2,class ...Args>
inline void Read_Binary(TYPED_ISTREAM& input,T1& d1,T2& d2,Args&& ...args)
{Read_Binary(input,d1);Read_Binary(input,d2,args...);}

//#####################################################################
// Multiple Argument Write_Binary
//#####################################################################
template<class RW,class T1,class T2,class ...Args>
inline void Write_Binary(std::ostream& output,const T1& d1,const T2& d2,Args&& ...args)
{Write_Binary<RW>(output,d1);Write_Binary<RW>(output,d2,args...);}

template<class T1,class T2,class ...Args>
inline void Write_Binary(TYPED_OSTREAM& output,const T1& d1,const T2& d2,Args&& ...args)
{Write_Binary(output,d1);Write_Binary(output,d2,args...);}

//#####################################################################
// Read/Write_Binary_Array
//#####################################################################
template<class RW,class T> inline typename enable_if<IS_BINARY_IO_SAFE<T,RW>::value>::type
Read_Binary_Array(std::istream& input,T* array,const int number_of_elements)
{input.read(reinterpret_cast<char*>(array),number_of_elements*sizeof(T));}

template<class RW,class T> inline typename enable_if<!IS_BINARY_IO_SAFE<T,RW>::value>::type
Read_Binary_Array(std::istream& input,T* array,const int number_of_elements)
{for(int i=0;i<number_of_elements;i++) Read_Binary<RW>(input,array[i]);}

template<class RW,class T> inline typename enable_if<IS_BINARY_IO_SAFE<T,RW>::value>::type
Write_Binary_Array(std::ostream& output,const T* array,const int number_of_elements)
{if(number_of_elements) output.write(reinterpret_cast<const char*>(array),number_of_elements*sizeof(T));}

template<class RW,class T> inline typename enable_if<!IS_BINARY_IO_SAFE<T,RW>::value>::type
Write_Binary_Array(std::ostream& output,const T* array,const int number_of_elements)
{for(int i=0;i<number_of_elements;i++) Write_Binary<RW>(output,array[i]);}
//#####################################################################
}
#include <Core/Read_Write/READ_WRITE.h>
#endif
