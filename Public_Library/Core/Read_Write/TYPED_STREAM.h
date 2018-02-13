//#####################################################################
// Copyright 2006, Geoffrey Irving, Eftychios Sifakis, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TYPED_ISTREAM and TYPED_OSTREAM
//#####################################################################
// STREAM_TYPE:
//   Indicates whether a stream should be read as float or double.
// 
// TYPED_ISTREAM and TYPED_OSTREAM
//   Wrapper for std::istream and std::ostream that also knows whether
//   it is writing floats or doubles.  Generally should be passed by value.
// 
// FILE_ISTREAM and FILE_OSTREAM
//   Stream that represents a binary file.  When writing a file, creating
//   this object writes the use_doubles flag to the file.  When reading,
//   this flag is read from file.  Can be moved but not copied.  Owns its
//   stream.
// 
#ifndef __TYPED_STREAM__
#define __TYPED_STREAM__

#include <Core/Read_Write/READ_WRITE_FORWARD.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <iosfwd>
#include <memory>
namespace PhysBAM{
//#####################################################################
// Class STREAM_TYPE
//#####################################################################
class STREAM_TYPE
{
public:
    bool use_doubles; // otherwise use floats

    explicit STREAM_TYPE(const bool use_doubles)
        :use_doubles(use_doubles)
    {}

    explicit STREAM_TYPE(const float)
        :use_doubles(false)
    {}

    explicit STREAM_TYPE(const double)
        :use_doubles(true)
    {}
};
//#####################################################################
// Class TYPED_ISTREAM
//#####################################################################
class TYPED_ISTREAM
{
public:
    std::istream& stream;
    STREAM_TYPE type;

    TYPED_ISTREAM(std::istream& stream,STREAM_TYPE type)
        :stream(stream),type(type)
    {}
};
//#####################################################################
// Class TYPED_OSTREAM
//#####################################################################
class TYPED_OSTREAM
{
public:
    std::ostream& stream;
    STREAM_TYPE type;

    TYPED_OSTREAM(std::ostream& stream,STREAM_TYPE type)
        :stream(stream),type(type)
    {}
};
//#####################################################################
// Class FILE_ISTREAM
//#####################################################################
class FILE_ISTREAM
{
public:
    std::istream* stream=0;
    STREAM_TYPE type=STREAM_TYPE(0.f);;

    FILE_ISTREAM(){}

    ~FILE_ISTREAM();

    FILE_ISTREAM(FILE_ISTREAM&& t)
        :stream(t.stream),type(t.type)
    {t.stream=0;}

    FILE_ISTREAM& operator=(FILE_ISTREAM&& t)
    {stream=t.stream;type=t.type;t.stream=0;return *this;}

    operator TYPED_ISTREAM ()
    {return TYPED_ISTREAM(*stream,type);}
    
    FILE_ISTREAM(const FILE_ISTREAM& t)=delete;

    FILE_ISTREAM& operator=(const FILE_ISTREAM& t)=delete;

    void Set(std::istream* stream_input);
};
//#####################################################################
// Class FILE_OSTREAM
//#####################################################################
class FILE_OSTREAM
{
public:
    std::ostream* stream=0;
    STREAM_TYPE type=STREAM_TYPE(0.f);

    FILE_OSTREAM(){}

    ~FILE_OSTREAM();

    FILE_OSTREAM(FILE_OSTREAM&& t)
        :stream(t.stream),type(t.type)
    {t.stream=0;}

    FILE_OSTREAM& operator=(FILE_OSTREAM&& t)
    {stream=t.stream;type=t.type;t.stream=0;return *this;}

    operator TYPED_OSTREAM ()
    {return TYPED_OSTREAM(*stream,type);}

    FILE_OSTREAM(const FILE_OSTREAM& t)=delete;

    FILE_OSTREAM& operator=(const FILE_OSTREAM& t)=delete;

    void Set(std::ostream* stream_input,STREAM_TYPE type_input);
};
//#####################################################################
// Detect whether a type has Read/Write taking typed streams
//#####################################################################
template<class T,class HAS> struct HAS_TYPED_READ{enum {value=false};};
template<class T> struct HAS_TYPED_READ<T,typename FIRST<void,typename T::HAS_TYPED_READ_WRITE>::TYPE>{enum {value=true};};
template<class T,class HAS> struct HAS_TYPED_WRITE{enum {value=false};};
template<class T> struct HAS_TYPED_WRITE<T,typename FIRST<void,typename T::HAS_TYPED_READ_WRITE>::TYPE>{enum {value=true};};
template<class T,class HAS> struct HAS_UNTYPED_READ{enum {value=false};};
template<class T> struct HAS_UNTYPED_READ<T,typename FIRST<void,typename T::HAS_UNTYPED_READ_WRITE>::TYPE>{enum {value=true};};
template<class T,class HAS> struct HAS_UNTYPED_WRITE{enum {value=false};};
template<class T> struct HAS_UNTYPED_WRITE<T,typename FIRST<void,typename T::HAS_UNTYPED_READ_WRITE>::TYPE>{enum {value=true};};
//#####################################################################
}
#endif
