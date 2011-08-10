//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ONSTREAM_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ONSTREAM_HPP

#include <ostream>
#include <streambuf>

namespace PhysBAM
{

template< class T_CHAR, class T_TRAITS = std::char_traits< T_CHAR > >
struct BASIC_ONSTREAM
    : std::basic_ostream< T_CHAR, T_TRAITS >
{
    BASIC_ONSTREAM()
        : std::basic_ostream< T_CHAR, T_TRAITS >(&m_streambuf)
    { }

private:
    struct STREAMBUF
        : std::basic_streambuf< T_CHAR, T_TRAITS >
    {
        typename T_TRAITS::int_type overflow(typename T_TRAITS::int_type c)
        { return T_TRAITS::not_eof(c); }
    } m_streambuf;
};

typedef BASIC_ONSTREAM<  char   >  ONSTREAM;
typedef BASIC_ONSTREAM< wchar_t > WONSTREAM;

extern  ONSTREAM  nout;
extern WONSTREAM wnout;

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_ONSTREAM_HPP
