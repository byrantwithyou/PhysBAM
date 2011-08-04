//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <string>

#include <boost/spirit/home/qi/char/char_class.hpp>

#include <Jeffrey_Utilities/Parsing/TRANSFORMED_LEVEL_SET_GRAMMAR.ipp>

namespace PhysBAM
{

#define EXPLICIT_INSTANTIATION( T, D ) \
template struct TRANSFORMED_LEVEL_SET_GRAMMAR< \
    T, D, \
    std::string::const_iterator, \
    boost::spirit::qi::standard::blank_type \
>;
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( double, 3 )
#endif // #ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#undef EXPLICIT_INSTANTIATION

} // namespace PhysBAM
