//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <string>

#include <boost/function.hpp>
#include <boost/spirit/home/qi/char/char_class.hpp>
#include <boost/spirit/home/qi/parse.hpp>

#include <Jeffrey_Utilities/Parsing/TRANSFORMED_LEVEL_SET_GRAMMAR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Parse_String_Level_Set.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
int Parse_String_Level_Set(
    const std::string& level_set_str,
    boost::function< T ( const VECTOR<T,D>& ) >& phi)
{
    namespace qi = boost::spirit::qi;

    static const TRANSFORMED_LEVEL_SET_GRAMMAR<
        T, D,
        std::string::const_iterator,
        qi::standard::blank_type
    > transformed_level_set;

    std::string::const_iterator it = level_set_str.begin();
    const bool b = qi::phrase_parse(
        it, level_set_str.end(),
        transformed_level_set,
        qi::standard::blank,
        phi
    );
    if(b && it == level_set_str.end())
        return 0;

    return 1;
}

#define EXPLICIT_INSTANTIATION( T, D ) \
template int Parse_String_Level_Set( \
    const std::string& level_set_str, \
    boost::function< T ( const VECTOR<T,D>& ) >& phi);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
