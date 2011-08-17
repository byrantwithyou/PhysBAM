//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_NULL_FUNCTION_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_NULL_FUNCTION_HPP

namespace PhysBAM
{

struct NULL_FUNCTION
{
    typedef void result_type;
    template< class T1 >
    void operator()(const T1&) const
    { }
    template< class T1, class T2 >
    void operator()(const T1&, const T2&) const
    { }
    template< class T1, class T2, class T3 >
    void operator()(const T1&, const T2&, const T3&) const
    { }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_FUNCTIONAL_NULL_FUNCTION_HPP
