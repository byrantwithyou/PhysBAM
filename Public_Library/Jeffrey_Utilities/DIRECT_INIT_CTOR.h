//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_DIRECT_INIT_CTOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_DIRECT_INIT_CTOR_HPP

#include <boost/call_traits.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/control/expr_iif.hpp>
#include <boost/preprocessor/detail/is_nullary.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/facilities/expand.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/tuple/elem.hpp>

#define PHYSBAM_HAS_KEYWORD_TYPENAME_PREFIX_typename ()
#define PHYSBAM_HAS_KEYWORD_TYPENAME_PREFIX( Type ) \
    BOOST_PP_IS_NULLARY( BOOST_PP_CAT( PHYSBAM_HAS_KEYWORD_TYPENAME_PREFIX_, Type ) )
#define PHYSBAM_REMOVE_KEYWORD_TYPENAME_PREFIX( Type ) \
    BOOST_PP_CAT( \
        PHYSBAM_REMOVE_KEYWORD_TYPENAME_PREFIX_dispatch, \
        PHYSBAM_HAS_KEYWORD_TYPENAME_PREFIX( Type ) \
    ) ( Type )
#define PHYSBAM_REMOVE_KEYWORD_TYPENAME_PREFIX_dispatch0( Type ) Type
#define PHYSBAM_REMOVE_KEYWORD_TYPENAME_PREFIX_dispatch1( Type ) \
    BOOST_PP_EXPAND( BOOST_PP_EMPTY BOOST_PP_CAT( PHYSBAM_HAS_KEYWORD_TYPENAME_PREFIX_, Type ) )

#define PHYSBAM_DIRECT_INIT_CTOR_comma_ctor_param( r, data, i, elem ) \
    BOOST_PP_COMMA_IF( i ) \
    BOOST_PP_EXPR_IIF( PHYSBAM_HAS_KEYWORD_TYPENAME_PREFIX( BOOST_PP_TUPLE_ELEM( 2, 0, elem ) ), typename ) \
    ::boost::call_traits< PHYSBAM_REMOVE_KEYWORD_TYPENAME_PREFIX( BOOST_PP_TUPLE_ELEM( 2, 0, elem ) ) >::param_type \
    BOOST_PP_CAT( _ ## i ## _, BOOST_PP_TUPLE_ELEM( 2, 1, elem ) )

#define PHYSBAM_DIRECT_INIT_CTOR_comma_ctor_init( r, data, i, elem ) \
    BOOST_PP_COMMA_IF( i ) \
    BOOST_PP_TUPLE_ELEM( 2, 1, elem ) ( \
        BOOST_PP_CAT( _ ## i ## _, BOOST_PP_TUPLE_ELEM( 2, 1, elem ) ) \
    )

#define PHYSBAM_DIRECT_INIT_CTOR_declare_member( r, data, i, elem ) \
    PHYSBAM_REMOVE_KEYWORD_TYPENAME_PREFIX( BOOST_PP_TUPLE_ELEM( 2, 0, elem ) ) \
    BOOST_PP_TUPLE_ELEM( 2, 1, elem ) ;
#define PHYSBAM_DIRECT_INIT_CTOR_declare_members( MemberSeq ) \
    BOOST_PP_SEQ_FOR_EACH_I( PHYSBAM_DIRECT_INIT_CTOR_declare_member, ~, MemberSeq )

#define PHYSBAM_DIRECT_INIT_CTOR( Class, MemberSeq ) \
    BOOST_PP_EXPR_IIF( BOOST_PP_EQUAL( BOOST_PP_SEQ_SIZE( MemberSeq ), 1 ), explicit ) \
    Class ( BOOST_PP_SEQ_FOR_EACH_I( PHYSBAM_DIRECT_INIT_CTOR_comma_ctor_param, ~, MemberSeq ) ) \
        : BOOST_PP_SEQ_FOR_EACH_I( PHYSBAM_DIRECT_INIT_CTOR_comma_ctor_init, ~, MemberSeq ) \
    { }

#define PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS( Class, MemberSeq ) \
    PHYSBAM_DIRECT_INIT_CTOR( Class, MemberSeq ) \
    PHYSBAM_DIRECT_INIT_CTOR_declare_members( MemberSeq )

#define PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PUBLIC_MEMBERS( Class, MemberSeq ) \
    PHYSBAM_DIRECT_INIT_CTOR( Class, MemberSeq ) \
    public: \
    PHYSBAM_DIRECT_INIT_CTOR_declare_members( MemberSeq )

#define PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PROTECTED_MEMBERS( Class, MemberSeq ) \
    PHYSBAM_DIRECT_INIT_CTOR( Class, MemberSeq ) \
    protected: \
    PHYSBAM_DIRECT_INIT_CTOR_declare_members( MemberSeq )

#define PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS( Class, MemberSeq ) \
    PHYSBAM_DIRECT_INIT_CTOR( Class, MemberSeq ) \
    private: \
    PHYSBAM_DIRECT_INIT_CTOR_declare_members( MemberSeq )

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_DIRECT_INIT_CTOR_HPP
