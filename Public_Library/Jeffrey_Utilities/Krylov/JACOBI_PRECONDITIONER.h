//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_JACOBI_PRECONDITIONER_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_JACOBI_PRECONDITIONER_HPP

#include <cassert>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/RESULT_OF.h>

namespace PhysBAM
{

template< class T_DIAG_OF_INDEX >
struct JACOBI_PRECONDITIONER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        JACOBI_PRECONDITIONER,
        (( typename T_DIAG_OF_INDEX const, diag_of_index ))
    )
public:
    typedef void result_type;
    template< class T_VECTOR >
    void operator()(const T_VECTOR& r, T_VECTOR& z) const
    {
        assert(r.Size() == z.Size());
        for(int i = 1; i <= r.Size(); ++i)
            z(i) = r(i) / diag_of_index(i);
    }
};

template< class T_DIAG_OF_INDEX >
struct SAFE_JACOBI_PRECONDITIONER
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SAFE_JACOBI_PRECONDITIONER,
        (( typename T_DIAG_OF_INDEX const, diag_of_index ))
    )
public:
    typedef void result_type;
    template< class T_VECTOR >
    void operator()(const T_VECTOR& r, T_VECTOR& z) const
    {
        typedef typename RESULT_OF< const T_DIAG_OF_INDEX ( int ) >::type DIAG_TYPE;
        assert(r.Size() == z.Size());
        for(int i = 1; i <= r.Size(); ++i) {
            const DIAG_TYPE diag = diag_of_index(i);
            z(i) = diag == 0 ? r(i) : r(i) / diag;
        }
    }
};

template< class T_DIAG_OF_INDEX >
inline JACOBI_PRECONDITIONER< T_DIAG_OF_INDEX >
Make_Jacobi_Preconditioner(const T_DIAG_OF_INDEX& diag_of_index)
{ return JACOBI_PRECONDITIONER< T_DIAG_OF_INDEX >(diag_of_index); }

template< class T_DIAG_OF_INDEX >
inline SAFE_JACOBI_PRECONDITIONER< T_DIAG_OF_INDEX >
Make_Safe_Jacobi_Preconditioner(const T_DIAG_OF_INDEX& diag_of_index)
{ return SAFE_JACOBI_PRECONDITIONER< T_DIAG_OF_INDEX >(diag_of_index); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_KRYLOV_JACOBI_PRECONDITIONER_HPP
