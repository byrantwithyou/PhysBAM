//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INTERFACE_INDEX_TRANSFORM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INTERFACE_INDEX_TRANSFORM_HPP

#include <cassert>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
struct INTERFACE_INDEX_TRANSFORM
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        INTERFACE_INDEX_TRANSFORM,
        (( /******/ int const, n_grid_index ))
        (( typename T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX const, virtual_index_offset_of_grid_index ))
        (( typename T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET const, grid_index_of_virtual_index_offset ))
        (( typename T_SIGN_OF_GRID_INDEX const, sign_of_grid_index ))
    )

    int Index_Of_Material_Grid_Index(const int grid_index) const;
    int Index_Of_Virtual_Grid_Index(const int grid_index) const;
    int Index_Of_Signed_Grid_Index(const int grid_index, const int sign) const;
    bool Index_Is_Material(const int index) const;
    bool Index_Is_Virtual(const int index) const;
    int Grid_Index_Of_Index(const int index) const;
    int Domain_Sign_Of_Index(const int index) const;

    class INDEX_OF_SIGNED_GRID_INDEX_FUNCTION;
    INDEX_OF_SIGNED_GRID_INDEX_FUNCTION
    Index_Of_Signed_Grid_Index_Function(const int sign) const;
};

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
inline INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>
Make_Interface_Index_Transform(
    const int n_grid_index,
    const T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX& virtual_index_offset_of_grid_index,
    const T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET& grid_index_of_virtual_index_offset,
    const T_SIGN_OF_GRID_INDEX& sign_of_grid_index)
{
    return INTERFACE_INDEX_TRANSFORM<
        T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
        T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
        T_SIGN_OF_GRID_INDEX
    >(
        n_grid_index,
        virtual_index_offset_of_grid_index,
        grid_index_of_virtual_index_offset,
        sign_of_grid_index
    );
}

//#####################################################################
//#####################################################################

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
inline int
INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::
Index_Of_Material_Grid_Index(const int grid_index) const
{
    assert(1 <= grid_index && grid_index <= n_grid_index);
    return grid_index;
}

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
inline int
INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::
Index_Of_Virtual_Grid_Index(const int grid_index) const
{
    assert(1 <= grid_index && grid_index <= n_grid_index);
    return n_grid_index + virtual_index_offset_of_grid_index(grid_index);
}

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
inline int
INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::
Index_Of_Signed_Grid_Index(const int grid_index, const int sign) const
{
    return sign_of_grid_index(grid_index) == sign ?
           Index_Of_Material_Grid_Index(grid_index) :
           Index_Of_Virtual_Grid_Index(grid_index);
}

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
inline bool
INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::
Index_Is_Material(const int index) const
{ return index <= n_grid_index; }

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
inline bool
INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::
Index_Is_Virtual(const int index) const
{ return index > n_grid_index; }

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
inline int
INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::
Grid_Index_Of_Index(const int index) const
{ return Index_Is_Material(index) ? index : grid_index_of_virtual_index_offset(index - n_grid_index); }

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
inline int
INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::
Domain_Sign_Of_Index(const int index) const
{
    if(Index_Is_Material(index))
        return sign_of_grid_index(index);
    const int grid_index = grid_index_of_virtual_index_offset(index - n_grid_index);
    return -sign_of_grid_index(grid_index);
}

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
class INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::INDEX_OF_SIGNED_GRID_INDEX_FUNCTION
{
    friend struct INTERFACE_INDEX_TRANSFORM;
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        INDEX_OF_SIGNED_GRID_INDEX_FUNCTION,
        (( typename INTERFACE_INDEX_TRANSFORM const &, this_ ))
        (( /******/ int const, sign ))
    )
public:
    typedef int result_type;
    int operator()(const int grid_index) const
    { return this_.Index_Of_Signed_Grid_Index(grid_index, sign); }
};

template<
    class T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    class T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    class T_SIGN_OF_GRID_INDEX
>
inline typename INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::INDEX_OF_SIGNED_GRID_INDEX_FUNCTION
INTERFACE_INDEX_TRANSFORM<
    T_VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX,
    T_GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET,
    T_SIGN_OF_GRID_INDEX
>::
Index_Of_Signed_Grid_Index_Function(const int sign) const
{ return INDEX_OF_SIGNED_GRID_INDEX_FUNCTION(*this, sign); }

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_INTERFACE_INDEX_TRANSFORM_HPP
