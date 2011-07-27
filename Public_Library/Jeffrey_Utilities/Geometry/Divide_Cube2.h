//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_DIVIDE_CUBE2_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_DIVIDE_CUBE2_HPP

#include <cassert>

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/BOUNDED_LIST.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Geometry/CUBE2_SIMPLEX_PARTITION.h>
#include <Jeffrey_Utilities/Geometry/Divide_Simplex.h>
#include <Jeffrey_Utilities/Geometry/Normal.h>
#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Detail_Divide_Cube2
{

template< class T, int D, class T_POLYTOPE_VISITOR >
struct DIVIDE_SIMPLEX_POLYTOPE_VISITOR;

} // namespace Detail_Divide_Cube2

template< class T, int D, class T_PHI_OF_CUBE2_VERTEX, class T_POLYTOPE_VISITOR >
void Divide_Cube2(
    const VECTOR<T,D> dx,
    T_PHI_OF_CUBE2_VERTEX phi_of_cube2_vertex,
    T_POLYTOPE_VISITOR polytope_visitor,
    const int sign_of_zero = -1)
{
    /*************************
     *          2
     * 1 o------o------o 3  o---x
     *
     * 
     *  3       6       9
     *   o------o------o
     *   |      |      |
     *   |      |      |    y
     *   |     5|      |    |
     * 2 o------o------o 8  o---x
     *   |      |      |
     *   |      |      |
     *   |      |      |
     *   o------o------o
     *  1       4       7
     *
     *
     *         3       6       9
     *          o------o------o
     *         /|     /|     /|
     *     12 / | 15 / | 18 / |
     *       o------o------o  |
     *      /|  *--/|--*--/|--o 8
     *  21 / | /| / | /| / | /|
     *    o------o------o  |/ |     z
     *    |  *---|--*---|--o  |     |
     *    | /|  *|-/|--*|-/|--o     *---y
     *    |/ | / |/ | / |/ | / 7   /
     * 20 o------o------o  |/     x
     *    |  *---|--*---|--o
     *    | / 10 | / 13 | / 16
     *    |/     |/     |/
     *    o------o------o
     *  19       22      25
     *************************/

    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef VECTOR< int, D+1 > SIMPLEX_TYPE;

    assert(sign_of_zero == -1 || sign_of_zero == +1);

    BOOST_FOREACH( const SIMPLEX_TYPE simplex, CUBE2_SIMPLEX_PARTITION<D>() ) {
        static const int CENTER_INDEX = CUBE2_SIMPLEX_PARTITION<D>::CENTER_INDEX;
        static_cast<void>(CENTER_INDEX);
        assert(simplex.Contains(CENTER_INDEX));
        static const int MAX_N_VERTEX = D == 1 ? 3 : D == 2 ? 5 : D == 3 ? 9 : 0;
        BOUNDED_LIST< VECTOR<T,D>, MAX_N_VERTEX > vertices_x;
        VECTOR< T, D+1 > phi_of_simplex_vertex;
        VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, D+1 > alignment_of_sub_simplex_opposite_simplex_vertex;
        for(int i = 1; i <= D+1; ++i) {
            const int linear_index = simplex[i];
            const MULTI_INDEX_TYPE multi_index = STATIC_MULTI_INDEX_CUBE<D,-1,+1>::Multi_Index(linear_index);
            vertices_x.Append(multi_index * dx);
            phi_of_simplex_vertex[i] = phi_of_cube2_vertex(linear_index);
            alignment_of_sub_simplex_opposite_simplex_vertex[i] =
                linear_index == CUBE2_SIMPLEX_PARTITION<D>::CENTER_INDEX ?
                SUB_CUBE_POLYTOPE_ALIGNMENT_CUBE_ALIGNED :
                SUB_CUBE_POLYTOPE_ALIGNMENT_SIMPLEX_ALIGNED;
        }
        VECTOR< VECTOR<T,D>, D+1 > normal_of_sub_simplex_opposite_simplex_vertex;
        for(int i = 1; i <= D+1; ++i) {
            VECTOR< VECTOR<T,D>, D > sub_simplex_vertices_x;
            for(int j = 1; j <= D; ++j)
                sub_simplex_vertices_x[j] = vertices_x(j + (i <= j));
            normal_of_sub_simplex_opposite_simplex_vertex[i] = (2 * ((i&1)^(D&1)) - 1) * Normal(sub_simplex_vertices_x);
        }
        typedef Detail_Divide_Cube2::DIVIDE_SIMPLEX_POLYTOPE_VISITOR<
            T, D, T_POLYTOPE_VISITOR
        > DIVIDE_SIMPLEX_POLYTOPE_VISITOR_;
        const ARRAY_VIEW< const VECTOR<T,D> > raw_vertices_x(MAX_N_VERTEX, vertices_x.Data());
        Divide_Simplex(
            vertices_x,
            phi_of_simplex_vertex,
            alignment_of_sub_simplex_opposite_simplex_vertex,
            normal_of_sub_simplex_opposite_simplex_vertex,
            DIVIDE_SIMPLEX_POLYTOPE_VISITOR_(-1, raw_vertices_x, polytope_visitor),
            DIVIDE_SIMPLEX_POLYTOPE_VISITOR_(+1, raw_vertices_x, polytope_visitor),
            sign_of_zero
        );
    }
}

namespace Detail_Divide_Cube2
{

template< class T, int D, class T_POLYTOPE_VISITOR >
struct DIVIDE_SIMPLEX_POLYTOPE_VISITOR
{
    typedef ARRAY_VIEW< const VECTOR<T,D> > VERTICES_X_TYPE;
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        DIVIDE_SIMPLEX_POLYTOPE_VISITOR,
        (( /******/ int const, sign ))
        (( typename VERTICES_X_TYPE const, vertices_x ))
        (( typename T_POLYTOPE_VISITOR const, polytope_visitor ))
    )
public:
    typedef void result_type;
    typedef SUB_CUBE_POLYTOPE<T,D> POLYTOPE_TYPE;
    void operator()(const POLYTOPE_TYPE& polytope) const
    { polytope_visitor(sign, vertices_x, polytope); }
};

} // namespace Detail_Divide_Cube2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_DIVIDE_CUBE2_HPP
