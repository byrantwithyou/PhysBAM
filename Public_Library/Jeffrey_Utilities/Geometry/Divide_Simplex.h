//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_DIVIDE_SIMPLEX_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_DIVIDE_SIMPLEX_HPP

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/BOUNDED_LIST.h>
#include <Jeffrey_Utilities/Geometry/Normal.h>
#include <Jeffrey_Utilities/Geometry/SUB_CUBE_POLYTOPE.h>

namespace PhysBAM
{

namespace Detail_Divide_Simplex
{

template< class T, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_0_3(
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 2 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,1>, 2 >& normal_of_sub_simplex_opposite_vertex,
    T_POLYTOPE_VISITOR polytope_visitor);

template< class T, int N, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_1_2(
    BOUNDED_LIST< VECTOR<T,1>, N >& vertices_x,
    const VECTOR<T,2>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 2 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,1>, 2 >& normal_of_sub_simplex_opposite_vertex,
    const int u, const int v,
    T_POLYTOPE_VISITOR u_polytope_visitor,
    T_POLYTOPE_VISITOR v_polytope_visitor);

template< class T, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_0_7(
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 3 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,2>, 3 >& normal_of_sub_simplex_opposite_vertex,
    T_POLYTOPE_VISITOR polytope_visitor);

template< class T, int N, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_1_2_3_4_5_6(
    BOUNDED_LIST< VECTOR<T,2>, N >& vertices_x,
    const VECTOR<T,3>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 3 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,2>, 3 >& normal_of_sub_simplex_opposite_vertex,
    const int u, const int v1, const int v2,
    T_POLYTOPE_VISITOR u_polytope_visitor,
    T_POLYTOPE_VISITOR v_polytope_visitor);

template< class T, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_0_15(
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 4 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,3>, 4 >& normal_of_sub_simplex_opposite_vertex,
    T_POLYTOPE_VISITOR polytope_visitor);

template< class T, int N, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_1_2_4_7_8_11_13_14(
    BOUNDED_LIST< VECTOR<T,3>, N >& vertices_x,
    const VECTOR<T,4>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 4 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,3>, 4 >& normal_of_sub_simplex_opposite_vertex,
    const int u, const int v1, const int v2, const int v3,
    T_POLYTOPE_VISITOR u_polytope_visitor,
    T_POLYTOPE_VISITOR v_polytope_visitor);

template< class T, int N, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_3_5_6_9_10_12(
    BOUNDED_LIST< VECTOR<T,3>, N >& vertices_x,
    const VECTOR<T,4>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 4 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,3>, 4 >& normal_of_sub_simplex_opposite_vertex,
    const int u1, const int u2, const int v1, const int v2,
    T_POLYTOPE_VISITOR u_polytope_visitor,
    T_POLYTOPE_VISITOR v_polytope_visitor);

} // namespace Detail_Divide_Simplex

template< class T, int N, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex(
    BOUNDED_LIST< VECTOR<T,1>, N >& vertices_x,
    const VECTOR<T,2>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 2 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,1>, 2 >& normal_of_sub_simplex_opposite_vertex,
    T_POLYTOPE_VISITOR negative_polytope_visitor,
    T_POLYTOPE_VISITOR positive_polytope_visitor,
    const int sign_of_zero = -1)
{
    BOOST_MPL_ASSERT_RELATION( N, >=, 3 );
    assert(vertices_x.Size() == 2);
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    /********
     * 2---1
     ********/
    const unsigned int case_ =
        ((phi_of_vertex[1] > 0 || (phi_of_vertex[1] == 0 && sign_of_zero == +1)) << 1)
      | ((phi_of_vertex[2] > 0 || (phi_of_vertex[2] == 0 && sign_of_zero == +1)) << 0);
    static const int _topological_case[] = { 0, 1, 1, 0 };
    const unsigned int topological_case = _topological_case[case_];
    if(topological_case == 1) {
        static const int _su[] = { 0,-1,+1, 0 };
        static const int _sv[] = { 0,+1,-1, 0 };
        static const int u = 1;
        const int su = _su[case_];
        static const int v = 2;
        const int sv = _sv[case_];
        assert(su + sv == 0);
        Detail_Divide_Simplex::Divide_Simplex_Cases_1_2(
            vertices_x, phi_of_vertex,
            alignment_of_sub_simplex_opposite_vertex,
            normal_of_sub_simplex_opposite_vertex,
            u, v,
            su == -1 ? negative_polytope_visitor : positive_polytope_visitor,
            sv == -1 ? negative_polytope_visitor : positive_polytope_visitor
        );
    }
    else {
        assert(topological_case == 0);
        Detail_Divide_Simplex::Divide_Simplex_Cases_0_3(
            alignment_of_sub_simplex_opposite_vertex,
            normal_of_sub_simplex_opposite_vertex,
            case_ == 0 ? negative_polytope_visitor : positive_polytope_visitor
        );
    }
}

template< class T, int N, class T_POLYTOPE_VISITOR >
void
Divide_Simplex(
    BOUNDED_LIST< VECTOR<T,2>, N >& vertices_x,
    const VECTOR<T,3>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 3 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,2>, 3 >& normal_of_sub_simplex_opposite_vertex,
    T_POLYTOPE_VISITOR negative_polytope_visitor,
    T_POLYTOPE_VISITOR positive_polytope_visitor,
    const int sign_of_zero = -1)
{
    BOOST_MPL_ASSERT_RELATION( N, >=, 5 );
    assert(vertices_x.Size() == 3);
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    /********
     *   1
     *  / \
     * 2---3
     ********/
    const unsigned int case_ =
        ((phi_of_vertex[1] > 0 || (phi_of_vertex[1] == 0 && sign_of_zero == +1)) << 2)
      | ((phi_of_vertex[2] > 0 || (phi_of_vertex[2] == 0 && sign_of_zero == +1)) << 1)
      | ((phi_of_vertex[3] > 0 || (phi_of_vertex[3] == 0 && sign_of_zero == +1)) << 0);
    static const unsigned int _topological_case[] = { 0, 1, 1, 1, 1, 1, 1, 0 };
    const unsigned int topological_case = _topological_case[case_];
    if(topological_case == 1) {
        static const int _u [] = { 0, 3, 2, 1, 1, 2, 3, 0 };
        static const int _su[] = { 0,+1,+1,-1,+1,-1,-1, 0 };
        static const int _v1[] = { 0, 1, 3, 2, 2, 3, 1, 0 };
        static const int _v2[] = { 0, 2, 1, 3, 3, 1, 2, 0 };
        static const int _sv[] = { 0,-1,-1,+1,-1,+1,+1, 0 };
        const int u  = _u [case_];
        const int su = _su[case_];
        const int v1 = _v1[case_];
        const int v2 = _v2[case_];
        const int sv = _sv[case_];
        assert(su + sv == 0);
        Detail_Divide_Simplex::Divide_Simplex_Cases_1_2_3_4_5_6(
            vertices_x, phi_of_vertex,
            alignment_of_sub_simplex_opposite_vertex,
            normal_of_sub_simplex_opposite_vertex,
            u, v1, v2,
            su == -1 ? negative_polytope_visitor : positive_polytope_visitor,
            sv == -1 ? negative_polytope_visitor : positive_polytope_visitor
        );
    }
    else {
        assert(topological_case == 0);
        Detail_Divide_Simplex::Divide_Simplex_Cases_0_7(
            alignment_of_sub_simplex_opposite_vertex,
            normal_of_sub_simplex_opposite_vertex,
            case_ == 0 ? negative_polytope_visitor : positive_polytope_visitor
        );
    }
}

template< class T, int N, class T_POLYTOPE_VISITOR >
void
Divide_Simplex(
    BOUNDED_LIST< VECTOR<T,3>, N >& vertices_x,
    const VECTOR<T,4>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 4 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,3>, 4 >& normal_of_sub_simplex_opposite_vertex,
    T_POLYTOPE_VISITOR negative_polytope_visitor,
    T_POLYTOPE_VISITOR positive_polytope_visitor,
    const int sign_of_zero = -1)
{
    BOOST_MPL_ASSERT_RELATION( N, >=, 9 );
    assert(vertices_x.Size() == 4);
    assert(sign_of_zero == -1 || sign_of_zero == +1);
    /**********
     *    1
     *   /|\
     *  / | \
     * 2--|--4
     *  \ | /
     *   \|/
     *    3
     **********/
    const unsigned int case_ =
        ((phi_of_vertex[1] > 0 || (phi_of_vertex[1] == 0 && sign_of_zero == +1)) << 3)
      | ((phi_of_vertex[2] > 0 || (phi_of_vertex[2] == 0 && sign_of_zero == +1)) << 2)
      | ((phi_of_vertex[3] > 0 || (phi_of_vertex[3] == 0 && sign_of_zero == +1)) << 1)
      | ((phi_of_vertex[4] > 0 || (phi_of_vertex[4] == 0 && sign_of_zero == +1)) << 0);
    static const unsigned int _topological_case[] = { 0, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 0 };
    const unsigned int topological_case = _topological_case[case_];
    if(topological_case == 1) {
        static const int _u [] = { 0, 4, 3, 0, 2, 0, 0, 1, 1, 0, 0, 2, 0, 3, 4, 0 };
        static const int _su[] = { 0,+1,+1, 0,+1, 0, 0,-1,+1, 0, 0,-1, 0,-1,-1, 0 };
        static const int _v1[] = { 0, 1, 1, 0, 1, 0, 0, 2, 2, 0, 0, 1, 0, 1, 1, 0 };
        static const int _v2[] = { 0, 3, 2, 0, 4, 0, 0, 3, 3, 0, 0, 4, 0, 2, 3, 0 };
        static const int _v3[] = { 0, 2, 4, 0, 3, 0, 0, 4, 4, 0, 0, 3, 0, 4, 2, 0 };
        static const int _sv[] = { 0,-1,-1, 0,-1, 0, 0,+1,-1, 0, 0,+1, 0,+1,+1, 0 };
        const int u  = _u [case_];
        const int su = _su[case_];
        const int v1 = _v1[case_];
        const int v2 = _v2[case_];
        const int v3 = _v3[case_];
        const int sv = _sv[case_];
        assert(su + sv == 0);
        Detail_Divide_Simplex::Divide_Simplex_Cases_1_2_4_7_8_11_13_14(
            vertices_x, phi_of_vertex,
            alignment_of_sub_simplex_opposite_vertex,
            normal_of_sub_simplex_opposite_vertex,
            u, v1, v2, v3,
            su == -1 ? negative_polytope_visitor : positive_polytope_visitor,
            sv == -1 ? negative_polytope_visitor : positive_polytope_visitor
        );
    }
    else if(topological_case == 2) {
        static const int _u1[] = { 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0 };
        static const int _u2[] = { 0, 0, 0, 2, 0, 3, 4, 0, 0, 4, 3, 0, 2, 0, 0, 0 };
        static const int _su[] = { 0, 0, 0,-1, 0,-1,-1, 0, 0,+1,+1, 0,+1, 0, 0, 0 };
        static const int _v1[] = { 0, 0, 0, 3, 0, 4, 2, 0, 0, 2, 4, 0, 3, 0, 0, 0 };
        static const int _v2[] = { 0, 0, 0, 4, 0, 2, 3, 0, 0, 3, 2, 0, 4, 0, 0, 0 };
        static const int _sv[] = { 0, 0, 0,+1, 0,+1,+1, 0, 0,-1,-1, 0,-1, 0, 0, 0 };
        const int u1 = _u1[case_];
        const int u2 = _u2[case_];
        const int su = _su[case_];
        const int v1 = _v1[case_];
        const int v2 = _v2[case_];
        const int sv = _sv[case_];
        assert(su + sv == 0);
        Detail_Divide_Simplex::Divide_Simplex_Cases_3_5_6_9_10_12(
            vertices_x, phi_of_vertex,
            alignment_of_sub_simplex_opposite_vertex,
            normal_of_sub_simplex_opposite_vertex,
            u1, u2, v1, v2,
            su == -1 ? negative_polytope_visitor : positive_polytope_visitor,
            sv == -1 ? negative_polytope_visitor : positive_polytope_visitor
        );
    }
    else {
        assert(topological_case == 0);
        Detail_Divide_Simplex::Divide_Simplex_Cases_0_15(
            alignment_of_sub_simplex_opposite_vertex,
            normal_of_sub_simplex_opposite_vertex,
            case_ == 0 ? negative_polytope_visitor : positive_polytope_visitor
        );
    }
}

namespace Detail_Divide_Simplex
{

template< class T, int N, int D >
inline VECTOR<T,D>
Intersection(
    const BOUNDED_LIST< VECTOR<T,D>, N >& vertices_x,
    const VECTOR< T, D+1 >& phi_of_vertex,
    const int u, const int v)
{
    const T phiu = phi_of_vertex[u];
    const T phiv = phi_of_vertex[v];
    assert(phiu * phiv <= 0);
    return phiu == phiv ?
           (vertices_x(v) + vertices_x(u)) / 2 :
           (phiv * vertices_x(u) - phiu * vertices_x(v)) / (phiv - phiu);
}

template< class T, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_0_3(
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 2 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,1>, 2 >& normal_of_sub_simplex_opposite_vertex,
    T_POLYTOPE_VISITOR polytope_visitor)
{
    polytope_visitor(SUB_CUBE_POLYTOPE<T,1>(
        alignment_of_sub_simplex_opposite_vertex[2],
           normal_of_sub_simplex_opposite_vertex[2],
        1
    ));
    polytope_visitor(SUB_CUBE_POLYTOPE<T,1>(
        alignment_of_sub_simplex_opposite_vertex[1],
           normal_of_sub_simplex_opposite_vertex[1],
        2
    ));
}

template< class T, int N, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_1_2(
    BOUNDED_LIST< VECTOR<T,1>, N >& vertices_x,
    const VECTOR<T,2>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 2 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,1>, 2 >& normal_of_sub_simplex_opposite_vertex,
    const int u, const int v,
    T_POLYTOPE_VISITOR u_polytope_visitor,
    T_POLYTOPE_VISITOR v_polytope_visitor)
{
    assert(vertices_x.Size() == 2);
    /********
     * v---u
     ********/
    assert(u == 1 && v == 2);
    assert(phi_of_vertex[u] * phi_of_vertex[v] <= 0);
    const VECTOR<T,1> x3 = Intersection(vertices_x, phi_of_vertex, u, v);
    vertices_x.Append(x3);
    const SUB_CUBE_POLYTOPE_ALIGNMENT au = alignment_of_sub_simplex_opposite_vertex[u];
    const SUB_CUBE_POLYTOPE_ALIGNMENT av = alignment_of_sub_simplex_opposite_vertex[v];
    static const SUB_CUBE_POLYTOPE_ALIGNMENT a3 = SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED;
    const VECTOR<T,1> nu = normal_of_sub_simplex_opposite_vertex[u];
    const VECTOR<T,1> nv = normal_of_sub_simplex_opposite_vertex[v];
    const VECTOR<T,1> n3 = Normal(x3);
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,1>(av, nv, u));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,1>(a3,-n3, 3));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,1>(au, nu, v));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,1>(a3, n3, 3));
}

template< class T, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_0_7(
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 3 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,2>, 3 >& normal_of_sub_simplex_opposite_vertex,
    T_POLYTOPE_VISITOR polytope_visitor)
{
    polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(
        alignment_of_sub_simplex_opposite_vertex[3],
           normal_of_sub_simplex_opposite_vertex[3],
        1, 2
    ));
    polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(
        alignment_of_sub_simplex_opposite_vertex[1],
           normal_of_sub_simplex_opposite_vertex[1],
        2, 3
    ));
    polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(
        alignment_of_sub_simplex_opposite_vertex[2],
           normal_of_sub_simplex_opposite_vertex[2],
        3, 1
    ));
}

template< class T, int N, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_1_2_3_4_5_6(
    BOUNDED_LIST< VECTOR<T,2>, N >& vertices_x,
    const VECTOR<T,3>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 3 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,2>, 3 >& normal_of_sub_simplex_opposite_vertex,
    const int u, const int v1, const int v2,
    T_POLYTOPE_VISITOR u_polytope_visitor,
    T_POLYTOPE_VISITOR v_polytope_visitor)
{
    assert(vertices_x.Size() == 3);
    /**********
     *    u
     *   / \
     * v1---v2
     **********/
    assert(
        (u == 1 && v1 == 2 && v2 == 3) ||
        (u == 2 && v1 == 3 && v2 == 1) ||
        (u == 3 && v1 == 1 && v2 == 2)
    );
    assert(phi_of_vertex[u ] * phi_of_vertex[v1] <= 0);
    assert(phi_of_vertex[u ] * phi_of_vertex[v2] <= 0);
    assert(phi_of_vertex[v1] * phi_of_vertex[v2] >= 0);
    const VECTOR<T,2> x4 = Intersection(vertices_x, phi_of_vertex, u, v1);
    const VECTOR<T,2> x5 = Intersection(vertices_x, phi_of_vertex, u, v2);
    vertices_x.Append(x4);
    vertices_x.Append(x5);
    const SUB_CUBE_POLYTOPE_ALIGNMENT au  = alignment_of_sub_simplex_opposite_vertex[u ];
    const SUB_CUBE_POLYTOPE_ALIGNMENT av1 = alignment_of_sub_simplex_opposite_vertex[v1];
    const SUB_CUBE_POLYTOPE_ALIGNMENT av2 = alignment_of_sub_simplex_opposite_vertex[v2];
    static const SUB_CUBE_POLYTOPE_ALIGNMENT a45 = SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED;
    const VECTOR<T,2>& nu  = normal_of_sub_simplex_opposite_vertex[u ];
    const VECTOR<T,2>& nv1 = normal_of_sub_simplex_opposite_vertex[v1];
    const VECTOR<T,2>& nv2 = normal_of_sub_simplex_opposite_vertex[v2];
    const VECTOR<T,2> n45 = Normal(x4, x5);
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(av2, nv2, u , 4 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(a45, n45, 4 , 5 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(av1, nv1, 5 , u ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(au , nu , v1, v2));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(av1, nv1, v2, 5 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(a45,-n45, 5 , 4 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,2>(av2, nv2, 4 , v1));
}

template< class T, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_0_15(
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 4 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,3>, 4 >& normal_of_sub_simplex_opposite_vertex,
    T_POLYTOPE_VISITOR polytope_visitor)
{
    polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(
        alignment_of_sub_simplex_opposite_vertex[4],
           normal_of_sub_simplex_opposite_vertex[4],
        1, 2, 3
    ));
    polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(
        alignment_of_sub_simplex_opposite_vertex[2],
           normal_of_sub_simplex_opposite_vertex[2],
        1, 3, 4
    ));
    polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(
        alignment_of_sub_simplex_opposite_vertex[3],
           normal_of_sub_simplex_opposite_vertex[3],
        1, 4, 2
    ));
    polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(
        alignment_of_sub_simplex_opposite_vertex[1],
           normal_of_sub_simplex_opposite_vertex[1],
        2, 4, 3
    ));
}

template< class T, int N, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_1_2_4_7_8_11_13_14(
    BOUNDED_LIST< VECTOR<T,3>, N >& vertices_x,
    const VECTOR<T,4>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 4 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,3>, 4 >& normal_of_sub_simplex_opposite_vertex,
    const int u, const int v1, const int v2, const int v3,
    T_POLYTOPE_VISITOR u_polytope_visitor,
    T_POLYTOPE_VISITOR v_polytope_visitor)
{
    assert(vertices_x.Size() == 4);
    /************
     *     u
     *    /|\
     *   / | \
     * v1--|--v3
     *   \ | /
     *    \|/
     *     v2
     ***********/
    assert(
        (u == 1 && v1 == 2 && v2 == 3 && v3 == 4) ||
        (u == 2 && v1 == 1 && v2 == 4 && v3 == 3) ||
        (u == 3 && v1 == 1 && v2 == 2 && v3 == 4) ||
        (u == 4 && v1 == 1 && v2 == 3 && v3 == 2)
    );
    assert(phi_of_vertex[u ] * phi_of_vertex[v1] <= 0);
    assert(phi_of_vertex[u ] * phi_of_vertex[v2] <= 0);
    assert(phi_of_vertex[u ] * phi_of_vertex[v3] <= 0);
    assert(phi_of_vertex[v1] * phi_of_vertex[v2] >= 0);
    assert(phi_of_vertex[v1] * phi_of_vertex[v3] >= 0);
    assert(phi_of_vertex[v2] * phi_of_vertex[v3] >= 0);
    const VECTOR<T,3> x5 = Intersection(vertices_x, phi_of_vertex, u, v1);
    const VECTOR<T,3> x6 = Intersection(vertices_x, phi_of_vertex, u, v2);
    const VECTOR<T,3> x7 = Intersection(vertices_x, phi_of_vertex, u, v3);
    vertices_x.Append(x5);
    vertices_x.Append(x6);
    vertices_x.Append(x7);
    const SUB_CUBE_POLYTOPE_ALIGNMENT au  = alignment_of_sub_simplex_opposite_vertex[u ];
    const SUB_CUBE_POLYTOPE_ALIGNMENT av1 = alignment_of_sub_simplex_opposite_vertex[v1];
    const SUB_CUBE_POLYTOPE_ALIGNMENT av2 = alignment_of_sub_simplex_opposite_vertex[v2];
    const SUB_CUBE_POLYTOPE_ALIGNMENT av3 = alignment_of_sub_simplex_opposite_vertex[v3];
    static const SUB_CUBE_POLYTOPE_ALIGNMENT a567 = SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED;
    const VECTOR<T,3>& nu  = normal_of_sub_simplex_opposite_vertex[u ];
    const VECTOR<T,3>& nv1 = normal_of_sub_simplex_opposite_vertex[v1];
    const VECTOR<T,3>& nv2 = normal_of_sub_simplex_opposite_vertex[v2];
    const VECTOR<T,3>& nv3 = normal_of_sub_simplex_opposite_vertex[v3];
    const VECTOR<T,3> n567 = Normal(x5, x6, x7);
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av3 , nv3 , u , 5 , 6 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av1 , nv1 , u , 6 , 7 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av2 , nv2 , u , 7 , 5 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a567,-n567, 5 , 7 , 6 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(au  , nu  , v1, v3, v2));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av3 , nv3 , v1, v2, 6 , 5 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av1 , nv1 , v2, v3, 7 , 6 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av2 , nv2 , v3, v1, 5 , 7 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a567, n567, 5 , 6 , 7 ));
}

template< class T, int N, class T_POLYTOPE_VISITOR >
inline void
Divide_Simplex_Cases_3_5_6_9_10_12(
    BOUNDED_LIST< VECTOR<T,3>, N >& vertices_x,
    const VECTOR<T,4>& phi_of_vertex,
    const VECTOR< SUB_CUBE_POLYTOPE_ALIGNMENT, 4 >& alignment_of_sub_simplex_opposite_vertex,
    const VECTOR< VECTOR<T,3>, 4 >& normal_of_sub_simplex_opposite_vertex,
    const int u1, const int u2, const int v1, const int v2,
    T_POLYTOPE_VISITOR u_polytope_visitor,
    T_POLYTOPE_VISITOR v_polytope_visitor)
{
    assert(vertices_x.Size() == 4);
    /*************
     *     u1
     *    /|\
     *   / | \
     * u2--|--v2
     *   \ | /
     *    \|/
     *     v1
     ************/
    assert(u1 == 1);
    assert(
        (u2 == 2 && v1 == 3 && v2 == 4) ||
        (u2 == 3 && v1 == 4 && v2 == 2) ||
        (u2 == 4 && v1 == 2 && v2 == 3)
    );
    assert(phi_of_vertex[u1] * phi_of_vertex[u2] >= 0);
    assert(phi_of_vertex[u1] * phi_of_vertex[v1] <= 0);
    assert(phi_of_vertex[u1] * phi_of_vertex[v2] <= 0);
    assert(phi_of_vertex[u2] * phi_of_vertex[v1] <= 0);
    assert(phi_of_vertex[u2] * phi_of_vertex[v2] <= 0);
    assert(phi_of_vertex[v1] * phi_of_vertex[v2] >= 0);
    const VECTOR<T,3> x5 = Intersection(vertices_x, phi_of_vertex, u1, v1);
    const VECTOR<T,3> x6 = Intersection(vertices_x, phi_of_vertex, u1, v2);
    const VECTOR<T,3> x7 = Intersection(vertices_x, phi_of_vertex, u2, v1);
    const VECTOR<T,3> x8 = Intersection(vertices_x, phi_of_vertex, u2, v2);
    const VECTOR<T,3> x9 = (x5 + x6 + x7 + x8) / 4;
    vertices_x.Append(x5);
    vertices_x.Append(x6);
    vertices_x.Append(x7);
    vertices_x.Append(x8);
    vertices_x.Append(x9);
    const SUB_CUBE_POLYTOPE_ALIGNMENT au1 = alignment_of_sub_simplex_opposite_vertex[u1];
    const SUB_CUBE_POLYTOPE_ALIGNMENT au2 = alignment_of_sub_simplex_opposite_vertex[u2];
    const SUB_CUBE_POLYTOPE_ALIGNMENT av1 = alignment_of_sub_simplex_opposite_vertex[v1];
    const SUB_CUBE_POLYTOPE_ALIGNMENT av2 = alignment_of_sub_simplex_opposite_vertex[v2];
    static const SUB_CUBE_POLYTOPE_ALIGNMENT a569 = SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED;
    static const SUB_CUBE_POLYTOPE_ALIGNMENT a689 = SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED;
    static const SUB_CUBE_POLYTOPE_ALIGNMENT a879 = SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED;
    static const SUB_CUBE_POLYTOPE_ALIGNMENT a759 = SUB_CUBE_POLYTOPE_ALIGNMENT_UNALIGNED;
    const VECTOR<T,3>& nu1 = normal_of_sub_simplex_opposite_vertex[u1];
    const VECTOR<T,3>& nu2 = normal_of_sub_simplex_opposite_vertex[u2];
    const VECTOR<T,3>& nv1 = normal_of_sub_simplex_opposite_vertex[v1];
    const VECTOR<T,3>& nv2 = normal_of_sub_simplex_opposite_vertex[v2];
    const VECTOR<T,3> n569 = Normal(x5, x6, x9);
    const VECTOR<T,3> n689 = Normal(x6, x8, x9);
    const VECTOR<T,3> n879 = Normal(x8, x7, x9);
    const VECTOR<T,3> n759 = Normal(x7, x5, x9);
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av2 , nv2 , u1, u2, 7 , 5 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av1 , nv1 , u2, u1, 6 , 8 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(au2 , nu2 , u1, 5 , 6 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(au1 , nu1 , u2, 8 , 7 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a759,-n759, 5 , 7 , 9 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a879,-n879, 7 , 8 , 9 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a689,-n689, 8 , 6 , 9 ));
    u_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a569,-n569, 6 , 5 , 9 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(au2 , nu2 , v1, v2, 6 , 5 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(au1 , nu1 , v2, v1, 7 , 8 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av2 , nv2 , v1, 5 , 7 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(av1 , nv1 , v2, 8 , 6 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a569, n569, 5 , 6 , 9 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a689, n689, 6 , 8 , 9 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a879, n879, 8 , 7 , 9 ));
    v_polytope_visitor(SUB_CUBE_POLYTOPE<T,3>(a759, n759, 7 , 5 , 9 ));
}

} // namespace Detail_Divide_Simplex

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_DIVIDE_SIMPLEX_HPP
