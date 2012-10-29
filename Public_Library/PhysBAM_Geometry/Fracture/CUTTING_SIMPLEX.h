//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CUTTING_SIMPLEX__
#define __CUTTING_SIMPLEX__
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_VECTOR_OPERATION_POLICY.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_FLOAT.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/ulp.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Barycentric_Coordinates.h>
#include <PhysBAM_Geometry/Fracture/CUTTING_SIMPLEX_COORDINATE.h>
#include <cassert>
#include <iostream>

namespace PhysBAM{

template<class T,int d> struct CUTTING_SIMPLEX;
template<class T,int d> struct HASH_REDUCE<CUTTING_SIMPLEX<T,d> >
{static int H(const CUTTING_SIMPLEX<T,d>& s);};

template<class T,int d>
struct CUTTING_SIMPLEX
{
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef CUTTING_SIMPLEX_COORDINATE<T,d> GET_ADAPTIVE_WEIGHTS_RESULT_TYPE;
    enum SIMPLEX_TYPE {GLOBAL_EMBEDDING_FACE,LOCAL_EMBEDDING_FACE,GLOBAL_CUT_FACE,LOCAL_CUT_FACE};

    SIMPLEX_TYPE type;
    int parent;
    int element_owner;
    VECTOR<int,d> nodes;
    VECTOR<VECTOR<T,d>,d> weights;
    T abs_tol;
    VECTOR<VECTOR<T,d>,d+1> element_original_coordinates;
    VECTOR<VECTOR<T,d>,d> simplex_original_coordinates;
    VECTOR<bool,d> node_in_embedded_simplex;

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,type,parent,element_owner,nodes,weights,abs_tol,element_original_coordinates,simplex_original_coordinates,node_in_embedded_simplex);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,type,parent,element_owner,nodes,weights,abs_tol,element_original_coordinates,simplex_original_coordinates,node_in_embedded_simplex);}

    friend struct HASH_REDUCE<CUTTING_SIMPLEX<T,d> >;

    CUTTING_SIMPLEX();
    CUTTING_SIMPLEX(const VECTOR<int,d>& nodes,SIMPLEX_TYPE type);
    CUTTING_SIMPLEX(const VECTOR<int,d>& nodes,const SIMPLEX_TYPE type,const VECTOR<VECTOR<T,d>,d>& weights,const int parent,const int element_owner);
    CUTTING_SIMPLEX(const VECTOR<int,d>& nodes,SIMPLEX_TYPE type,T abs_tol,int parent,int element_owner,const VECTOR<VECTOR<T,d>,d+1>& element_original_coordinates,
        const VECTOR<VECTOR<T,d>,d>& simplex_original_coordinates);
    void Get_Adaptive_Weights(VECTOR<VECTOR<GET_ADAPTIVE_WEIGHTS_RESULT_TYPE,d>,d>& adaptive_weights) const;
    template<int N>
    void Get_Adaptive_Weights(VECTOR<VECTOR<GET_ADAPTIVE_WEIGHTS_RESULT_TYPE,d>,N>& adaptive_weights,const VECTOR<int,N>& node_indices) const;
    void Get_Adaptive_Weights(VECTOR<GET_ADAPTIVE_WEIGHTS_RESULT_TYPE,d>& adaptive_weights,int node_index) const;
private:
    void Get_Adaptive_Weights(VECTOR<GET_ADAPTIVE_WEIGHTS_RESULT_TYPE,d>& adaptive_weights,int node_index,CUTTING_SIMPLEX_COORDINATE_SHARED<T,d>* shared) const;
//#####################################################################
};
}
#endif
