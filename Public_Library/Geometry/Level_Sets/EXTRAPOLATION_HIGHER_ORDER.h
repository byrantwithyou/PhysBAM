//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_HIGHER_ORDER
//#####################################################################
#ifndef __EXTRAPOLATION_HIGHER_ORDER__
#define __EXTRAPOLATION_HIGHER_ORDER__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <functional>
namespace PhysBAM{

// TODO: Limit to near interface
template<class TV,class T2>
class EXTRAPOLATION_HIGHER_ORDER
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
public:
    struct STENCIL {VECTOR<int,4> nodes;T scale;};

    GRID<TV> grid;
    const LEVELSET<TV>& phi;
    int iterations;
    int order;
    int fill_width;
    ARRAY<int,TV_INT> node_to_index;
    ARRAY<TV_INT> index_to_node;
    VECTOR<INTERVAL<int>,3> fill_indices,solve_indices;
    ARRAY<TV> normal;
    bool periodic;
    VECTOR<bool,TV::m> combine_ends;

    EXTRAPOLATION_HIGHER_ORDER(const GRID<TV>& grid,const LEVELSET<TV>& phi,int iterations,int order,int fill_width);
    EXTRAPOLATION_HIGHER_ORDER(const EXTRAPOLATION_HIGHER_ORDER&) = delete;
    void operator=(const EXTRAPOLATION_HIGHER_ORDER&) = delete;
    ~EXTRAPOLATION_HIGHER_ORDER();

    void Extrapolate_Node(std::function<bool(const TV_INT& index)> inside_mask,ARRAYS_ND_BASE<T2,TV_INT>& x);
    void Extrapolate_Cell(std::function<bool(const TV_INT& index)> inside_mask,ARRAYS_ND_BASE<T2,TV_INT>& x);
    void Extrapolate_Face(std::function<bool(const FACE_INDEX<TV::m>& index)> inside_mask,ARRAY<T2,FACE_INDEX<TV::m> >& x);
    int Lookup_Index(TV_INT index) const;
    int Register_Index(TV_INT index,int only_neg=0);
    void Periodic_Index(TV_INT& index) const;
protected:
    void Add_Neighbors(ARRAY<TV_INT>& next,const ARRAY<TV_INT>& neighbors,const TV_INT& index,int unregistered,int registered);
    void Register_Nodes(std::function<bool(const TV_INT& index)> inside_mask,ARRAY<VECTOR<STENCIL,TV::m> >& stencil);
    void Extrapolate_FE(const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,const ARRAY<T2>& x,ARRAY<T2>& y,const ARRAY<T2>* z,int o,T dt,T alpha);
    void Extrapolate_RK2(const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,ARRAY<T2>& x,const ARRAY<T2>* z,ARRAY<T2>& tmp,int o,T dt);
    void Fill_un(const TV& one_over_dx,const ARRAY<T2>& x,ARRAY<T2>& xn,int o,int mo);
//#####################################################################
};
}
#endif
