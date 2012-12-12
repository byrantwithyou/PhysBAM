//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_HIGHER_ORDER
//#####################################################################
#ifndef __EXTRAPOLATION_HIGHER_ORDER__
#define __EXTRAPOLATION_HIGHER_ORDER__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <boost/function.hpp>
namespace PhysBAM{

// TODO: Limit to near interface
template<class TV,class T2>
class EXTRAPOLATION_HIGHER_ORDER:public NONCOPYABLE
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    struct MAPPING
    {
        ARRAY<int,TV_INT> node_to_index;
        ARRAY<TV_INT> index_to_node;
        VECTOR<INTERVAL<int>,3> fill_indices,solve_indices;
    };
    struct STENCIL {VECTOR<int,4> nodes;T scale;};

public:
    EXTRAPOLATION_HIGHER_ORDER();
    ~EXTRAPOLATION_HIGHER_ORDER();

    static void Extrapolate_Node(const GRID<TV>& grid,const LEVELSET<TV>& phi,boost::function<bool(const TV_INT& index)> inside_mask,int ghost,ARRAYS_ND_BASE<T2,TV_INT>& x,int iterations,int order,int fill_width);
    static void Extrapolate_Cell(const GRID<TV>& grid,const LEVELSET<TV>& phi,boost::function<bool(const TV_INT& index)> inside_mask,int ghost,ARRAYS_ND_BASE<T2,TV_INT>& x,int iterations,int order,int fill_width);
    static void Extrapolate_Face(const GRID<TV>& grid,const LEVELSET<TV>& phi,boost::function<bool(const FACE_INDEX<TV::m>& index)> inside_mask,int ghost,ARRAY<T2,FACE_INDEX<TV::m> >& x,int iterations,int order,int fill_width);
protected:
    static void Add_Neighbors(MAPPING& m,ARRAY<TV_INT>& next,const ARRAY<TV_INT>& neighbors,const TV_INT& index,int unregistered,int registered);
    static void Register_Nodes(const GRID<TV>& grid,const LEVELSET<TV>& phi,boost::function<bool(const TV_INT& index)> inside_mask,int ghost,MAPPING& m,ARRAY<TV>& normal,
        ARRAY<VECTOR<STENCIL,TV::m> >& stencil,int order,int fill_width);
    static void Extrapolate_FE(const MAPPING& m,const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,const ARRAY<T2>& x,ARRAY<T2>& y,const ARRAY<T2>* z,int o,T dt,T alpha);
    static void Extrapolate_RK2(const MAPPING& m,const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,ARRAY<T2>& x,const ARRAY<T2>* z,ARRAY<T2>& tmp,int o,T dt);
    static void Fill_un(const MAPPING& m,const TV& one_over_dx,const ARRAY<TV>& normal,const ARRAY<T2>& x,ARRAY<T2>& xn,int o,int mo);
//#####################################################################
};
}
#endif
