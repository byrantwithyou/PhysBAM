//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM__
#define __LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM__

#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<TV>
class LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM:public INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP> interpolation;

    LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM()
    {}

    TV From_Block_Face(const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const override
    {u.Set_Reference_Point(X);TV result=interpolation.From_Block_Face(grid,block,u,X);u.Clear_Reference_Point();return result;}

    T From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const override
    {u.Set_Reference_Point(X);T result=interpolation.From_Block_Face_Component(axis,grid,block,u,X);u.Clear_Reference_Point();return result;}

    VECTOR<TV,2> Extrema_From_Block_Face(const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const override
    {u_min.Set_Reference_Point(X);u_max.Set_Reference_Point(X);VECTOR<TV,2> result=interpolation.Extrema_From_Block_Face(grid,block,u_min,u_max,X);
    u_min.Clear_Reference_Point();u_max.Clear_Reference_Point();return result;}

    VECTOR<T,2> Extrema_From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,
        const typename T_FACE_LOOKUP::LOOKUP& u_min,const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const override
    {u_min.Set_Reference_Point(X);u_max.Set_Reference_Point(X);VECTOR<T,2> result=interpolation.Extrema_From_Block_Face_Component(axis,grid,block,u_min,u_max,X);
    u_min.Clear_Reference_Point();u_max.Clear_Reference_Point();return result;}

//#####################################################################
};
}
#endif
