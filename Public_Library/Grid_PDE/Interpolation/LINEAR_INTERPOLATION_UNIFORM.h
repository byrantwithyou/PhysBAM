//#####################################################################
// Copyright 2005-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_UNIFORM__
#define __LINEAR_INTERPOLATION_UNIFORM__

#include <Core/Data_Structures/PAIR.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/Componentwise_Min_Max.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC_1D_HELPER.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC_HELPER.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV>
class LINEAR_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    template<class T3> struct REBIND{typedef LINEAR_INTERPOLATION_UNIFORM<TV,T3,T_FACE_LOOKUP> TYPE;};

    TV From_Block_Face(const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const override
    {return LINEAR_INTERPOLATION_MAC_HELPER<TV>::Interpolate_Face(block,u,X);}

    VECTOR<TV,2> Extrema_From_Block_Face(const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const override
    {return LINEAR_INTERPOLATION_MAC_HELPER<TV>::Extrema_Face(block,u_min,u_max,X);}

    VECTOR<T,2> Extrema_From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,
        const typename T_FACE_LOOKUP::LOOKUP& u_min,const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const override
    {return LINEAR_INTERPOLATION_MAC_HELPER<TV>::Extrema_Face_Component(axis,block,u_min,u_max,X);}
    
    T2 Clamped_To_Array_No_Extrema(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    VECTOR<T2,TV::m> Clamped_To_Array_Gradient(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    SYMMETRIC_MATRIX<T2,TV::m> Clamped_To_Array_Hessian(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const override;
    ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const override;
    VECTOR<T2,2> Extrema_Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X) const  override;
    VECTOR<T2,2> Extrema_From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X,const TV_INT& index) const override;
    T From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const override;
    ARRAY<PAIR<FACE_INDEX<TV::m>,T> > From_Block_Face_Component_Weights(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const override;
//#####################################################################
};
}
#endif
