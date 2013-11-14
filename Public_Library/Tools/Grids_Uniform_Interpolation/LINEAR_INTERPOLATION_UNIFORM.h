//#####################################################################
// Copyright 2005-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_UNIFORM__
#define __LINEAR_INTERPOLATION_UNIFORM__

#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_1D_HELPER.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_HELPER.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Math_Tools/Componentwise_Min_Max.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV>
class LINEAR_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    template<class T3> struct REBIND{typedef LINEAR_INTERPOLATION_UNIFORM<TV,T3,T_FACE_LOOKUP> TYPE;};

    TV From_Block_Face(const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {return LINEAR_INTERPOLATION_MAC_HELPER<TV>::Interpolate_Face(block,u,X);}

    VECTOR<TV,2> Extrema_From_Block_Face(const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {return LINEAR_INTERPOLATION_MAC_HELPER<TV>::Extrema_Face(block,u_min,u_max,X);}

    VECTOR<T,2> Extrema_From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,
        const typename T_FACE_LOOKUP::LOOKUP& u_min,const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {return LINEAR_INTERPOLATION_MAC_HELPER<TV>::Extrema_Face_Component(axis,block,u_min,u_max,X);}
    
    T2 Clamped_To_Array_No_Extrema(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const PHYSBAM_OVERRIDE;
    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const PHYSBAM_OVERRIDE;
    VECTOR<T2,TV::m> Clamped_To_Array_Gradient(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const PHYSBAM_OVERRIDE;
    SYMMETRIC_MATRIX<T2,TV::m> Clamped_To_Array_Hessian(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const PHYSBAM_OVERRIDE;
    T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE;
    ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const PHYSBAM_OVERRIDE;
    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE;
    VECTOR<T2,2> Extrema_Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X) const  PHYSBAM_OVERRIDE;
    VECTOR<T2,2> Extrema_From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE;
    T From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE;
    ARRAY<PAIR<FACE_INDEX<TV::m>,T> > From_Block_Face_Component_Weights(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
