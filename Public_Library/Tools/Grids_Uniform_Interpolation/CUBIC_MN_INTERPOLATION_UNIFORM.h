//#####################################################################
// Copyright 2003-2010, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUBIC_MN_INTERPOLATION_UNIFORM 
//#####################################################################
#ifndef __CUBIC_MN_INTERPOLATION_UNIFORM__
#define __CUBIC_MN_INTERPOLATION_UNIFORM__

#include <Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <Tools/Interpolation/CUBIC_MN_INTERPOLATION.h>
namespace PhysBAM{

template<class TV,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV>
class CUBIC_MN_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP> BASE;
    using BASE::Clamped_Index_Interior_End_Minus_One;

    T b,c; // two parameter family
    CUBIC_MN_INTERPOLATION<T,T2> cubic_mn_interpolation;

    CUBIC_MN_INTERPOLATION_UNIFORM();
    ~CUBIC_MN_INTERPOLATION_UNIFORM();

    void Set_Sharpness(T sharpness_input=2./3)
    {Set_Parameters(1-sharpness_input,(T).5*sharpness_input);}

    void Set_Parameters(const T b_input=1./3,const T c_input=1./3)
    {b=b_input;c=c_input;cubic_mn_interpolation.Set_Parameters(b,c);}

    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    T2 Periodic(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const;
    // T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;
    // T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
    // T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const;
    // ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;
    // ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
    // ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const;
    // only works for power of two grids!
//    T2 From_Base_Node_Periodic(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
//#####################################################################
};
}
#endif
