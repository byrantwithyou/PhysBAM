//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EIGENSYSTEM  
//##################################################################### 
#ifndef __EIGENSYSTEM__
#define __EIGENSYSTEM__   

#include <Core/Vectors/VECTOR_3D.h>
#include <Grid_Tools/Arrays/ARRAYS_UNIFORM_FORWARD.h>
namespace PhysBAM{

template<class T,int m,int n> class MATRIX;

template<class T,class TV_DIMENSION>
class EIGENSYSTEM
{
public:
    VECTOR<int,3> slice_index;
protected:
    enum WORKAROUND1 {d=TV_DIMENSION::m};
    EIGENSYSTEM()
    {}
public:
    virtual ~EIGENSYSTEM()
    {}

    virtual bool All_Eigenvalues_Same()
    {return false;}

//#####################################################################
    virtual void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,
        ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0)=0;
    virtual void Flux_Divided_By_Velocity(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,
        ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Get_Face_Velocity_Component(const int face_index,const bool use_standard_average,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,
        ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0)
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right)=0;
    virtual void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)=0;
//#####################################################################
};   
}
#endif

