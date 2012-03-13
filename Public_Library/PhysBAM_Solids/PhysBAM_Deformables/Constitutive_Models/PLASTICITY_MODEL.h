//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLASTICITY_MODEL
//##################################################################### 
#ifndef __PLASTICITY_MODEL__
#define __PLASTICITY_MODEL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T,int d>
class PLASTICITY_MODEL:public NONCOPYABLE
{
    typedef VECTOR<T,d> TV;
public:
    ARRAY<SYMMETRIC_MATRIX<T,d> > Fp_inverse;

    PLASTICITY_MODEL(const int elements)
    {
        Fp_inverse.Exact_Resize(elements);
        Fp_inverse.Fill(SYMMETRIC_MATRIX<T,d>::Identity_Matrix());
    }
    
    virtual ~PLASTICITY_MODEL()
    {}

    virtual void Read_State(TYPED_ISTREAM& input)
    {Read_Binary(input,Fp_inverse);}

    virtual void Write_State(TYPED_OSTREAM& output) const
    {Write_Binary(output,Fp_inverse);}

//#####################################################################
    virtual bool Project_Fe(const DIAGONAL_MATRIX<T,d>& Fe_trial,DIAGONAL_MATRIX<T,d>& Fe_project) const=0;
    virtual void Project_Fp(const int tetrahedron,const MATRIX<T,d>& Fp_trial)=0;
//#####################################################################
};
}
#endif
