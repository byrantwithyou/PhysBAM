//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_PROGRAMMING
//#####################################################################
#ifndef __LINEAR_PROGRAMMING__
#define __LINEAR_PROGRAMMING__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Matrices/MATRIX_MXN.h>

namespace PhysBAM{

template<class T>
class LINEAR_PROGRAMMING
{
public:
    LINEAR_PROGRAMMING()
    {}

    ~LINEAR_PROGRAMMING()
    {}

//##################################################################### 
private:
    static void Move_Column_From_B_To_N_And_Shift_Down(MATRIX_MXN<T>& B,ARRAY<int>& permute_B,const int b_column,MATRIX_MXN<T>& N,ARRAY<int>& permute_N,const int n_column);
    static void Update_Upper_Triangular_Matrix_After_Column_Shift(MATRIX_MXN<T>& A,MATRIX_MXN<T>& N,ARRAY<T>& b,const int column,const T tolerance,const bool check_last_column=true);
public:
    static void Find_Feasible_Solution(MATRIX_MXN<T>& B,MATRIX_MXN<T>& N,ARRAY<T>& x_B,ARRAY<T>& b,ARRAY<T>& b_N,ARRAY<int>& permute_B,ARRAY<int>& permute_N,
        ARRAY<T>& x_unpermuted,const ARRAY<PAIR<bool,T> >& x_min,const ARRAY<PAIR<bool,T> >& x_max,const T tolerance,const T step_tolerance,bool debug_optimization=false);
//##################################################################### 
};
}
#endif
