//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADRATIC_PROGRAMMING
//#####################################################################
#ifndef __QUADRATIC_PROGRAMMING__
#define __QUADRATIC_PROGRAMMING__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Matrices/MATRIX_MXN.h>

namespace PhysBAM{

template<class T>
class QUADRATIC_PROGRAMMING
{
public:
    QUADRATIC_PROGRAMMING()
    {}

    ~QUADRATIC_PROGRAMMING()
    {}

//##################################################################### 
private:
    static void Move_Column_From_B_To_N_And_Shift_Down(MATRIX_MXN<T>& B,ARRAY<int>& permute_B,const int b_column,MATRIX_MXN<T>& N,ARRAY<int>& permute_N,const int n_column);
    static void Remove_Column(const int column,MATRIX_MXN<T>& A,ARRAY<int>& permute,ARRAY<T>* x=0);
    static void Update_Upper_Triangular_Matrix_After_Column_Shift(MATRIX_MXN<T>& A,MATRIX_MXN<T>& S,MATRIX_MXN<T>& N,ARRAY<T>& b,const int column,const T tolerance,const bool check_last_column=true);
public:
    static void Find_Optimal_Solution(MATRIX_MXN<T>& B,MATRIX_MXN<T>& S,MATRIX_MXN<T>& N,ARRAY<T>& x_B,ARRAY<T>& x_S,ARRAY<T>& b,ARRAY<T>& b_N,ARRAY<int>& permute_B,
                                      ARRAY<int>& permute_S,ARRAY<int>& permute_N,const MATRIX_MXN<T>& D_unpermuted,const MATRIX_MXN<T>& epsilon_hat_unpermuted,const ARRAY<T>& f_hat,
                                      ARRAY<T>& x_unpermuted,const ARRAY<PAIR<bool,T> >& x_min,const ARRAY<PAIR<bool,T> >& x_max,const T tolerance,const T step_tolerance,
                                      const bool debug_optimization=false);
//##################################################################### 
};
}
#endif
