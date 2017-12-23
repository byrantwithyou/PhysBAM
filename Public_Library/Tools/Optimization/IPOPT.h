//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IPOPT
//##################################################################### 
#ifndef __IPOPT__
#define __IPOPT__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T>
class IPOPT
{
public:
    T tolerance;
    
    ARRAY<INTERVAL<T> > dof_range;
    ARRAY<INTERVAL<T> > constraint_range;
    ARRAY<T> initial_guess;
    ARRAY<T>& solution=initial_guess;

    // Must be initialized with the correct maximal sparsity pattern.  This may not be changed later.
    // Should be filled in when Compute_Hessian and Compute_Constraint_Jacobian are called.
    HASHTABLE<VECTOR<int,2>,T> H; // objective Hessian (only lower triangular will be accessed)
    HASHTABLE<VECTOR<int,2>,T> J; // constraint Jacobian

    IPOPT()=default;
    IPOPT(const IPOPT&)=delete;
    virtual ~IPOPT()=default;
    IPOPT& operator=(const IPOPT&)=delete;

    bool Solve();

    virtual T Compute_Objective(const ARRAY<T>& x)=0;
    virtual void Compute_Gradient(ARRAY<T>& g,const ARRAY<T>& x)=0;
    virtual void Compute_Hessian(const ARRAY<T>& x)=0;
    virtual void Compute_Constraints(ARRAY<T>& f,const ARRAY<T>& x)=0;
    virtual void Compute_Constraint_Jacobian(const ARRAY<T>& x)=0;
};
}
#endif
