//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IPOPT
//##################################################################### 
#ifdef USE_IPOPT
#ifndef __IPOPT__
#define __IPOPT__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Vectors/VECTOR.h>
#define HAVE_CSTDDEF
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>
namespace PhysBAM{

template<class T>
class IPOPT : public Ipopt::TNLP
{
public:
    T tolerance;
    
    ARRAY<INTERVAL<T> > dof_range;
    ARRAY<INTERVAL<T> > constraint_range;
    ARRAY<T> initial_guess;
    ARRAY<T>& solution=tmp_dofs2;

    // Must be initialized with the correct maximal sparsity pattern.  This may not be changed later.
    // Should be filled in when Compute_Hessian and Compute_Constraint_Jacobian are called.
    HASHTABLE<VECTOR<int,2>,T> H; // objective Hessian (only lower triangular will be accessed)
    HASHTABLE<VECTOR<int,2>,T> J; // constraint Jacobian

protected:
    ARRAY<T> tmp_dofs,tmp_dofs2,tmp_const;
    ARRAY<VECTOR<int,2> > nz_H,nz_J;
public:

    IPOPT()=default;
    virtual ~IPOPT()=default;

    bool Solve();

    virtual T Compute_Objective(const ARRAY<T>& x)=0;
    virtual void Compute_Gradient(ARRAY<T>& g,const ARRAY<T>& x)=0;
    virtual void Compute_Hessian(const ARRAY<T>& x)=0;
    virtual void Compute_Constraints(ARRAY<T>& f,const ARRAY<T>& x)=0;
    virtual void Compute_Constraint_Jacobian(const ARRAY<T>& x)=0;

protected:
    typedef Ipopt::Index Index;typedef Ipopt::Number Number;
    typedef Ipopt::SolverReturn SolverReturn;typedef Ipopt::IpoptData IpoptData;
    typedef Ipopt::IpoptCalculatedQuantities IpoptCalculatedQuantities;
    virtual bool get_nlp_info(Index& n,Index& m,Index& nnz_jac_g,
        Index& nnz_h_lag,IndexStyleEnum& index_style) override;
    virtual bool get_bounds_info(Index n,Number* x_l,Number* x_u,Index m,
        Number* g_l,Number* g_u) override;
    virtual bool get_starting_point(Index n,bool init_x,Number* x,bool init_z,
        Number* z_L,Number* z_U,Index m,bool init_lambda,Number* lambda) override;
    virtual bool eval_f(Index n,const Number* x,bool new_x,Number& obj_value) override;
    virtual bool eval_grad_f(Index n,const Number* x,bool new_x,Number* grad_f) override;
    virtual bool eval_g(Index n,const Number* x,bool new_x,Index m,Number* g) override;
    virtual bool eval_jac_g(Index n,const Number* x,bool new_x,Index m,
        Index nele_jac,Index* iRow,Index *jCol,Number* values) override;
    virtual bool eval_h(Index n,const Number* x,bool new_x,Number obj_factor,
        Index m,const Number* lambda,bool new_lambda,Index nele_hess,
        Index* iRow,Index* jCol,Number* values) override;
    virtual void finalize_solution(SolverReturn status,Index n,const Number* x,
        const Number* z_L,const Number* z_U,Index m,const Number* g,
        const Number* lambda,Number obj_value,const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq) override;

private:
    IPOPT(const IPOPT&);
    IPOPT& operator=(const IPOPT&);
};
}
#endif
#endif
