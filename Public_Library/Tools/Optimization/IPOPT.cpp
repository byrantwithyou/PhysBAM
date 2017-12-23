//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Tools/Optimization/IPOPT.h>
#ifdef USE_IPOPT
#define HAVE_CSTDDEF
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>
#endif

namespace PhysBAM{
#ifdef USE_IPOPT
template<class T>
class IPOPT_CORE : public Ipopt::TNLP
{
public:
    IPOPT<T>& ipopt;

    ARRAY<T> tmp_dofs,tmp_dofs2,tmp_const;
    ARRAY<VECTOR<int,2> > nz_H,nz_J;

    IPOPT_CORE(IPOPT<T>& ipopt):ipopt(ipopt) {}
    IPOPT_CORE(const IPOPT_CORE&)=delete;
    virtual ~IPOPT_CORE()=default;
    IPOPT_CORE& operator=(const IPOPT_CORE&)=delete;

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
};

//#####################################################################
// Function get_nlp_info
//#####################################################################
template<class T> bool IPOPT_CORE<T>::
get_nlp_info(Index& n,Index& m,Index& nnz_jac_g,
    Index& nnz_h_lag,IndexStyleEnum& index_style)
{
    n=ipopt.dof_range.m;
    m=ipopt.constraint_range.m;
    tmp_dofs.Resize(ipopt.dof_range.m);
    tmp_dofs2.Resize(ipopt.dof_range.m);
    tmp_const.Resize(ipopt.constraint_range.m);

    for(auto it:ipopt.H) if(it.key.x>=it.key.y) nz_H.Append(it.key);
    nz_H.Sort(LEXICOGRAPHIC_COMPARE());
    nnz_h_lag=nz_H.m;

    for(auto it:ipopt.J) nz_J.Append(it.key);
    nz_J.Sort(LEXICOGRAPHIC_COMPARE());
    nnz_jac_g=nz_J.m;

    index_style=TNLP::C_STYLE;
    return true;
}
//#####################################################################
// Function get_bounds_info
//#####################################################################
template<class T> bool IPOPT_CORE<T>::
get_bounds_info(Index n,Number* x_l,Number* x_u,
    Index m,Number* g_l,Number* g_u)
{
    for(int i=0; i<ipopt.dof_range.m; i++){
        x_l[i]=ipopt.dof_range(i).min_corner;
        x_u[i]=ipopt.dof_range(i).max_corner;}

    for(int i=0; i<ipopt.constraint_range.m; i++){
        g_l[i]=ipopt.constraint_range(i).min_corner;
        g_u[i]=ipopt.constraint_range(i).max_corner;}

    return true;
}
//#####################################################################
// Function get_starting_point
//#####################################################################
template<class T> bool IPOPT_CORE<T>::
get_starting_point(Index n,bool init_x,Number* x,
    bool init_z,Number* z_L,Number* z_U,
    Index m,bool init_lambda,
    Number* lambda)
{
    // Here,we assume we only have starting values for x,ifyou code
    // your own NLP,you can provide starting values for the dual variables
    // ifyou wish
    assert(init_x==true);
    assert(init_z==false);
    assert(init_lambda==false);

    for(int i=0;i<ipopt.initial_guess.m;i++)
        x[i]=ipopt.initial_guess(i);

    return true;
}
//#####################################################################
// Function eval_f
//#####################################################################
template<class T> bool IPOPT_CORE<T>::
eval_f(Index n,const Number* x,bool new_x,Number& obj_value)
{
    for(int i=0;i<tmp_dofs.m;i++) tmp_dofs(i)=x[i];
    obj_value=ipopt.Compute_Objective(tmp_dofs);
    return true;
}
//#####################################################################
// Function eval_grad_f
//#####################################################################
template<class T> bool IPOPT_CORE<T>::
eval_grad_f(Index n,const Number* x,bool new_x,Number* grad_f)
{
    for(int i=0;i<tmp_dofs.m;i++) tmp_dofs(i)=x[i];
    ipopt.Compute_Gradient(tmp_dofs2,tmp_dofs);
    for(int i=0;i<tmp_dofs2.m;i++) grad_f[i]=tmp_dofs2(i);
    return true;
}
//#####################################################################
// Function eval_g
//#####################################################################
template<class T> bool IPOPT_CORE<T>::
eval_g(Index n,const Number* x,bool new_x,Index m,Number* g)
{
    for(int i=0;i<tmp_dofs.m;i++) tmp_dofs(i)=x[i];
    ipopt.Compute_Constraints(tmp_const,tmp_dofs);
    for(int i=0;i<tmp_const.m;i++) g[i]=tmp_const(i);
    return true;
}
//#####################################################################
// Function eval_jac_g
//#####################################################################
template<class T> bool IPOPT_CORE<T>::
eval_jac_g(Index n,const Number* x,bool new_x,
    Index m,Index nele_jac,Index* iRow,Index *jCol,
    Number* values)
{
    if(values==NULL){
        for(int i=0;i<nz_J.m;i++){
            iRow[i]=nz_J(i).x;
            jCol[i]=nz_J(i).y;}}
    else{
        for(int i=0;i<tmp_dofs.m;i++) tmp_dofs(i)=x[i];
        ipopt.J.Remove_All();
        ipopt.Compute_Constraint_Jacobian(tmp_dofs);
        for(int i=0;i<nz_J.m;i++)
            values[i]=ipopt.J.Get_Or_Insert(nz_J(i));}
    return true;
}
//#####################################################################
// Function eval_h
//#####################################################################
template<class T> bool IPOPT_CORE<T>::
eval_h(Index n,const Number* x,bool new_x,
    Number obj_factor,Index m,const Number* lambda,
    bool new_lambda,Index nele_hess,Index* iRow,
    Index* jCol,Number* values)
{
    if(values==NULL){
        for(int i=0;i<nz_H.m;i++){
            iRow[i]=nz_H(i).x;
            jCol[i]=nz_H(i).y;}}
    else{
        for(int i=0;i<tmp_dofs.m;i++) tmp_dofs(i)=x[i];
        ipopt.H.Remove_All();
        ipopt.Compute_Hessian(tmp_dofs);
        for(int i=0;i<nz_H.m;i++)
            values[i]=ipopt.H.Get_Or_Insert(nz_H(i));}
    return true;
}
//#####################################################################
// Function finalize_solution
//#####################################################################
template<class T> void IPOPT_CORE<T>::
finalize_solution(SolverReturn status,
    Index n,const Number* x,const Number* z_L,const Number* z_U,
    Index m,const Number* g,const Number* lambda,
    Number obj_value,
    const IpoptData* ip_data,
    IpoptCalculatedQuantities* ip_cq)
{
    for(int i=0;i<ipopt.solution.m;i++) ipopt.solution(i)=x[i];
}
#endif
//#####################################################################
// Function Solve
//#####################################################################
template<class T> bool IPOPT<T>::
Solve()
{
#ifndef USE_IPOPT
    Ipopt::IpoptApplication app;
    app.RethrowNonIpoptException(true);
    app.Options()->SetNumericValue("tol",tolerance);
    app.Options()->SetStringValue("mu_strategy","adaptive");
    app.Options()->SetIntegerValue("print_level",0);

    Ipopt::ApplicationReturnStatus status;
    status=app.Initialize();
    PHYSBAM_ASSERT(status==Ipopt::Solve_Succeeded);

    // Ask Ipopt to solve the problem
    Ipopt::SmartPtr<Ipopt::TNLP> ptr(new IPOPT_CORE<T>(*this));
    status=app.OptimizeTNLP(ptr);
    return status==Ipopt::Solve_Succeeded;
#else
    PHYSBAM_FATAL_ERROR("Must compile with IPOPT support to use this routine.");
#endif
}
template class IPOPT<float>;
template class IPOPT<double>;
}
