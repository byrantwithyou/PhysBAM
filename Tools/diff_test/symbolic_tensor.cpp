//#####################################################################
// Copyright 2015, Daniel Ram, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <cctype>
#include "symbolic_tensor.h"
symbolic_tensor::symbolic_tensor()
{
}
symbolic_tensor::symbolic_tensor(const std::string& indices,const vector<int>& size)
{
    assert((int)indices.size()==size.size());
    set(indices,size);
}
symbolic_tensor::symbolic_tensor(const std::string& indices,int fixed_size)
{
    set(indices,fixed_size);
}
symbolic_tensor::symbolic_tensor(const std::string& indices,const symbolic_tensor& a)
{
    set(indices,a);
}
void symbolic_tensor::set(const std::string& ind,const vector<int>& new_size)
{
    indices=ind;
    size=new_size;
    compute_strides();
}
void symbolic_tensor::set(const std::string& ind,int fixed_size)
{
    indices=ind;
    size.resize(indices.size());
    fill(size.begin(),size.end(),fixed_size);
    compute_strides();
}
static void reordered_strides(vector<int>& strides,const symbolic_tensor& a,const std::string& indices)
{
    strides.resize(indices.size());
    for(int i=0;i<(int)indices.size();i++){
        int f=a.indices.find(indices[i]);
        strides[i]=f>=0?a.strides[f]:0;}
}
static void recurse_assign(symbolic_tensor& r,const symbolic_tensor& a,const vector<int>& as,int i,int ir,int ia)
{
    if(i>=r.size.size()){
        r.x[ir]=a.x[ia];
        return;}
    for(int s=0;s<r.size[i];s++)
        recurse_assign(r,a,as,i+1,ir+s*r.strides[i],ia+s*as[i]);
}
void symbolic_tensor::set(const std::string& ind,const symbolic_tensor& a)
{
    indices=ind;
    size.resize(indices.size());
    for(int i=0;i<(int)indices.size();i++)
        size[i]=a.size[a.indices.find(indices[i])];
    compute_strides();
    vector<int> as;
    reordered_strides(as,a,indices);
    recurse_assign(*this,a,as,0,0,0);
}
void symbolic_tensor::set_id(const std::string& ind,int new_size)
{
    indices=ind;
    size={new_size,new_size};
    compute_strides();
    for(int i=0;i<new_size;i++)
        (*this)(i,i)=1;
}
void symbolic_tensor::set_perm(const std::string& ind)
{
    indices=ind;
    size={3,3,3};
    compute_strides();
    (*this)(0,1,2)=1;
    (*this)(2,0,1)=1;
    (*this)(1,2,0)=1;
    (*this)(0,2,1)=-1;
    (*this)(2,1,0)=-1;
    (*this)(1,0,2)=-1;
}
void symbolic_tensor::compute_strides()
{
    if(size.size()){
        strides.resize(size.size());
        strides.back()=1;
        for(int i=strides.size()-1;i>0;i--)
            strides[i-1]=strides[i]*size[i];
        x.resize(strides[0]*size[0]);}
    else x.resize(1);
    fill(x.begin(),x.end(),0);
}
template<class f> static void recurse_op(symbolic_tensor& r,const symbolic_tensor& a,const symbolic_tensor& b,
    const vector<int>& as,const vector<int>& bs,f func,int i,int ir,int ia,int ib)
{
    if(i>=r.size.size()){
        r.x[ir]=func(a.x[ia],b.x[ib]);
        return;}
    for(int s=0;s<r.size[i];s++)
        recurse_op(r,a,b,as,bs,func,i+1,ir+s*r.strides[i],ia+s*as[i],ib+s*bs[i]);
}
symbolic_tensor operator+ (const symbolic_tensor& a,const symbolic_tensor& b)
{
    symbolic_tensor r(a);
    vector<int> bs;
    reordered_strides(bs,b,r.indices);
    recurse_op(r,a,b,a.strides,bs,[](double u,double v){return u+v;},0,0,0,0);
    return r;
}
symbolic_tensor operator- (const symbolic_tensor& a,const symbolic_tensor& b)
{
    symbolic_tensor r(a);
    vector<int> bs;
    reordered_strides(bs,b,r.indices);
    recurse_op(r,a,b,a.strides,bs,[](double u,double v){return u-v;},0,0,0,0);
    return r;
}
symbolic_tensor operator* (const symbolic_tensor& a,double b)
{
    symbolic_tensor r(a);
    for(auto& it:r.x) it*=b;
    return r;
}
static double recurse_mul(const symbolic_tensor& a,const symbolic_tensor& b,
    const std::string& sum,const vector<int>& sum_size,const vector<int>& as_sum,
    const vector<int>& bs_sum,int i,int ia,int ib)
{
    if(i>=(int)sum.size()) return a.x[ia]*b.x[ib];
    double x=0;
    for(int s=0;s<sum_size[i];s++)
        x+=recurse_mul(a,b,sum,sum_size,as_sum,bs_sum,i+1,ia+s*as_sum[i],ib+s*bs_sum[i]);
    return x;
}
static void recurse_mul(symbolic_tensor& r,const symbolic_tensor& a,const symbolic_tensor& b,
    const std::string& sum,const vector<int>& sum_size,const vector<int>& as,const vector<int>& bs,
    const vector<int>& as_sum,const vector<int>& bs_sum,int i,int ir,int ia,int ib)
{
    if(i>=r.size.size()){
        r.x[ir]=recurse_mul(a,b,sum,sum_size,as_sum,bs_sum,0,ia,ib);
        return;}
    for(int s=0;s<r.size[i];s++)
        recurse_mul(r,a,b,sum,sum_size,as,bs,as_sum,bs_sum,i+1,ir+s*r.strides[i],ia+s*as[i],ib+s*bs[i]);
}
symbolic_tensor operator* (const symbolic_tensor& a,const symbolic_tensor& b)
{
    symbolic_tensor r;
    int mp[26]={};
    for(int i=0;i<26;i++) mp[i]=-2;
    std::string sum;
    vector<int> sum_size;
    for(int i=0;i<(int)a.indices.size();i++)
        mp[a.indices[i]-'a']=i;
    for(int i=0;i<(int)b.indices.size();i++){
        int& x=mp[b.indices[i]-'a'];
        if(x!=-2){
            sum+=b.indices[i];
            sum_size.push_back(b.size[i]);
            x=-1;}
        else x=i;}
    for(int i=0;i<(int)a.indices.size();i++){
        int& x=mp[a.indices[i]-'a'];
        if(x<0) continue;
        r.indices+=a.indices[i];
        r.size.push_back(a.size[i]);}
    for(int i=0;i<(int)b.indices.size();i++){
        int& x=mp[b.indices[i]-'a'];
        if(x<0) continue;
        r.indices+=b.indices[i];
        r.size.push_back(b.size[i]);}
    r.compute_strides();
    vector<int> as,bs,as_sum,bs_sum;
    reordered_strides(as,a,r.indices);
    reordered_strides(bs,b,r.indices);
    reordered_strides(as_sum,a,sum);
    reordered_strides(bs_sum,b,sum);
    recurse_mul(r,a,b,sum,sum_size,as,bs,as_sum,bs_sum,0,0,0,0);
    return r;
}
symbolic_tensor operator+ (const symbolic_tensor& a,double b)
{
    assert(a.size.size()==0 && a.x.size()==1);
    symbolic_tensor r(a);
    r.x[0]+=b;
    return r;
}
symbolic_tensor operator/ (double b,const symbolic_tensor& a)
{
    assert(a.size.size()==0 && a.x.size()==1);
    symbolic_tensor r(a);
    r.x[0]=b/r.x[0];
    return r;
}
static double recurse_contract(const symbolic_tensor& a,const std::string& sum,const vector<int>& sum_size,
    const vector<int>& as_sum,int i,int ia)
{
    if(i>=(int)sum.size()) return a.x[ia];
    double x=0;
    for(int s=0;s<sum_size[i];s++)
        x+=recurse_contract(a,sum,sum_size,as_sum,i+1,ia+s*as_sum[i]);
    return x;
}
static void recurse_contract(symbolic_tensor& r,const symbolic_tensor& a,const std::string& sum,
    const vector<int>& sum_size,const vector<int>& as,const vector<int>& as_sum,int i,int ir,int ia)
{
    if(i>=r.size.size()){
        r.x[ir]=recurse_contract(a,sum,sum_size,as_sum,0,ia);
        return;}
    for(int s=0;s<r.size[i];s++)
        recurse_contract(r,a,sum,sum_size,as,as_sum,i+1,ir+s*r.strides[i],ia+s*as[i]);
}
symbolic_tensor symbolic_tensor::operator()(const std::string& new_indices) const
{
    assert(indices.size()==new_indices.size());
    int mp[26]={};
    for(int i=0;i<26;i++) mp[i]=-2;
    std::string sum,reduced_indices;
    vector<int> sum_size,sum_strides,new_size,as,as_sum;
    for(int i=0;i<(int)new_indices.size();i++){
        assert(islower(new_indices[i]));
        int& x=mp[new_indices[i]-'a'];
        if(x!=-2){
            assert(x>=0 && size[i]==size[x]);
            sum+=new_indices[i];
            sum_size.push_back(size[i]);
            sum_strides.push_back(strides[i]+strides[x]);
            x=-1;}
        else x=i;}
    for(int i=0;i<(int)new_indices.size();i++){
        int& x=mp[new_indices[i]-'a'];
        if(x<0) continue;
        reduced_indices+=new_indices[i];
        new_size.push_back(size[i]);}
    if(reduced_indices==new_indices){
        symbolic_tensor r(*this);
        r.indices=new_indices;
        return r;}
    symbolic_tensor r(reduced_indices,new_size);
    reordered_strides(as,*this,r.indices);
    recurse_contract(r,*this,sum,sum_size,as,sum_strides,0,0,0);
    return r;
}
void symbolic_tensor::random(random_type& rand,double lo,double hi)
{
    uniform_real_distribution<> dist(lo,hi);
    for(size_t i=0;i<x.size();i++)
        x[i]=dist(rand);
}
symbolic_tensor operator/ (const symbolic_tensor& a,const symbolic_tensor& b)
{
    assert(b.size.size()==0);
    return a/b.x[0];
}
symbolic_tensor scalar_op(const symbolic_tensor& a,function<double(double)> f)
{
    assert(a.size.size()==0);
    symbolic_tensor r(a);
    r.x[0]=f(a.x[0]);
    return r;
}
