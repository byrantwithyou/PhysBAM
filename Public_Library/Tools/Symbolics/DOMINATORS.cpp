//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <Tools/Symbolics/DOMINATORS.h>
using namespace PhysBAM;
namespace PhysBAM{
//#####################################################################
// Function DFS
//#####################################################################
void DOMINATORS::
DFS(const ARRAY<ARRAY<int> >& out,ARRAY<int>& semi,int v)
{
    semi(v)=vertex.Append(v);
    for(int i=0;i<out(v).m;i++)
        if(semi(out(v)(i))==-1){
            parent(out(v)(i))=v;
            DFS(out,semi,out(v)(i));}
}
//#####################################################################
// Function Compress
//#####################################################################
void DOMINATORS::
Compress(const ARRAY<int>& semi,ARRAY<int>& ancestor,ARRAY<int>& label,int v)
{
    if(ancestor(ancestor(v))!=-1){
        Compress(semi,ancestor,label,ancestor(v));
        if(semi(label(ancestor(v)))<semi(label(v)))
            label(v)=label(ancestor(v));
        ancestor(v)=ancestor(ancestor(v));}
}
//#####################################################################
// Function Eval
//#####################################################################
int DOMINATORS::
Eval(const ARRAY<int>& semi,ARRAY<int>& ancestor,ARRAY<int>& label,int v)
{
    if(ancestor(v)==-1) return v;
    Compress(semi,ancestor,label,v);
    return label(v);
}
//#####################################################################
// Function Compute_Idom
//#####################################################################
void DOMINATORS::
Compute_Idom()
{
    ARRAY<int> ancestor(in.m),label(IDENTITY_ARRAY<>(in.m));
    ARRAY<int> semi(in.m);
    ARRAY<HASHTABLE<int> > bucket(in.m);
    idom.Resize(in.m);
    idom.Fill(-1);
    semi.Fill(-1);
    ancestor.Fill(-1);
    parent.Resize(in.m);

    DFS(out,semi,0);
    parent(0)=-1;

    for(int i=vertex.m-1;i>=1;i--){
        int w=vertex(i);
        for(int j=0;j<in(w).m;j++){
            int u=Eval(semi,ancestor,label,in(w)(j));
            if(semi(u)<semi(w)) semi(w)=semi(u);}
        bucket(vertex(semi(w))).Set(w);
        Link(ancestor,parent(w),w);
        for(HASHTABLE<int>::ITERATOR it(bucket(parent(w)));it.Valid();it.Next()){
            int u=Eval(semi,ancestor,label,it.Key());
            idom(it.Key())=(semi(u)<semi(it.Key()))?u:parent(w);}
        bucket(parent(w)).Remove_All();}

    for(int i=1;i<vertex.m;i++){
        int w=vertex(i);
        if(idom(w)!=vertex(semi(w)))
            idom(w)=idom(idom(w));}
    idom(vertex(0))=-1;

    children.Resize(out.m);
    for(int i=1;i<out.m;i++)
        children(parent(i)).Append(i);
}
//#####################################################################
// Function Compute_Frontier
//#####################################################################
void DOMINATORS::
Compute_Frontier()
{
    if(!parent.m) Compute_Idom();

    frontier.Resize(out.m);
    for(int i=vertex.m-1;i>=0;i--){
        HASHTABLE<int> hash;
        int x=vertex(i);
        for(int j=0;j<out(x).m;j++){
            int y=out(x)(j);
            if(idom(y)!=x)
                hash.Set(y);}
        for(int j=0;j<children(x).m;j++){
            int z=children(x)(j);
            for(int k=0;k<frontier(z).m;k++){
                int y=frontier(z)(k);
                if(idom(y)!=x)
                    hash.Set(y);}}
        hash.Get_Keys(frontier(x));}
}
}
