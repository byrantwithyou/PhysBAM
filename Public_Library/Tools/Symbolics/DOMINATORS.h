//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
// Class DOMINATORS
//#####################################################################
#ifndef __DOMINATORS__
#define __DOMINATORS__

#include <Tools/Arrays/ARRAY.h>
namespace PhysBAM{
struct DOMINATORS
{
    const ARRAY<ARRAY<int> > &in,&out;
    ARRAY<int> idom,parent,vertex;
    ARRAY<ARRAY<int> > frontier,children;

    DOMINATORS(const ARRAY<ARRAY<int> >& in,const ARRAY<ARRAY<int> >& out)
        :in(in),out(out)
    {}

    void Compute_Idom();
    void Compute_Frontier();

protected:
    void DFS(const ARRAY<ARRAY<int> >& out,ARRAY<int>& semi,int v);
    void Link(ARRAY<int>& ancestor,int u,int w) {ancestor(w)=u;}
    int Eval(const ARRAY<int>& semi,ARRAY<int>& ancestor,ARRAY<int>& label,int v);
    void Compress(const ARRAY<int>& semi,ARRAY<int>& ancestor,ARRAY<int>& label,int v);
};
}
#endif
