//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include "CANONICAL_BLOCK.h"

namespace PhysBAM{

template<class T> void CANONICAL_BLOCK<T>::Compute_Element_Edges()
{
    element_edges.Resize(E.m);
    HASHTABLE<IV,PAIR<int,bool> > h;
    for(int i=0;i<S.m;i++)
    {
        h.Set(S(i),{i,0});
        h.Set(S(i).Reversed(),{i,1});
    }
    for(int i=0;i<E.m;i++)
        for(int j=0;j<3;j++)
            element_edges(i)(j)=h.Get({E(i)((j+1)%3),E(i)((j+2)%3)});
}
template class CANONICAL_BLOCK<double>;
}
