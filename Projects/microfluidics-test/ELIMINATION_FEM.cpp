//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "CACHED_ELIMINATION_MATRIX.h"
#include "ELIMINATION_FEM.h"

namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class T> ELIMINATION_FEM<T>::
ELIMINATION_FEM(COMPONENT_LAYOUT_FEM<T>& cl)
    :cl(cl)
{
}
//#####################################################################
// Function Eliminate_Strip
//#####################################################################
template<class T> void ELIMINATION_FEM<T>::
Eliminate_Strip(CACHED_ELIMINATION_MATRIX<T>& cem,const ARRAY<BLOCK_ID>& a)
{
    for(int sep=1;sep-1<a.m;sep*=2)
        for(int i=sep-1;i<a.m;i+=sep*2)
            if(cem.valid_row(Value(a(i))))
                cem.Eliminate_Row(Value(a(i)));
}
//#####################################################################
// Function Eliminate_Strip
//#####################################################################
template<class T> void ELIMINATION_FEM<T>::
Eliminate_Simple(CACHED_ELIMINATION_MATRIX<T>& cem,BLOCK_ID first,CON_ID con_id_source)
{
    ARRAY<BLOCK_ID> list;
    BLOCK_ID b=first;
    CON_ID con_id=con_id_source;
    while(1)
    {
        const auto& bl=cl.blocks(b);
        if(bl.connections.m>CON_ID(2)) break;
        if(bl.flags&1) break;
        list.Append(b);
        if(bl.connections.m==CON_ID(1)) break;
        CON_ID o(1-Value(con_id));
        if(!bl.connections(o).is_regular) break;
        b=bl.connections(o).id;
        con_id=bl.connections(o).con_id;
    }
    Eliminate_Strip(cem,list);
}
//#####################################################################
// Function Eliminate_Rows
//#####################################################################
template<class T> void ELIMINATION_FEM<T>::
Eliminate_Irregular_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem)
{
    // Get rid of irregular connections first.
    for(const auto& ic:cl.irregular_connections)
    {
        BLOCK_ID prev(-1);
        ARRAY<BLOCK_ID> found;
        for(const auto& p:ic.edge_on)
        {
            if(p.b!=prev)
            {
                found.Append(p.b);
                prev=p.b;
            }
        }
        Eliminate_Strip(cem,found);
    }
}
//#####################################################################
// Function Eliminate_Irregular_Blocks
//#####################################################################
template<class T> void ELIMINATION_FEM<T>::
Eliminate_Non_Seperators(CACHED_ELIMINATION_MATRIX<T>& cem)
{
    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
    {
        if(!cem.valid_row(Value(b))) continue;
        const auto& bl=cl.blocks(b);
        if(bl.connections.m>CON_ID(2) || bl.flags&1)
        {
            for(const auto& c:bl.connections)
                if(c.is_regular)
                    Eliminate_Simple(cem,c.id,c.con_id);
            if(bl.connections.m>CON_ID(2))
                cem.Eliminate_Row(Value(b));
        }
    }
    for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        assert(!cem.valid_row(Value(b)) || cl.blocks(b).flags&1);
}

template class ELIMINATION_FEM<double>;
}
