//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ELIMINATION_FEM__
#define __ELIMINATION_FEM__
#include "BLOCK_MATRIX.h"
#include "BLOCK_VECTOR.h"
#include "COMPONENT_LAYOUT_FEM.h"

namespace PhysBAM{

template<class T> struct CACHED_ELIMINATION_MATRIX;

template<class T>
struct ELIMINATION_FEM
{
    COMPONENT_LAYOUT_FEM<T>& cl;

    ELIMINATION_FEM(COMPONENT_LAYOUT_FEM<T>& cl);
    
    void Eliminate_Irregular_Blocks(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Eliminate_Non_Seperators(CACHED_ELIMINATION_MATRIX<T>& cem);
    void Eliminate_Strip(CACHED_ELIMINATION_MATRIX<T>& cem,const ARRAY<BLOCK_ID>& a);
    void Eliminate_Simple(CACHED_ELIMINATION_MATRIX<T>& cem,BLOCK_ID first,CON_ID con_id_source);
};

}
#endif
