//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Flatten_Index
//#####################################################################
#ifndef __FLATTEN_INDEX__
#define __FLATTEN_INDEX__
namespace PhysBAM{
static int Flatten_Index(const VECTOR<int,2>& ind,const VECTOR<int,2>& count)
{
    return count(1)*ind(0)+ind(1);
}
static int Flatten_Index(const VECTOR<int,3>& ind,const VECTOR<int,3>& count)
{
    return count(1)*count(2)*ind(0)+count(2)*ind(1)+ind(2);
}
}
#endif
