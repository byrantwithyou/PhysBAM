//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINTER_POOL
//##################################################################### 
#ifndef __POINTER_POOL__
#define __POINTER_POOL__

#include <Core/Arrays/ARRAY.h>
namespace PhysBAM{

template<class T,int block_size=4096>
class POINTER_POOL
{
private:
    ARRAY<T*> pools;
    int current;
public:

    POINTER_POOL()
        :current(0)
    {pools.Append(new T[block_size]);}
    
    POINTER_POOL(const POINTER_POOL&) = delete;
    void operator=(const POINTER_POOL&) = delete;

    ~POINTER_POOL()
    {Delete_All();}

    T* New() // memory will be freed by Delete_All (do not use delete!)
    {if(current>=block_size){pools.Append(new T[block_size]);current=0;}
    return pools.Last()+(current++);}

    void Delete_All()
    {for(int i=0;i<pools.m;i++) delete[] pools(i);
    pools.Remove_All();current=block_size;}

//#####################################################################
};
}
#endif
