//#####################################################################
// Copyright 2004-2006, Zhaosheng Bao, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KD_TREE_NODE
//#####################################################################
#ifndef __KD_TREE_NODE__
#define __KD_TREE_NODE__

namespace PhysBAM{

template<class T>
class KD_TREE_NODE
{
public:
    int split_axis; // -1 means leaf
    T split_value;
    int node_index;
    KD_TREE_NODE<T>* left;
    KD_TREE_NODE<T>* right;

    KD_TREE_NODE()
        :split_axis(-1),split_value(0),left(0),right(0)
    {}

    KD_TREE_NODE(const KD_TREE_NODE&) = delete;
    void operator=(const KD_TREE_NODE&) = delete;

    ~KD_TREE_NODE()
    {}

//#####################################################################
};
}
#endif
