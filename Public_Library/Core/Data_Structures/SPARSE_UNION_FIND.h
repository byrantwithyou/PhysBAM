//#####################################################################
// Copyright 2008, Don Hatch, Geoffrey Irving, Michael Lentine, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_UNION_FIND
//#####################################################################
#ifndef __SPARSE_UNION_FIND__
#define __SPARSE_UNION_FIND__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/CONSTANT_ARRAY.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class ID> // ID=int
class SPARSE_UNION_FIND
{
    typedef unsigned char T_RANK;
public:
    mutable HASHTABLE<ID,ID> parents;
    HASHTABLE<ID,T_RANK> ranks;
    ID size;

    explicit SPARSE_UNION_FIND(const ID entries=ID()) 
        :size(entries)
    {}

    void Initialize(const ID entries)
    {Clear_Connectivity();size=entries;}

    ID Size() const
    {return size;}

    void Clear_Connectivity()
    {parents.Remove_All();ranks.Remove_All();}

    ID Add_Entry()
    {return ++size;}

    bool Is_Root(const ID i) const
    {return !parents.Contains(i);}

    ID Find(const ID i) const
    {ID root=Find_Without_Path_Compression(i);Path_Compress(i,root);return root;}

    ID Union(const ID i,const ID j)
    {ID root_i=Find_Without_Path_Compression(i),root_j=Find_Without_Path_Compression(j);
    T_RANK rank_i(0),rank_j(0);ranks.Get(root_i,rank_i);ranks.Get(root_j,rank_j);
    ID root=rank_i>=rank_j?root_i:root_j;
    Path_Compress(i,root);Path_Compress(j,root);
    if(rank_i==rank_j && root_i!=root_j) ranks.Get_Or_Insert(root)++;
    return root;}

    template<class T_ARRAY>
    int Union(const T_ARRAY& array)
    {int root=-1;typename T_ARRAY::ELEMENT i(0);for(;i<array.Size();i++){root=Find(array(i));break;}if(root<0) return -1;
    for(;i<array.Size();i++) Union(root,array(i));
    return Find(root);}

    template<int d>
    ID Union(const VECTOR<ID,d>& indices)
    {T_RANK max_rank(0);ID root=Find_Without_Path_Compression(indices[0]);bool max_tie=false;ranks.Get(root,max_rank);
    for(int i=1;i<d;i++) {
        ID root_i=Find_Without_Path_Compression(indices[i]);
        T_RANK tmp_rank(0);ranks.Get(root_i,tmp_rank);
        if(max_rank<tmp_rank){max_rank=tmp_rank;root=root_i;max_tie=false;}
        else if(max_rank==tmp_rank && root!=root_i) max_tie=true;}
    for(int i=0;i<d;i++) Path_Compress(indices[i],root);
    if(max_tie) ranks.Get_Or_Insert(root)++;
    return root;}

    void Merge(const SPARSE_UNION_FIND<ID>& union_find)
    {assert(Size()==union_find.Size());for(const auto& it:union_find.parents) Union(it.Key(),it.Data());}

    // Okay for map to yield invalid indices for isolated elements
    template<class ID2,class T_ARRAY>
    void Mapped_Merge(const SPARSE_UNION_FIND<ID2>& union_find,const T_ARRAY& map)
    {for(const auto& it:union_find.parents) Union(map(it.Key()),map(it.Data()));}

    void Forest_Edges(ARRAY<PAIR<ID,ID> >& pairs) const
    {pairs.Remove_All();for(const auto& it:parents) pairs.Append(PAIR<ID,ID>(it.Key(),it.Data()));}

    void Merge_Forest_Edges(const ARRAY<PAIR<ID,ID> >& pairs)
    {for(int i=0;i<pairs.m;i++) Union(pairs(i).x,pairs(i).y);}

private:
    ID Find_Without_Path_Compression(const ID i) const
    {ID k,j=i;while(parents.Get(j,k)) j=k;return j;}

    void Path_Compress(const ID i,const ID root) const
    {ID j=i;while(ID *k=parents.Get_Pointer(j)){j=*k;*k=root;}if(j!=root) parents.Set(j,root);}

//#####################################################################
};
}
#endif
