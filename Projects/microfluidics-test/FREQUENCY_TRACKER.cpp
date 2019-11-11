//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "FREQUENCY_TRACKER.h"
#include <cassert>
namespace PhysBAM{
void FREQUENCY_TRACKER::Add(int id, int row)
{
    auto it=id_map.find(id);
    pr p;
    if(it!=id_map.end())
    {
        p=it->second->second;
        freq_map.erase(it->second);
    }
    else p=std::make_pair(new std::set<int>,id);
    p.first->insert(row);
    id_map[id]=freq_map.insert(std::make_pair(p.first->size(),p));
}
void FREQUENCY_TRACKER::Remove(int id, int row)
{
    auto it=id_map.find(id);
    if(it==id_map.end()) return;
    pr p=it->second->second;
    freq_map.erase(it->second);
    p.first->erase(row);
    if(p.first->empty())
    {
        delete p.first;
        id_map.erase(it);
    }
    else id_map[id]=freq_map.insert(std::make_pair(p.first->size(),p));
}
int FREQUENCY_TRACKER::Most_Frequent(std::vector<int>& v) const
{
    v.clear();
    auto it=freq_map.end();
    it--;
    for(auto a:*it->second.first)
        v.push_back(a);
    return it->second.second;
}
bool FREQUENCY_TRACKER::Empty() const
{
    return freq_map.empty();
}
void FREQUENCY_TRACKER::Dump() const
{
    auto dump_mm=[](mm::const_iterator m)
        {
            printf("freq: %i  rows:",m->first);
            for(auto a:*m->second.first) printf(" %i",a);
            printf(" id: %i",m->second.second);
        };

    printf("id map:\n");
    for(auto a:id_map)
    {
        printf("%i -> [",a.first);
        dump_mm(a.second);
        printf("]\n");
    }

    printf("freq map:\n");
    for(auto it=freq_map.begin();it!=freq_map.end();it++)
    {
        dump_mm(it);
    printf("\n");
    }

}
}
