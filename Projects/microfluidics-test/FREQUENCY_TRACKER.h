//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FREQUENCY_TRACKER__
#define __FREQUENCY_TRACKER__
#include <map>
#include <set>
#include <vector>

namespace PhysBAM{

struct FREQUENCY_TRACKER
{
    typedef std::pair<std::set<int>*,int> pr;
    typedef std::multimap<int,pr> mm;
    std::map<int,mm::iterator> id_map;
    mm freq_map;

    void Add(int id, int row);
    void Remove(int id, int row);
    int Most_Frequent(std::vector<int>& v) const;
    bool Empty() const;
    void Dump() const;
};

}
#endif
