//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MANAGED_MATRIX__
#define __MANAGED_MATRIX__
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Vectors/VECTOR.h>
#include "COMMON.h"
#include <atomic>
#include <functional>
#include <mutex>

namespace PhysBAM{

template<class T,class OBJ>
struct MANAGED_MATRIX
{
    enum class state_type {need_fill, in_cache, on_disk, freed};

    std::string pattern;

    struct ENTRY
    {
        OBJ* obj=0;
        std::function<void(OBJ& M)> fill_func=0;
        int num_users=0;
        bool recent=0;
        state_type state=state_type::need_fill;
        std::mutex* mtx=0;
    };

    ARRAY<ENTRY> entries;
    int entries_size=0;
    // alias state on_disk to need_fill
    bool refill_on_pagefault=false;

    std::mutex clock_mtx;
    ARRAY<int> cache_entries;
    int head=0;

    MANAGED_MATRIX()
    {
    }

    ~MANAGED_MATRIX()
    {
        for(auto& e:entries) delete e.mtx;
    }

    void Init(const std::string& pat,int cache_size=-1)
    {
        if(cache_size<0) cache_size=entries.m;
        pattern=pat;
        cache_entries.Resize(cache_size,use_init,-1);
        for(auto& e:entries) e.mtx=new std::mutex;
    }

    OBJ& Use(int i)
    {
        ENTRY& e=entries(i);
        std::unique_lock<std::mutex> lck(*e.mtx);
        if(e.state==state_type::in_cache)
        {
            e.recent=true;
            e.num_users++;
            e.state=state_type::in_cache;
            return *e.obj;
        }

        Reserve_Cache_Entry(i);
        Load(i);
        e.recent=false;
        e.num_users++;
        return *e.obj;
    }

    void Reserve_Cache_Entry(int k)
    {
        std::unique_lock<std::mutex> clock_lck(clock_mtx);
        while(1)
        {
            int &i=cache_entries(head);
            if(++head==cache_entries.m)
                head=0;

            if(i<0)
            {
                i=k;
                return;
            }

            ENTRY& e=entries(i);
            std::unique_lock<std::mutex> lck(*e.mtx,std::try_to_lock);
            if(lck.owns_lock())
            {
                if(e.state==state_type::freed)
                {
                    i=k;
                    return;
                }
                assert(e.state==state_type::in_cache);
                if(e.recent)
                {
                    e.recent=false;
                    continue;
                }
                if(e.num_users) continue;
                int old_i=i;
                i=k;
                clock_lck.unlock();

                Evict(old_i);
                return;
            }
        }
    }

    void Release(int i, bool remove)
    {
        ENTRY& e=entries(i);
        std::unique_lock<std::mutex> lck(*e.mtx);
        assert(e.state==state_type::in_cache);
        --e.num_users;
        if(remove)
        {
            assert(!e.num_users);
            delete e.obj;
            e.obj=0;
            e.state=state_type::freed;
            e.recent=false;
        }
    }

    // must own mtx to call
    void Load(int i)
    {
        ENTRY& e=entries(i);
        e.obj=new OBJ;
        if(e.state==state_type::need_fill)
        {
            if(e.fill_func) e.fill_func(*e.obj);
        }
        else
        {
            assert(e.state==state_type::on_disk);
            if(!refill_on_pagefault)
            {
                std::string file=LOG::sprintf(pattern.c_str(),i);
                Read_From_File(file,*e.obj);
                Remove_File(file);
            }
            else
            {
                if(e.fill_func) e.fill_func(*e.obj);
            }
        }
        e.state=state_type::in_cache;
    }

    // must own mtx to call
    void Evict(int i)
    {
        ENTRY& e=entries(i);
        if(!refill_on_pagefault)
        {
            std::string file=LOG::sprintf(pattern.c_str(),i);
            Write_To_File<T>(file,*e.obj);
        }
        delete e.obj;
        e.obj=0;
        e.state=state_type::on_disk;
    }
};
}
#endif
