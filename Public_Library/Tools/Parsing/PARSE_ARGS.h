//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARSE_ARGS  
//##################################################################### 
#ifndef __PARSE_ARGS__
#define __PARSE_ARGS__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

class PARSE_ARGS
{
public:
    struct OPTION
    {
        ~OPTION();
        std::string opt;
        std::string name;
        std::string desc;
        bool* found;
        bool found_value;
        void* store;
        bool (*store_func)(void*,const std::string&);
        void (*print_default_func)(const void*);
    };

    struct EXTRA
    {
        ~EXTRA();
        std::string name;
        std::string desc;
        bool required;
        bool exhaust;
        bool* found;
        void* store;
        bool (*store_func)(void*,const std::string&);
        void (*print_default_func)(const void*);
    };

    int argc;
    char** argv;
    std::string program_name;
    HASHTABLE<std::string,OPTION> options;
    ARRAY<EXTRA> extras;
    bool unclaimed_arguments;

    PARSE_ARGS(int argc_input,char** argv_input);
    PARSE_ARGS(const PARSE_ARGS&) = delete;
    void operator=(const PARSE_ARGS&) = delete;
    ~PARSE_ARGS();

    template<class T> static bool store_impl(void* store,const std::string& s)
    {std::stringstream ss(s);ss>>*(T*)store;return !ss.fail();}

    template<class T,int d> static bool store_vec_impl(void* store,const std::string& s)
    {std::stringstream ss(s);VECTOR<T,d>& v=*(VECTOR<T,d>*)store;
    for(int i=0;i<d;i++) ss>>v(i);
    return !ss.fail();}

    template<class T> static bool store_multi_impl(void* store,const std::string& s)
    {ARRAY<T>* a=(ARRAY<T>*)store;a->Append(T());return store_impl<T>(&a->Last(),s);}

    template<class T> static void print_default_impl(const void* store)
    {LOG::cerr<<*(T*)store;}

    template<class T> void Add(const std::string& arg_str,T* store,bool* found,const std::string& name,const std::string& desc)
    {
        OPTION o={arg_str,name,desc,found,true,store,&store_impl<T>,&print_default_impl<T>};
        options.Set(arg_str,o);
    }

    template<class T,int d> void Add(const std::string& arg_str,VECTOR<T,d>* store,bool* found,const std::string& name,const std::string& desc)
    {
        OPTION o={arg_str,name,desc,found,true,store,&store_vec_impl<T,d>,&print_default_impl<VECTOR<T,d> >};
        options.Set(arg_str,o);
    }

    template<class T> void Add(const std::string& arg_str,ARRAY<T>* store,bool* found,const std::string& name,const std::string& desc)
    {
        OPTION o={arg_str,name,desc,found,true,store,&store_multi_impl<T>,&print_default_impl<ARRAY<T> >};
        options.Set(arg_str,o);
    }

    template<class T> void Add(const std::string& arg_str,T* store,const std::string& name,const std::string& desc)
    {Add(arg_str,store,0,name,desc);}

    template<class T> void Extra_Optional(T* store,bool* found,const std::string& name,const std::string& desc)
    {
        EXTRA e={name,desc,false,false,found,store,&store_impl<T>,&print_default_impl<T>};
        extras.Append(e);
    }

    template<class T,int d> void Extra_Optional(VECTOR<T,d>* store,bool* found,const std::string& name,const std::string& desc)
    {
        EXTRA e={name,desc,false,false,found,store,&store_vec_impl<T,d>,&print_default_impl<VECTOR<T,d> >};
        extras.Append(e);
    }

    template<class T> void Extra_Optional(ARRAY<T>* store,bool* found,const std::string& name,const std::string& desc)
    {
        EXTRA e={name,desc,false,true,found,store,&store_multi_impl<T>,&print_default_impl<ARRAY<T> >};
        extras.Append(e);
    }

    template<class T> void Extra_Optional(T* store,const std::string& name,const std::string& desc)
    {Extra_Optional(store,0,name,desc);}

    template<class T> void Extra(T* store,const std::string& name,const std::string& desc)
    {
        Extra_Optional(store,0,name,desc);
        extras.Last().required=true;
    }

//#####################################################################
    void Add(const std::string& arg_str,bool* found,const std::string& desc);
    void Add_Not(const std::string& arg_str,bool* found,const std::string& desc);
    void Parse(bool partial=false);
    const std::string& Get_Program_Name() const;
    void Print_Usage(bool do_exit=false) const;
    std::string Print_Arguments() const;
//#####################################################################
};
extern template void PARSE_ARGS::print_default_impl<int>(const void*);
extern template void PARSE_ARGS::print_default_impl<std::string>(const void*);
extern template void PARSE_ARGS::print_default_impl<float>(const void*);
extern template void PARSE_ARGS::print_default_impl<double>(const void*);
extern template bool PARSE_ARGS::store_impl<int>(void*, const std::string&);
extern template bool PARSE_ARGS::store_impl<std::string>(void*, const std::string&);
extern template bool PARSE_ARGS::store_impl<float>(void*, const std::string&);
extern template bool PARSE_ARGS::store_impl<double>(void*, const std::string&);
extern template void PARSE_ARGS::Add<int>(const std::string&, int*, bool*, const std::string&, const std::string&);
extern template void PARSE_ARGS::Add<std::string>(const std::string&, std::string*, bool*, const std::string&, const std::string&);
extern template void PARSE_ARGS::Add<float>(const std::string&, float*, bool*, const std::string&, const std::string&);
extern template void PARSE_ARGS::Add<double>(const std::string&, double*, bool*, const std::string&, const std::string&);

extern template ARRAY<HASHTABLE_ENTRY_TEMPLATE<std::string, PARSE_ARGS::OPTION>, int>::~ARRAY();
extern template bool HASHTABLE<std::string, PARSE_ARGS::OPTION>::Set(const std::string&, const PARSE_ARGS::OPTION&);
extern template void HASHTABLE<std::string, PARSE_ARGS::OPTION>::Resize_Table(int);
}
#endif
