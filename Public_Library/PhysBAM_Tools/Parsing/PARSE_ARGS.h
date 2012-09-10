//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARSE_ARGS  
//##################################################################### 
#ifndef __PARSE_ARGS__
#define __PARSE_ARGS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Parsing/ARG_DATA.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

class PARSE_ARGS:public NONCOPYABLE
{
private:
    ARRAY<ARG_DATA> arg_data_list;
    int num_expected_extra_args;
    std::string extra_args_synopsis,extra_args_desc,program_name;
    bool use_help_option;
    void (*extra_usage_callback)();
public:

    struct OPTION
    {
        std::string opt;
        std::string name;
        std::string desc;
        bool* found;
        bool found_value;
        void* store;
        bool (*store_func)(void*,const std::string&);
        void (*print_default_func)(const void*);
    };

    int argc;
    char** argv;
    HASHTABLE<std::string,OPTION> options;
    ARRAY<std::string> extra_arg_list;

    PARSE_ARGS(int argc_input,char** argv_input);
    ~PARSE_ARGS();

    template<class T> static bool store_impl(void* store,const std::string& s)
    {std::stringstream ss(s);return ss>>*(T*)store;}

    template<class T,int d> static bool store_vec_impl(void* store,const std::string& s)
    {std::stringstream ss(s);VECTOR<T,d>& v=*(VECTOR<T,d>*)store;for(int i=0;i<d;i++) ss>>v(i);return ss;}

    template<class T> static bool store_multi_impl(void* store,const std::string& s)
    {ARRAY<T>* a=(ARRAY<T>*)store;a->Append(T());return store_impl<T>(&a->Last(),s);}

    template<class T> static void print_default_impl(const void* store)
    {LOG::cerr<<*(T*)store;}

    template<class T> void Add(const std::string& arg_str,T* store,bool* found,const std::string& name,const std::string& desc)
    {
        OPTION o={arg_str,name,desc,found,true,store,&store_impl<T>,&print_default_impl<T>};
        options.Set(arg_str,o);
    }

    template<class T,int d> void Add(const std::string& arg_str,bool* found,VECTOR<T,d>* store,const std::string& name,const std::string& desc)
    {
        OPTION o={arg_str,name,desc,found,true,store,&store_vec_impl<T,d>,&print_default_impl<VECTOR<T,d> >};
        options.Set(arg_str,o);
    }

    template<class T> void Add(const std::string& arg_str,bool* found,ARRAY<T>* store,const std::string& name,const std::string& desc)
    {
        OPTION o={arg_str,name,desc,found,true,store,&store_multi_impl<T>,&print_default_impl<ARRAY<T> >};
        options.Set(arg_str,o);
    }

    void Add(const std::string& arg_str,bool* found,const std::string& desc)
    {Add(arg_str,(int*)0,found,"",desc);}

    void Add_Not(const std::string& arg_str,bool* found,const std::string& desc)
    {
        OPTION o={arg_str,"",desc,found,false,(int*)0,&store_impl<int>,&print_default_impl<int>};
        options.Set(arg_str,o);
    }

    template<class T> void Add(const std::string& arg_str,T* store,const std::string& name,const std::string& desc)
    {Add(arg_str,store,0,name,desc);}

//#####################################################################
    void Use_Help_Option(bool use_it);
    void Add_Option_Argument(const std::string& arg_str,const std::string& desc="");
    void Add_Integer_Argument(const std::string& arg_str,int default_value,const std::string& val_name="",const std::string& desc="");
    void Add_Double_Argument(const std::string& arg_str,double default_value,const std::string& val_name="",const std::string& desc="");
    void Add_String_Argument(const std::string& arg_str,const std::string& default_value,const std::string& val_name="",const std::string& desc="");
    void Set_Extra_Arguments(int num,const std::string& synopsis="",const std::string& desc="");
    void Set_Extra_Usage_Callback(void (*extra_usage_callback_input)());
    void Parse(bool partial=false);
    bool Get_Option_Value(const std::string& arg_str) const;
    int Get_Integer_Value(const std::string& arg_str) const;
    double Get_Double_Value(const std::string& arg_str) const;
    VECTOR<double,2> Get_Vector_2D_Value(const std::string& arg_str) const;
    VECTOR<double,3> Get_Vector_3D_Value(const std::string& arg_str) const;
    const std::string& Get_String_Value(const std::string& arg_str) const;
    bool Is_Value_Set(const std::string& arg_str) const;
    void Override_String_Value(const std::string& arg_str,const std::string& value);
    int Num_Extra_Args() const;
    const std::string& Extra_Arg(int i) const;
    const std::string& Get_Program_Name() const;
    int Find_Match(const std::string& str) const;
    int Find_Match(const std::string& str,const ARG_DATA::TYPE& type) const;
    static bool Find_And_Remove(const char *str,int& argc,char** argv);
    static int Find_And_Remove_Integer(const char *str,int& argc,char** argv);
    static double Find_And_Remove_Double(const char *str,int& argc,char** argv);
    void Print_Usage(bool do_exit=false) const;
    std::string Print_Arguments() const;
//#####################################################################
};
}
#endif
