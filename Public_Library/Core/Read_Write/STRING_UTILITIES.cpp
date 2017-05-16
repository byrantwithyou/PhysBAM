//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Read_Write/STRING_UTILITIES.h>
#include <cctype>
#include <cstdio>
namespace PhysBAM{
std::string Stripped_Whitespace(const std::string& str)
{
    int first=0,last=(int)str.length()-1;
    while(first<=last && isspace(str[first])) first++;
    while(last>=0 && isspace(str[last])) last--;
    return str.substr(first,last-first+1);
}

bool Is_Number(const std::string& str) // integer or floating point
{
    const char *nptr=str.c_str(),*last_character=nptr+str.size();
    char *endptr;
    if(strtod(nptr,&endptr)){}
    return endptr==last_character;
}

int Compare_Strings(const char* str1,int len1,const char* str2,int len2,bool case_sensitive)
{
    int min_len=std::min(len1,len2);
    if(!case_sensitive){
        for(int i=0;i<min_len;i++)
            if(int r=::tolower(str1[i])-::tolower(str2[i]))
                return r<0?-1:1;}
    else if(int r=memcmp(str1,str2,min_len)) return r;
    if(len1<len2) return -1;
    return len1>len2;
}

int Compare_Strings(const std::string &str1,const std::string &str2,bool case_sensitive)
{
    return Compare_Strings(str1.c_str(),str1.length(),str2.c_str(),str2.length(),case_sensitive);
}

std::string toupper(const std::string& str)
{
    std::string str_copy=str;
    for(std::string::iterator iter=str_copy.begin();iter!=str_copy.end();iter++)
        *iter=::toupper(*iter);
    return str_copy;
}

bool Parse_Integer_Range(const std::string& str,ARRAY<int>& integer_list)
{
    const char* c_str=str.c_str();
    const char* endptr=c_str;
    int start_val=(int)strtol(c_str,&const_cast<char*&>(endptr),10);
    if(endptr==c_str) return false;
    while(isspace(*endptr)) endptr++; // skip whitespace
    if(*endptr=='\0'){
        integer_list.Append(start_val);
        return true;} // single value
    if(*endptr!='-') return false;
    const char* nextptr=endptr+1;
    int end_val=(int)strtol(nextptr,&const_cast<char*&>(endptr),10);
    if(endptr==nextptr) return false;
    while(isspace(*endptr)) endptr++; // skip whitespace
    int step;
    if(*endptr=='\0')
        step=1;
    else if(*endptr==':'){
        nextptr=endptr+1;
        step=(int)strtol(nextptr,&const_cast<char*&>(endptr),10);
        if(endptr==nextptr) return false;
        while(isspace(*endptr))endptr++; // skip whitespace
        if(*endptr) return false;}
    else return false;
    for(int i=start_val;i<=end_val;i+=step) integer_list.Append(i);
    return true;
}

// integer list format: [range],[range],[range],... where each [range] is either a single <number> or <number>-<number>, e.g. "1-3,4,7,10-11"
bool Parse_Integer_List(const std::string& str,ARRAY<int>& integer_list)
{
    integer_list.Remove_All();
    std::string remaining_string=str;
    while(!str.empty()){
        std::string token;
        std::string::size_type comma_pos=remaining_string.find(",");
        if(comma_pos!=std::string::npos){
            token=remaining_string.substr(0,comma_pos);
            remaining_string=remaining_string.substr(comma_pos+1);}
        else{
            token=remaining_string;
            remaining_string="";}
        if(!Parse_Integer_Range(token,integer_list))
            return false;}
    return true;
}

std::string Join(const std::string& separator,const ARRAY<std::string>& tokens)
{
    if(tokens.m==0) return "";
    std::string str=tokens(0);
    for(int i=1;i<tokens.m;i++) str+=separator+tokens(i);
    return str;
}

void Split(const std::string& str,const std::string& separator,ARRAY<std::string>& tokens)
{
    std::string remaining_string=str;
    std::string::size_type separator_length=separator.size();
    while(!remaining_string.empty()){
        std::string token;
        std::string::size_type separator_position=remaining_string.find(separator);
        if(separator_position!=std::string::npos){
            tokens.Append(remaining_string.substr(0,separator_position));
            remaining_string=remaining_string.substr(separator_position+separator_length);}
        else{
            tokens.Append(remaining_string);
            remaining_string="";}}
}

bool Ends_With(const std::string& input,const std::string& test)
{
    int input_index=input.length(),test_index=test.length();
    for(;test_index>0 && input_index>0;)
        if(input[--input_index]!=test[--test_index])
            return false;
    return test_index==0;
}

bool IEnds_With(const std::string& input,const std::string& test)
{
    int input_index=input.length(),test_index=test.length();
    for(;test_index>0 && input_index>0;)
        if(std::toupper(input[--input_index])!=std::toupper(test[--test_index]))
            return false;
    return test_index==0;
}

//#####################################################################
}
