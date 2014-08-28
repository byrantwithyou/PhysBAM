//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOG_PRINTF
//##################################################################### 
#ifndef __LOG_PRINTF__
#define __LOG_PRINTF__
#include <Tools/Log/LOG.h>
#include <climits>
#include <cstdio>
namespace PhysBAM{
namespace LOG_REAL{

struct PRINTF_FORMAT_FLAGS
{
    int width,precision;
    bool just_right;
};

template<typename T>
int fprintf_formatted_item_builtin(std::ostream& out,const char *format,int len,T value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,int value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned int value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,long value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned long value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,long long value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned long long value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,short value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned short value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,char value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned char value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,double value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,long double value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,const char* value);
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,const void* value);

template<typename T>
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,const T* value)
{
    return fprintf_formatted_item(out,format,len,flags,(const void*)value);
}

template<typename T>
int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,T* value)
{
    return fprintf_formatted_item(out,format,len,flags,(const void*)value);
}

extern int fprintf_formatted_item_string(std::ostream& out,const PRINTF_FORMAT_FLAGS& flags,const std::string& str);

template<typename T> typename DISABLE_IF<IS_FUNDAMENTAL<typename REMOVE_REFERENCE<T>::TYPE>::value,int>::TYPE
fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,T&& value)
{
    std::ostringstream stream;
    if(flags.precision!=INT_MAX){
        int old_prec=stream.precision(flags.precision);
        stream<<value;
        stream.precision(old_prec);}
    else stream<<value;
    return fprintf_formatted_item_string(out,flags,stream.str());
}

extern void fprintf_parse_flags(const char *format,int len,PRINTF_FORMAT_FLAGS& flags);
extern int fprintf_rewrite_flags(const char *format,char* new_format,int len,const PRINTF_FORMAT_FLAGS& flags);

int fprintf(std::ostream& out,const char *format);
template<typename T,typename... Args> int fprintf(std::ostream& out,const char *format,T&& value,Args&&... args);

template<typename T> typename DISABLE_IF<IS_FUNDAMENTAL<typename REMOVE_REFERENCE<T>::TYPE>::value,int>::TYPE
fprintf_with_format(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,T&& value)
{
    return fprintf_formatted_item(out,format,len,flags,value);
}

template<typename T> typename ENABLE_IF<IS_FUNDAMENTAL<typename REMOVE_REFERENCE<T>::TYPE>::value,int>::TYPE
fprintf_with_format(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,T&& value)
{
    if(format[len-1]!='P')
        return fprintf_formatted_item_builtin(out,format,len,value);
    return fprintf_formatted_item(out,format,len,flags,value);
}

template<typename T,typename... Args> int
fprintf_fill_format(std::ostream& out,const char *format,int len,PRINTF_FORMAT_FLAGS& flags)
{
    throw std::runtime_error("invalid format string.");
}

template<typename T,typename... Args> int
fprintf_fill_format(std::ostream& out,const char *format,int len,PRINTF_FORMAT_FLAGS& flags,T&& param)
{
    throw std::runtime_error("invalid format string.");
}

template<typename T,typename U,typename... Args> typename DISABLE_IF<IS_INTEGRAL<typename REMOVE_REFERENCE<T>::TYPE>::value,int>::TYPE
fprintf_fill_format(std::ostream& out,const char *format,int len,PRINTF_FORMAT_FLAGS& flags,T&& param,
    U&& value,Args&&... args)
{
    throw std::runtime_error("invalid format string.");
}

template<typename T,typename... Args> int
fprintf_fill_format(std::ostream& out,const char *format,int len,PRINTF_FORMAT_FLAGS& flags,int param,
    T&& value,Args&&... args)
{
    if(flags.width<0){
        flags.width=param;
        if(flags.precision<0)
            return fprintf_fill_format(out,format,len,flags,value,args...);}
    else flags.precision=param;
    char new_format[len+50];
    int new_len=fprintf_rewrite_flags(format,new_format,len,flags);
    int n=fprintf_with_format(out,new_format,new_len,flags,value);
    return n+fprintf(out,format+len,args...);
}

template<typename T,typename... Args>
int fprintf(std::ostream& out,const char *format,T&& value,Args&&... args)
{
    int n=0;
    //while(*format){
    //    if(*format=='%'){
    //        if(format[1]=='%') format++;
    //        else{
    //            int option_len=strspn(format+1,"0123456789#-+ 'I*.hlLqjzt")+2;
    //            PRINTF_FORMAT_FLAGS flags;
    //            fprintf_parse_flags(format,option_len,flags);
    //            if(flags.width>=0 && flags.precision>=0){
    //                char term_format[option_len+3];
    //                memcpy(term_format,format,option_len);
    //                term_format[option_len]=0;
    //                n+=fprintf_with_format(out,term_format,option_len,flags,value);
    //                return n+fprintf(out,format+option_len,args...);}
    //            else
    //                return n+fprintf_fill_format(out,format,option_len,flags,value,args...);}}
    //    n++;
    //    out<<*format++;}
    return n;
}

template<typename... Args>
int printf(const char *format,Args&&... args)
{
    return fprintf(LOG::cout,format,args...);
}

template<typename... Args>
std::string sprintf(const char *format,Args&&... args)
{
    std::ostringstream stream;
    fprintf(stream,format,args...);
    return stream.str();
}
extern int fputc(int c, std::ostream& out);
extern int fputs(const char *s, std::ostream& out);
extern int puts(const char *s);
inline int putc(int c, std::ostream& out){return fputc(c,out);}
inline int putchar(int c){return putc(c,LOG::cout);}
}
}
#endif
