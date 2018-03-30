//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOG_PRINTF
//##################################################################### 
#ifndef __LOG_PRINTF__
#define __LOG_PRINTF__
#include <Core/Log/LOG.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <climits>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <iomanip>
#include <type_traits>
namespace PhysBAM{
namespace LOG{
using std::enable_if;
using std::is_integral;

struct OUT_STATE
{
    int width,precision,opts;
    std::ios::fmtflags flags;
    char fill;
    std::ostream::pos_type pos;
    std::ostream& out;
    
    OUT_STATE(std::ostream& out)
        :out(out)
    {
        width=out.width();
        precision=out.precision();
        opts=0;
        fill=out.fill();
        flags=out.flags();
        pos=out.tellp();
    }

    void Restore()
    {
        out.width(width);
        out.precision(precision);
        out.fill(fill);
        out.flags(flags);
    }

    int Number_Output()
    {
        return out.tellp()-pos;
    }
};

extern int fprintf_parse_flags(std::ostream& out,const char *format,int len);

void fprintf_rec(const char *format,OUT_STATE& state);
template<typename T,typename... Args> void fprintf_rec(const char *format,OUT_STATE& state,T&& value,Args&&... args);

template<typename T> typename enable_if<is_integral<T>::value && !is_const<T>::value>::type
capture_num_written(int num,T* value)
{
    *value=num;
}

template<typename T> void
capture_num_written(int num,T&& value)
{
    throw std::runtime_error("invalid format string.");
}

inline void
fprintf_fill_format_p(const char *format,OUT_STATE& state)
{
    throw std::runtime_error("invalid format string.");
}

template<typename T> void
fprintf_fill_format_p(const char *format,OUT_STATE& state,T&& param)
{
    throw std::runtime_error("invalid format string.");
}

template<typename T,typename U,typename... Args> typename enable_if<!is_integral<typename remove_reference<T>::type>::value>::type
fprintf_fill_format_p(const char *format,OUT_STATE& state,T&& param,
    U&& value,Args&&... args)
{
    throw std::runtime_error("invalid format string.");
}

template<typename T> void
fprintf_fill_format_p_safe(std::ostream& out,const T& value)
{
    out<<value;
}

inline void
fprintf_fill_format_p_safe(std::ostream& out,const char* value)
{
    if(value) out<<value;
    else out<<"(null)";
}

template<typename T,typename... Args> void
fprintf_fill_format_p(const char *format,OUT_STATE& state,int param,
    T&& value,Args&&... args)
{
    if(state.opts&4){
        state.out<<std::setw(std::max((int)state.out.width(),param));
        state.out<<std::setfill('0');}
    else state.out<<std::setprecision(param);
    int n=0;
    if(state.opts&8) capture_num_written(n,value);
    fprintf_fill_format_p_safe(state.out,value);
    state.Restore();
    fprintf_rec(format,state,args...);
}

inline void
fprintf_fill_format_w(const char *format,OUT_STATE& state)
{
    throw std::runtime_error("invalid format string.");
}

template<typename T> void
fprintf_fill_format_w(const char *format,OUT_STATE& state,T&& param)
{
    throw std::runtime_error("invalid format string.");
}

template<typename T,typename U,typename... Args> typename enable_if<!is_integral<typename remove_reference<T>::type>::value>::type
fprintf_fill_format_w(const char *format,OUT_STATE& state,T&& param,
    U&& value,Args&&... args)
{
    throw std::runtime_error("invalid format string.");
}

template<typename T,typename... Args> void
fprintf_fill_format_w(const char *format,OUT_STATE& state,int param,
    T&& value,Args&&... args)
{
    if(state.opts&32) state.out<<std::setw(std::max(param,(int)state.out.width()));
    else state.out<<std::setw(param);
    if(state.opts&6) return fprintf_fill_format_p(format,state,value,args...);
    int n=0;
    if(state.opts&8) capture_num_written(n,value);
    else fprintf_fill_format_p_safe(state.out,value);
    state.Restore();
    fprintf_rec(format,state,args...);
}

template<typename T,typename... Args>
void fprintf_rec(const char *format,OUT_STATE& state,T&& value,Args&&... args)
{
    while(*format){
        if(*format=='%'){
            if(format[1]=='%') format++;
            else{
                int option_len=strspn(format+1,"0123456789#-+ 'I*.hlLqjzt")+2;
                state.opts=fprintf_parse_flags(state.out,format,option_len);
                if(state.opts&1) return fprintf_fill_format_w(format+option_len,state,value,args...);
                if(state.opts&6) return fprintf_fill_format_p(format+option_len,state,value,args...);
                if(state.opts&8) capture_num_written(state.Number_Output(),value);
                else fprintf_fill_format_p_safe(state.out,value);
                state.Restore();
                return fprintf_rec(format+option_len,state,args...);}}
        state.out<<*format++;}
}

template<typename... Args>
int fprintf(std::ostream& out,const char *format,Args&&... args)
{
    OUT_STATE state(out);
    fprintf_rec(format,state,args...);
    out<<std::flush;
    return state.Number_Output();
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
