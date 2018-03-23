//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG_PRINTF.h>
namespace PhysBAM{
namespace LOG{

// 1 = need width, 2 = need precision, 4 = need precision for int, 8 = type n
int fprintf_parse_flags(std::ostream& out,const char *format,int len)
{
    int ret_flags=0;
    char *end=0;
    const char *last=format+len;

    bool is_int=false;
    switch(last[-1]){
        case 'd':
        case 'i':
        case 'u':
            is_int=true;
            break;
        case 's':
        case 'c':
        case 'p':
        case 'P':
            break;
        case 'F':
            out<<std::uppercase;
        case 'f':
            out.setf(std::ios_base::fixed);
            break;
        case 'E':
            out<<std::uppercase;
        case 'e':
            out.setf(std::ios_base::scientific);
            break;
        case 'G':
            out<<std::uppercase;
        case 'g':
            out<<std::defaultfloat;
            break;
        case 'X':
            out<<std::uppercase;
        case 'x':
            out<<std::hex;
            is_int=true;
            break;
        case 'o':
            out<<std::oct;
            is_int=true;
            break;
        case 'A':
            out<<std::uppercase;
        case 'a':
            out<<std::hexfloat;
            break;
        case 'n': ret_flags|=8; break;
            break;
        default:
            throw std::runtime_error(std::string("invalid format string: ")+format);}

    format++; // skip %

    // Parse "flags"
    bool done=false; 
    while(!done)
        switch(*format++){
            case '+': out.setf(std::ios_base::showpos); break;
            case '-': out.setf(std::ios_base::left); break;
            case ' ': break; // TODO what do we do here?
            case '0': out<<std::setfill('0'); break;
            case '#':
                out<<std::showpoint;
                out<<std::showbase;
            // switch(last[-1]){
            //     case 'g':
            //     case 'G':
            //         out<<std::showpoint;
            //         break; // TODO what do we do here?   (For g and G types, trailing zeros are not removed.)
            //     case 'f':
            //     case 'F':
            //     case 'e':
            //     case 'E':
            //         out<<std::showpoint;
            //         break;
            //     case 'o':
            //     case 'x':
            //     case 'X':
            //         out<<std::showbase;
            //         break;
            // }
            break;
            default: done=true;break;}
    format--;

    // Parse width
    if(isdigit(*format)){
        out<<std::setw(strtol(format,&end,10));
        format=end;}
    else if(*format=='*'){
        ret_flags|=1;
        format++;}

    // Parse precision
    if(*format=='.'){
        if(format[1]=='*'){
            ret_flags|=2<<is_int;
            format+=2;}
        else if(format[1]=='-' || isdigit(format[1])){
            int p=strtol(format+1,&end,10);
            if(p>=0){
                if(is_int){
                    ret_flags|=32;
                    out<<std::setw(std::max((int)out.width(),p));
                    out<<std::setfill('0');}
                else out<<std::setprecision(p);}
            format=end;}
        else{
            out<<std::setprecision(0);
            format++;}}

    // Parse length
    switch(*format){
        case 'h':
        case 'l':
            format++;
        case 'L':
        case 'z':
        case 'j':
        case 't':
            format++;}

    if(format!=last-1)
        throw std::runtime_error(std::string("invalid format string: ")+format);

    return ret_flags;
}

void fprintf_rec(const char *format,OUT_STATE& state)
{ 
    while(*format){
        if(*format=='%'){
            if(format[1]=='%') format++;
            else throw std::runtime_error("invalid format string: missing arguments");}
        state.out<<*format++;}
}

int fputc(int c, std::ostream& out)
{
    if(out<<c) return (unsigned char)c;
    return EOF;
}

int fputs(const char *s, std::ostream& out)
{
    if(out<<s) return 0;
    return EOF;
}

int puts(const char *s)
{
    if(LOG::cout<<s<<std::endl) return 0;
    return EOF;
}
}
}
