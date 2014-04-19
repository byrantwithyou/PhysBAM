//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG_PRINTF.h>
namespace PhysBAM{
namespace LOG_REAL{

void fprintf_parse_flags(const char *format,int len,PRINTF_FORMAT_FLAGS& flags)
{
    flags.width=0;
    flags.precision=INT_MAX;
    flags.just_right=true;
    char *end=0;
    const char *last=format+len;
    format++;
    while(format<last-1){
        switch(*format){
            case '0':
            case '#':
            case '+':
            case '\'':
            case 'I':
            case 'h':
            case 'l':
            case 'L':
            case 'q':
            case 'j':
            case 'z':
            case 't':
            case ' ':
                format++;
                continue;

            case '-':
                flags.just_right=false;
                format++;
                continue;

            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                flags.width=strtol(format,&end,10);
                format=end;
                continue;

            case '*':
                flags.width=-1;
                format++;
                continue;

            case '.':
                if(format[1]=='*'){
                    flags.precision=-1;
                    format+=2;}
                else{
                    flags.precision=strtol(format+1,&end,10);
                    format=end;}
                continue;

            default:
                throw std::runtime_error(std::string("invalid format string: ")+format);}}
}

int fprintf_rewrite_flags(const char *format,char* new_format,int len,const PRINTF_FORMAT_FLAGS& flags)
{
    char *start=new_format;
    const char *last=format+len;
    *new_format++=*format++;
    while(format<last){
        if(*format=='*'){
            new_format+=::std::sprintf(new_format,"%d",flags.width);
            format++;}
        else if(*format=='.' && format[1]=='*'){
            new_format+=::std::sprintf(new_format,".%d",flags.precision);
            format+=2;}
        else *new_format++=*format++;}
    *new_format=0;
    return new_format-start;
}

template<typename T>
int fprintf_formatted_item_builtin(std::ostream& out,const char *format,int len,T value)
{
    const int buff_size=1000;
    char buff[buff_size];
    int r=snprintf(buff,buff_size,format,value);
    if(r>=buff_size){
        char buff2[r+1];
        int s=snprintf(buff2,r+1,format,value);
        assert(s<=r);
        out<<buff2;
        return s;}
    out<<buff;
    return r;
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,int value)
{
    format[len-1]='d';
    return fprintf_formatted_item_builtin(out,format,len,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned int value)
{
    format[len-1]='u';
    return fprintf_formatted_item_builtin(out,format,len,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,long value)
{
    format[len-1]='l';
    format[len]='d';
    format[len+1]=0;
    return fprintf_formatted_item_builtin(out,format,len+1,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned long value)
{
    format[len-1]='l';
    format[len]='u';
    format[len+1]=0;
    return fprintf_formatted_item_builtin(out,format,len+1,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,long long value)
{
    format[len-1]='l';
    format[len]='l';
    format[len+1]='d';
    format[len+2]=0;
    return fprintf_formatted_item_builtin(out,format,len+2,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned long long value)
{
    format[len-1]='l';
    format[len]='l';
    format[len+1]='u';
    format[len+2]=0;
    return fprintf_formatted_item_builtin(out,format,len+2,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,short value)
{
    format[len-1]='h';
    format[len]='d';
    format[len+1]=0;
    return fprintf_formatted_item_builtin(out,format,len+1,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned short value)
{
    format[len-1]='h';
    format[len]='u';
    format[len+1]=0;
    return fprintf_formatted_item_builtin(out,format,len+1,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,char value)
{
    format[len-1]='h';
    format[len]='h';
    format[len+1]='d';
    format[len+2]=0;
    return fprintf_formatted_item_builtin(out,format,len+2,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,unsigned char value)
{
    format[len-1]='h';
    format[len]='h';
    format[len+1]='u';
    format[len+2]=0;
    return fprintf_formatted_item_builtin(out,format,len+2,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,double value)
{
    format[len-1]='g';
    return fprintf_formatted_item_builtin(out,format,len,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,long double value)
{
    format[len-1]='L';
    format[len]='g';
    format[len+1]=0;
    return fprintf_formatted_item_builtin(out,format,len+1,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,const char* value)
{
    format[len-1]='s';
    return fprintf_formatted_item_builtin(out,format,len,value);
}

int fprintf_formatted_item(std::ostream& out,char *format,int len,const PRINTF_FORMAT_FLAGS& flags,const void* value)
{
    format[len-1]='p';
    return fprintf_formatted_item_builtin(out,format,len,value);
}

int fprintf_formatted_item_string(std::ostream& out,const PRINTF_FORMAT_FLAGS& flags,const std::string& str)
{
    if(flags.width>(int)str.size() && flags.just_right)
        for(int i=str.size();i<flags.width;i++)
            out<<' ';
    out<<str;
    if(flags.width>(int)str.size() && !flags.just_right)
        for(int i=str.size();i<flags.width;i++)
            out<<' ';
    return std::max((int)str.size(),flags.width);
}

int fprintf(std::ostream& out,const char *format)
{ 
    while(*format){
        if(*format=='%'){
            if(format[1]=='%') format++;
            else throw std::runtime_error("invalid format string: missing arguments");}
        out<<*format++;}
    out<<std::flush;
    return 0;
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
template int fprintf_formatted_item_builtin<float>(std::ostream&,char const*,int,float);
template int fprintf_formatted_item_builtin<double>(std::ostream&,char const*,int,double);
template int fprintf_formatted_item_builtin<int>(std::ostream&,char const*,int,int);
template int fprintf_formatted_item_builtin<char>(std::ostream&,char const*,int,char);
template int fprintf_formatted_item_builtin<long>(std::ostream&,char const*,int,long);
template int fprintf_formatted_item_builtin<long long>(std::ostream&,char const*,int,long long);
template int fprintf_formatted_item_builtin<short>(std::ostream&,char const*,int,short);
template int fprintf_formatted_item_builtin<unsigned int>(std::ostream&,char const*,int,unsigned int);
template int fprintf_formatted_item_builtin<unsigned char>(std::ostream&,char const*,int,unsigned char);
template int fprintf_formatted_item_builtin<unsigned long>(std::ostream&,char const*,int,unsigned long);
template int fprintf_formatted_item_builtin<unsigned long long>(std::ostream&,char const*,int,unsigned long long);
template int fprintf_formatted_item_builtin<unsigned short>(std::ostream&,char const*,int,unsigned short);
template int fprintf_formatted_item_builtin<bool>(std::ostream&,char const*,int,bool);
}
}
