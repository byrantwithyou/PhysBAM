//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Read_Write/ZIP.h>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <fstream>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif
#if defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__)
#include <sys/stat.h>
#endif
namespace PhysBAM{

//###################################################################
// Win32 Specific Function Definitions
//###################################################################
#if defined(_WIN32)

bool Directory_Exists(const std::string& dirname)
{
    DWORD attr=GetFileAttributes(dirname.c_str());
    return((attr!=-1)&&(attr&FILE_ATTRIBUTE_DIRECTORY));
}

bool Create_Directory(const std::string& dirname)
{
    if(!Directory_Exists(dirname)){
        LOG::cout<<"Creating directory using CreateDirectory...";
        CreateDirectory(dirname.c_str(),0);
        if(!Directory_Exists(dirname)){
            LOG::cerr<<"Failed!"<<std::endl;
            throw FILESYSTEM_ERROR("Create_Directory failed");
            return false;}
        LOG::cout<<"Successful!"<<std::endl;}
    return true;
}

bool Remove_File(const std::string& filename,bool check_compressed)
{
    PHYSBAM_NOT_IMPLEMENTED();
}

bool Remove_Directory(const std::string& filename)
{
    PHYSBAM_NOT_IMPLEMENTED();
}

std::string Real_Path(const std::string& path)//TODO: Implement this for windows
{
    PHYSBAM_NOT_IMPLEMENTED();
}

int Compare_File_Times_Ignoring_Compression_Suffix(const std::string& filename1,const std::string& filename2)
{
    HANDLE handle1=CreateFile(filename1.c_str(),0,0,0,OPEN_EXISTING,0,0);
    if(handle1==INVALID_HANDLE_VALUE){
        LOG::cerr<<"Compare_File_Times: can't CreateFile "<<filename1<<std::endl;
        return 0;}
    HANDLE handle2=CreateFile(filename2.c_str(),0,0,0,OPEN_EXISTING,0,0);
    if(handle2==INVALID_HANDLE_VALUE){
        LOG::cerr<<"Compare_File_Times: can't CreateFile "<<filename2<<std::endl;
        return 0;}
    FILETIME time1,time2;
    if(!GetFileTime(handle1,0,&time1,0)||!GetFileTime(handle2,0,&time2,0)){
        LOG::cerr<<"Compare_File_Times: error with GetFileTime"<<std::endl;
        return 0;}
    CloseHandle(handle1);
    CloseHandle(handle2);
    return CompareFileTime(&time1,&time2);
}

FILE* Temporary_File()
{
    return fopen(_tempnam(getenv("TEMP"),"pb"),"w");
}

std::string Get_Working_Directory()
{
    return ".";
}

//###################################################################
// Linux Specific Function Definitions
//###################################################################
#elif defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__)

bool Directory_Exists(const char* dirname)
{
    struct stat buf;
    int st_ret=stat(dirname,&buf);
    if(st_ret) return false;
    return S_ISDIR(buf.st_mode);
}

bool Directory_Exists(const std::string& dirname)
{
    return Directory_Exists(dirname.c_str());
}

bool Create_Directory(char* dirname)
{
    for(char* s=dirname+1; *s; s++)
        if(*s=='/')
        {
            *s=0;
            if(mkdir(dirname, 0775) && errno!=EEXIST) return false;
            *s='/';
        }

    if(!mkdir(dirname, 0775)) return true;
    if(errno!=EEXIST) return false;
    return Directory_Exists(dirname);
}

bool Create_Directory(const char* dirname)
{
    int len=strlen(dirname)+1;
    char copy[len];
    memcpy(copy,dirname,len);
    return Create_Directory(copy);
}

bool Create_Directory(const std::string& dirname)
{
    return Create_Directory(dirname.c_str());
}

bool Remove_File(const std::string& filename,bool check_compressed)
{
    if(!unlink(filename.c_str())) return true;
    if(check_compressed)
    {
        std::string comp=filename+".gz";
        return !unlink(comp.c_str());
    }
    return false;
}

bool Remove_Directory(const std::string& filename)
{
    return !rmdir(filename.c_str());
}

std::string Real_Path(const std::string& path)
{
    char real_path[PATH_MAX];
    PHYSBAM_ASSERT(realpath(path.c_str(),real_path));
    return real_path;
}

int Compare_File_Times_Ignoring_Compression_Suffix(const std::string& filename1,const std::string& filename2)
{
    struct stat stat1,stat2;
    if(stat(filename1.c_str(),&stat1)!=0) LOG::cerr<<"Compare_File_Times: can't stat "<<filename1<<std::endl;
    if(stat(filename2.c_str(),&stat2)!=0) LOG::cerr<<"Compare_File_Times: can't stat "<<filename2<<std::endl;
    if(stat1.st_mtime<stat2.st_mtime) return -1;
    else if(stat1.st_mtime>stat2.st_mtime) return 1;
    return 0;
}
FILE* Temporary_File()
{
    return tmpfile();
}

std::string Get_Working_Directory()
{
    ARRAY<char> buffer(128,no_init);
    for(;;){
        if(getcwd(buffer.Get_Array_Pointer(),buffer.m-1))
            return std::string(buffer.Get_Array_Pointer());
        else if(errno==ERANGE){
            if(buffer.m>=4096) PHYSBAM_FATAL_ERROR("refusing to allocate more than 4k to return working directory");
            buffer.Resize(2*buffer.m,no_init);}
        else PHYSBAM_FATAL_ERROR();}
}

//###################################################################
// Default (Unimplemented) Function Definitions
//###################################################################
#else

bool Directory_Exists(const std::string& dirname)
{
    return true;
}  // always return true on unsupported platforms

std::string Real_Path(const std::string& path)
{
    PHYSBAM_NOT_IMPLEMENTED();
}

int Compare_File_Times_Ignoring_Compression_Suffix(const std::string& filename1,const std::string& filename2)
{
    PHYSBAM_NOT_IMPLEMENTED();
}

FILE* Temporary_File()
{
    PHYSBAM_NOT_IMPLEMENTED();
}

std::string Get_Working_Directory()
{
    return ".";
}

#endif
//###################################################################
// File open with compression capability
//###################################################################
std::istream* Safe_Open_Input_Raw(const std::string& filename,bool binary)
{
    bool compressed=File_Is_Compressed(filename);
    std::ios_base::openmode flags=std::ios::in;if(binary) flags|=std::ios::binary;
    std::string filename_compressed=compressed?filename:filename+".gz";
#ifndef COMPILE_WITHOUT_ZLIB_SUPPORT
    if(File_Exists(filename_compressed))
        return Gzip_In(filename_compressed,flags);
#else
    if(File_Exists(filename_compressed))
        PHYSBAM_FATAL_ERROR("Must compile with zlib to read compressed files");
#endif
    if(!compressed){
        std::istream* input=new std::ifstream(filename.c_str(),flags);if(*input) return input;
        filename_compressed=filename+"(.gz)";}
    LOG::cerr<<"Can't open "<<filename_compressed<<" for read "<<(binary?"(binary)":"")<<std::endl;
    throw FILESYSTEM_ERROR("Safe_Open_Input failed");
}
std::ostream* Safe_Open_Output_Raw(const std::string& filename,bool binary,bool write_compressed_if_possible)
{
    bool compressed=File_Is_Compressed(filename);
    std::ios_base::openmode flags=std::ios::out;if(binary) flags|=std::ios::binary;
    if(!write_compressed_if_possible && File_Exists_Ignoring_Compression_Suffix(filename+".gz")){
        LOG::cerr<<"Refusing to write "<<filename<<" uncompressed when compressed version already exists\n";
        throw FILESYSTEM_ERROR("Safe_Open_Output failed (compressed already exists)");}
    if(compressed && !binary){LOG::cerr<<"Refusing to open compressed file "<<filename<<"in text mode\n";
        throw FILESYSTEM_ERROR("Safe_Open_Output failed (tried to make text compressed)");}
    std::string actual_filename;
#ifndef COMPILE_WITHOUT_ZLIB_SUPPORT
    if(binary && write_compressed_if_possible){actual_filename=compressed?filename:filename+".gz";
        if(File_Writable(actual_filename))
            return Gzip_Out(actual_filename,flags);
    }
    else{actual_filename=compressed?Strip_Compression_Suffix(filename):filename;
        std::ostream* output=new std::ofstream(actual_filename.c_str(),flags);if(*output) return output;}
#else
    actual_filename=compressed?Strip_Compression_Suffix(filename):filename;
    std::ostream* output=new std::ofstream(actual_filename.c_str(),flags);if(*output) return output;
#endif
    LOG::cerr<<"Can't open "<<actual_filename<<" for write "<<(binary?"(binary)":"")<<std::endl;
    throw FILESYSTEM_ERROR("Safe_Open_Input failed");
}
//###################################################################
// Function Compare_File_Times
//###################################################################
int Compare_File_Times(const std::string& filename1,const std::string& filename2)
{
    return Compare_File_Times_Ignoring_Compression_Suffix(Real_File(filename1),Real_File(filename2));
}
//###################################################################
// Function File_Exists_Ignoring_Compression_Suffix
//###################################################################
bool File_Exists_Ignoring_Compression_Suffix(const std::string& filename)
{
    return !std::ifstream(filename.c_str()).fail();
}
//###################################################################
// Function File_Writable_Ignoring_Compression_Suffix
//###################################################################
bool File_Writable_Ignoring_Compression_Suffix(const std::string& filename)
{
    return !std::ofstream(filename.c_str(),std::ios::out).fail();
} // TODO: make this not create the file
//###################################################################
// Function Directory_Writable
//###################################################################
bool Directory_Writable(const std::string& dirname) // TODO: make this nicer
{
    static const char* dummy_filename="_PHYSBAM_FILE_UTILITIES_DUMMY_";
    std::string filename=dirname+"/"+dummy_filename;
    bool success=!std::ofstream(filename.c_str()).fail();
    remove(filename.c_str());
    return success;
}
//###################################################################
// Function Find_First_Nonexistent_File_In_Sequence
//###################################################################
std::string Find_First_Nonexistent_File_In_Sequence(std::string filename_pattern,const int id_start,int* id_result)
{
    int id=id_start;
    while(File_Exists(LOG::sprintf(filename_pattern.c_str(),id))) id++;
    if(id_result) *id_result=id;
    return LOG::sprintf(filename_pattern.c_str(),id);
}
//###################################################################
// Function Find_First_Nonexistent_Directory_In_Sequence
//###################################################################
std::string Find_First_Nonexistent_Directory_In_Sequence(std::string directory_pattern,const int id_start,int* id_final)
{
    int id=id_start;
    while(Directory_Exists(LOG::sprintf(directory_pattern.c_str(),id))) id++;
    if(id_final) *id_final=id;
    return LOG::sprintf(directory_pattern.c_str(),id);
}
//###################################################################
// Function Make_First_Nonexistent_Directory_In_Sequence
//###################################################################
std::string Make_First_Nonexistent_Directory_In_Sequence(std::string directory_pattern,const int id_start,int* id_final)
{
    std::string output_directory=Find_First_Nonexistent_Directory_In_Sequence(directory_pattern,id_start,id_final);
    Create_Directory(output_directory);
    return output_directory;
}
//#####################################################################
std::string Get_File_Extension_Ignoring_Compression_Suffix(const std::string &filename)
{
    std::string::size_type lastperiod=filename.rfind('.'),lastslash=filename.rfind('/');
    if(lastperiod!=std::string::npos&&(lastslash==std::string::npos||lastperiod>lastslash))
        return filename.substr(lastperiod+1);
    return "";
}

std::string Get_Basename_Ignoring_Compression_Suffix(const std::string& filename)
{
    std::string::size_type lastperiod=filename.rfind(".");std::string::size_type lastslash=filename.rfind("/");
    if(lastperiod!=std::string::npos&&(lastslash==std::string::npos||lastperiod>lastslash))
        return filename.substr(0,lastperiod);
    return filename;
}

std::string Get_Short_Name_Ignoring_Compression_Suffix(const std::string& filename)
{
    size_t lastslash=filename.rfind('/');
    if(lastslash==std::string::npos) return filename;
    return filename.substr(lastslash+1);
}

bool File_Extension_Matches_Ignoring_Compression_Suffix(const std::string& filename,const char* ext,int len_ext,const bool case_sensitive)
{
    int len=filename.length();
    const char* end=filename.c_str()+len;
    if(len<=len_ext+1) return false;
    if(len>4+len_ext && !Compare_Strings(end-3,3,".gz",3,case_sensitive)) end-=3;
    if(Compare_Strings(end-len_ext,len_ext,ext,len_ext,case_sensitive)) return false;
    return end[-len_ext-1]=='.';
}

bool File_Is_Compressed(const std::string& filename)
{
    int len=filename.length();
    const char* end=filename.c_str()+len;
    return len>4 && !Compare_Strings(end-3,3,".gz",3,true);
}

bool File_Exists(const std::string& filename)
{
    return File_Exists_Ignoring_Compression_Suffix(filename) || (!File_Is_Compressed(filename) && File_Exists_Ignoring_Compression_Suffix(filename+".gz"));
}

bool File_Writable(const std::string& filename)
{
    return File_Writable_Ignoring_Compression_Suffix(filename) || (!File_Is_Compressed(filename) && File_Writable_Ignoring_Compression_Suffix(filename+".gz"));
}

std::string Strip_Compression_Suffix(const std::string& filename)
{
    if(File_Is_Compressed(filename)) return Get_Basename_Ignoring_Compression_Suffix(filename);
    return filename;
}

std::string Real_File(const std::string& filename)
{
    if(File_Exists_Ignoring_Compression_Suffix(filename))return filename;
    else if(!File_Is_Compressed(filename) && File_Exists_Ignoring_Compression_Suffix(filename+".gz"))return filename+".gz";
    return "";
}

std::string Get_File_Extension(const std::string &filename)
{
    return Get_File_Extension_Ignoring_Compression_Suffix(Strip_Compression_Suffix(filename));
}

std::string Get_Basename(const std::string& filename)
{
    return Get_Basename_Ignoring_Compression_Suffix(Strip_Compression_Suffix(filename));
}

std::string Get_Base_Directory_Name(const std::string& path)
{
    size_t lastslash=path.rfind('/'); if(lastslash==std::string::npos) return ".";
    return path.substr(0,lastslash);
}

std::string Get_Short_Name(const std::string& filename)
{
    return Get_Short_Name_Ignoring_Compression_Suffix(Strip_Compression_Suffix(filename));
}

bool Is_Animated(const std::string &filename)
{
    return filename.find("%d") != std::string::npos;
}

std::string Get_Frame_Filename(const std::string &filename,int frame)
{
    return Is_Animated(filename)?LOG::sprintf(filename.c_str(),frame):filename;
}

bool Frame_File_Exists(const std::string &filename,int frame)
{
    return File_Exists(Get_Frame_Filename(filename,frame));
}

void Ignore(std::istream& input,char c)
{
    while(isspace(input.peek())) input.get();
    if(input.peek()==c) input.get();
}
}
