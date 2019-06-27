//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FILE_UTILITIES__
#define __FILE_UTILITIES__

#include <Core/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Core/Read_Write/STRING_UTILITIES.h>
#include <Core/Utilities/EXCEPTIONS.h>
#include <cstring>
namespace PhysBAM{

//###################################################################
// Platform Specific Function Definitions
//###################################################################
bool Directory_Exists(const std::string& dirname);
bool Directory_Exists(const char* dirname);
bool Create_Directory(char* dirname);
bool Create_Directory(const char* dirname);
bool Create_Directory(const std::string& dirname);
bool Remove_File(const std::string& filename,bool check_compressed=true);
bool Remove_Directory(const std::string& filename);
std::string Real_Path(const std::string& path);
int Compare_File_Times_Ignoring_Compression_Suffix(const std::string& filename1,const std::string& filename2);
std::istream* Safe_Open_Input_Raw(const std::string& filename,bool binary=true);
std::ostream* Safe_Open_Output_Raw(const std::string& filename,bool binary=true,bool write_compressed_if_possible=true);
inline void Safe_Open_Input(FILE_ISTREAM& stream,const std::string& filename,bool binary=true)
{return stream.Set(Safe_Open_Input_Raw(filename,binary));}
inline void Safe_Open_Output(FILE_OSTREAM& stream,STREAM_TYPE stream_type,const std::string& filename,bool binary=true,bool write_compressed_if_possible=true)
{return stream.Set(Safe_Open_Output_Raw(filename,binary,write_compressed_if_possible),stream_type);}
FILE* Temporary_File();
//###################################################################
// Platform Non-Specific Function Declarations
//###################################################################
int Compare_File_Times(const std::string& filename1,const std::string& filename2);
bool File_Exists_Ignoring_Compression_Suffix(const std::string& filename);
bool File_Writable_Ignoring_Compression_Suffix(const std::string& filename);
bool Directory_Writable(const std::string& dirname);
std::string Get_Working_Directory();
//###################################################################
// Platform Non-Specific Function Definitions
//###################################################################
std::string Get_File_Extension_Ignoring_Compression_Suffix(const std::string &filename);
std::string Get_Basename_Ignoring_Compression_Suffix(const std::string& filename);
std::string Get_Short_Name_Ignoring_Compression_Suffix(const std::string& filename);
bool File_Extension_Matches_Ignoring_Compression_Suffix(const std::string& filename,const char* ext,int len_ext,const bool case_sensitive);
inline bool File_Extension_Matches_Ignoring_Compression_Suffix(const std::string& filename,const char* ext,const bool case_sensitive=true)
{return File_Extension_Matches_Ignoring_Compression_Suffix(filename,ext,strlen(ext),case_sensitive);}
inline bool File_Extension_Matches_Ignoring_Compression_Suffix(const std::string& filename,const std::string& ext,const bool case_sensitive=true)
{return File_Extension_Matches_Ignoring_Compression_Suffix(filename,ext.c_str(),ext.length(),case_sensitive);}
bool File_Is_Compressed(const std::string& filename);
bool File_Exists(const std::string& filename);
bool File_Writable(const std::string& filename);
std::string Strip_Compression_Suffix(const std::string& filename);
std::string Real_File(const std::string& filename);

std::string Get_File_Extension(const std::string &filename);
std::string Get_Basename(const std::string& filename);
std::string Get_Base_Directory_Name(const std::string& path);
std::string Get_Short_Name(const std::string& filename);

inline std::string Number_To_String(const int i)
{return Value_To_String(i);}

std::string Find_First_Nonexistent_File_In_Sequence(std::string filename_pattern,const int id_start=0,int* id_result=0);
std::string Find_First_Nonexistent_Directory_In_Sequence(std::string directory_pattern,const int id_start=0,int* id_final=0);
std::string Make_First_Nonexistent_Directory_In_Sequence(std::string directory_pattern,const int id_start=0,int* id_final=0);
//###################################################################
// Utilities for "animated" files (files with %d for frame number)
//###################################################################
bool Is_Animated(const std::string &filename);
std::string Get_Frame_Filename(const std::string &filename,int frame);
bool Frame_File_Exists(const std::string &filename,int frame);
//#####################################################################
// Read_From_File
//#####################################################################
// Convenience functions
template<class T1,class ...Args>
inline void Read_From_File(const std::string& filename,T1& d1,Args&& ...args)
{FILE_ISTREAM input;Safe_Open_Input(input,filename);Read_Binary(input,d1,args...);}
// runtime float/double versions
//#####################################################################
// Write_To_File
//#####################################################################
// Convenience functions
// runtime float/double versions
template<class T1,class ...Args>
inline void Write_To_File(const STREAM_TYPE stream_type,const std::string& filename,const T1& d1,Args&& ...args)
{FILE_OSTREAM output;Safe_Open_Output(output,stream_type,filename);Write_Binary(output,d1,args...);}
template<class RW,class T1,class ...Args>
inline void Write_To_File(const std::string& filename,const T1& d1,Args&& ...args)
{Write_To_File(STREAM_TYPE((RW)0),filename,d1,args...);}
//#####################################################################
// Read_From_Text_File
//#####################################################################
// Convenience function
inline void Read_From_Text_File_Helper(std::istream& in){}
template<class T1,class ...Args>
inline void Read_From_Text_File_Helper(std::istream& in,T1& d1,Args&& ...args)
{in>>d1;Read_From_Text_File_Helper(in,args...);}
template<class ...Args>
inline void Read_From_Text_File(const std::string& filename,Args&& ...args)
{std::istream* input=Safe_Open_Input_Raw(filename,false);Read_From_Text_File_Helper(*input,args...);delete input;}
//#####################################################################
// Write_To_Text_File
//#####################################################################
// Convenience functions
inline void Write_To_Text_File_Helper(std::ostream& out){}
template<class T1,class ...Args>
inline void Write_To_Text_File_Helper(std::ostream& out,T1& d1,Args&& ...args)
{out<<d1;Write_To_Text_File_Helper(out,args...);}
template<class ...Args>
inline void Write_To_Text_File(const std::string& filename,Args&& ...args)
{std::ostream* output=Safe_Open_Output_Raw(filename,false);Write_To_Text_File_Helper(*output,args...);delete output;}
//#####################################################################
// Create_From_File
//#####################################################################
template<class T>
inline void Create_From_File(const std::string& filename,T*& d)
{typename remove_const<T>::type* read=T::Create();Read_From_File(filename,*read);d=read;}

template<class T>
inline T* Create_From_File(const std::string& filename)
{typename remove_const<T>::type* d=T::Create();Read_From_File(filename,*d);return d;}

void Ignore(std::istream& input,char c);
//#####################################################################
}
#endif
