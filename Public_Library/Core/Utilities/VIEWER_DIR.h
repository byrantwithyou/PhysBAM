//#####################################################################
// Copyright 2019, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VIEWER_DIR
//#####################################################################
#ifndef __VIEWER_DIR__
#define __VIEWER_DIR__

#include <Core/Arrays/ARRAY.h>
#include <string>
namespace PhysBAM{
struct VIEWER_DIR
{
    static const int max_substep_level = 10;
    std::string output_directory;
    std::string current_directory;
    ARRAY<int> frame_stack;
    bool no_log=false;
    bool made_common=false;
    
    explicit VIEWER_DIR(const std::string& output_directory);
    VIEWER_DIR(const VIEWER_DIR&)=delete;
    VIEWER_DIR&operator=(const VIEWER_DIR&)=delete;
    ~VIEWER_DIR();

    // Write
    void Make_Common_Directory(bool append_log=false);
    void Start_Directory(int substep_level,const char* title);
    void Start_Directory(int substep_level,const std::string& title)
    {Start_Directory(substep_level,title.c_str());}
    void Finish_Directory();
    bool First_Frame() const
    {return frame_stack.m==1 && frame_stack(0)==0;}

    // Read
    bool Find_Next_Directory(int substep_level);
    bool Find_Prev_Directory(int substep_level);
    void Read_Last_Frame(ARRAY<int>& stack) const;
    void Read_Last_Frame(int substep_level=max_substep_level);
    void Set(int frame);
    void Set(const std::string& frame_string)
    {Set(frame_string.c_str());}
    void Set(const char* frame_string);

    // Helper utilities
    void Advance_Directory(int substep_level);
    void Update_Current_Directory();
    void Parse_Frame(ARRAY<int>& stack,const std::string& frame_string) const
    {Parse_Frame(stack,frame_string.c_str());}
    void Parse_Frame(ARRAY<int>& stack,const char* frame_string) const;
};
}
#endif
