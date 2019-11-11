//#####################################################################
// Copyright 2019, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <fstream>
namespace PhysBAM{

VIEWER_DIR::VIEWER_DIR(const std::string& output_directory)
    :output_directory(output_directory)
{
}

VIEWER_DIR::~VIEWER_DIR()
{
}

void VIEWER_DIR::Advance_Directory(int substep_level)
{
    if(!frame_stack.m)
    {
        assert(substep_level==0);
        frame_stack.Append(0);
    }
    else
    {
        assert((unsigned)substep_level<10);
        frame_stack.Resize(substep_level+1);
        frame_stack.Last()++;
    }
    Update_Current_Directory();
}

void VIEWER_DIR::Update_Current_Directory()
{
    if(!frame_stack.m) current_directory="";
    else{
        current_directory=output_directory+"/"+std::to_string(frame_stack(0));
        for(int i=1;i<frame_stack.m;i++)
            current_directory+="."+std::to_string(frame_stack(i));}
}

void VIEWER_DIR::Make_Common_Directory()
{
    if(made_common) return;
    Create_Directory(output_directory);
    Create_Directory(output_directory+"/common");
    if(!no_log)
        LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    made_common=true;
}

void VIEWER_DIR::Start_Directory(int substep_level,const char* title)
{
    Advance_Directory(substep_level);
    if(First_Frame()) Make_Common_Directory();
    Create_Directory(current_directory);
    if(title) Write_To_Text_File(current_directory+"/frame_title",title);
}

void VIEWER_DIR::Finish_Directory()
{
    std::ofstream lf(output_directory+"/common/last_frame");
    frame_stack.Write_Raw(lf);
    lf<<std::endl;
}

bool VIEWER_DIR::Find_Next_Directory(int substep_level)
{
    ARRAY<int> stack_save=frame_stack;
    std::string dir_save=current_directory;
    assert(substep_level>=0 && substep_level<10);
    int last_level=frame_stack.m-1;
    for(int i=frame_stack.m;i<=substep_level;i++)
    {
        Advance_Directory(i);
        if(Directory_Exists(current_directory))
            return true;
    }
    for(int i=last_level;i>=0;i--)
    {
        Advance_Directory(i);
        if(Directory_Exists(current_directory))
            return true;
    }
    frame_stack.Exchange(stack_save);
    current_directory=std::move(dir_save);
    return false;
}

bool VIEWER_DIR::Find_Prev_Directory(int substep_level)
{
    if(frame_stack.m==1 && frame_stack(0)==0) return false;
    assert(frame_stack.Last()>0);
    frame_stack.Last()--;
    while(frame_stack.m>1 && frame_stack.Last()==0) frame_stack.Pop();
    Update_Current_Directory();
    assert(Directory_Exists(current_directory));
    for(int s=frame_stack.m;s<=substep_level;s++)
    {
        frame_stack.Append(0);
        while(Directory_Exists(current_directory+"."+std::to_string(frame_stack.Last()+1)))
            frame_stack.Last()++;
        current_directory+="."+std::to_string(frame_stack.Last());
    }
    assert(Directory_Exists(current_directory));
    return true;
}

void VIEWER_DIR::Set(int frame)
{
    frame_stack.Resize(1);
    frame_stack(0)=frame;
    current_directory=output_directory+"/"+std::to_string(frame_stack(0));
}

void VIEWER_DIR::Set(const char* frame_string)
{
    Parse_Frame(frame_stack, frame_string);
    Update_Current_Directory();
}

void VIEWER_DIR::Read_Last_Frame(ARRAY<int>& stack) const
{
    std::string line;
    std::ifstream lf(output_directory+"/common/last_frame");
    getline(lf,line);
    Parse_Frame(stack, line);
}

void VIEWER_DIR::Parse_Frame(ARRAY<int>& stack, const char* frame_string) const
{
    stack.Remove_All();
    char* t=0;
    while(1)
    {
        frame_string+=strcspn(frame_string,"0123456789");
        if(!*frame_string) break;
        stack.Append(strtol(frame_string, &t, 10));
        frame_string=t;
    }
}

}
