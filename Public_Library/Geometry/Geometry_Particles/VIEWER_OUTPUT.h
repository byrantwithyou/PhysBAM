//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VIEWER_OUTPUT
//#####################################################################
#ifndef __VIEWER_OUTPUT__
#define __VIEWER_OUTPUT__

#include <Core/Utilities/VIEWER_DIR.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

class VIEWER_OUTPUT
{
public:
    VIEWER_DIR& viewer_dir;
    STREAM_TYPE stream_type;

    ARRAY<std::function<void()> > common_entries;
    ARRAY<std::function<void()> > entries;

    VIEWER_OUTPUT(STREAM_TYPE stream_type,VIEWER_DIR& viewer_dir);
    ~VIEWER_OUTPUT();

    static VIEWER_OUTPUT* Singleton(VIEWER_OUTPUT* vo=0);

    void Flush_Frame(const char* title);
    void Flush_Frame(const std::string& title)
    {Flush_Frame(title.c_str());}

    template<class OBJ>
    void Add_Common(const std::string& name,const OBJ& obj)
    {
        common_entries.Append([this,name,&obj]()
            {Write_To_File(stream_type,viewer_dir.output_directory+"/common/"+name,obj);});
    }

    template<class OBJ>
    void Add(const std::string& name,const OBJ& obj)
    {
        entries.Append([this,name,&obj]()
            {Write_To_File(stream_type,viewer_dir.current_directory+"/"+name,obj);});
    }

    template<class TV> void Use_Debug_Particles();
};
template<class TV> void Use_Debug_Particles();
inline void Flush_Frame(const char* title="")
{VIEWER_OUTPUT::Singleton()->Flush_Frame(title);}
}
#endif
